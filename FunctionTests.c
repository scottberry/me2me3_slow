#include "definitions.h"

void fprintf_K27(FILE *fptr, chromatin *c) {
  int i;
  for (i=0;i<c->K27->len;i++) {
    fprintf(fptr,"%ld ",c->K27->el[i]);
  }
  fprintf(fptr,"\n");
}

int main(int argc, char *argv[]) {
  FILE *fptr;
  char avgfile[128]="", fname[128]="", buffer[128]="";
  char id[16]="";
  chromatin c;
  parameters p;
  gillespie g;
  record r;
  signed char initial;
  double gap, Mavg, tTot, tTotM, tTotU, tM, tU, lifetime;
  double firstPassage, firstPassageM, firstPassageU, fpU, fpM;
  long i, j, fh, initM, initU, seed;
  double probM, probU, bistability;
  double FIRING, P_DEMETHYLATE, P_METHYLATE;
  logical startM = FALSE, startU = FALSE, randomSeed = TRUE;
  
  c.sites = 8;
  c.controlSites = c.sites; // can be replaced via command line  
  p.maxReact = 10;
  p.loci = 1;
  p.samples = 10; 
  p.sampleFreq = p.maxReact/p.samples;

  p.cellCycles = 6;
  p.cellCycleDuration = 16.0; // (hours)
  p.alpha = 0.0; // can be replaced via command line
  p.beta = 1.0; // can be replaced via command line
  
  /* Parse command line */
  opterr = 0;
  while ((j = getopt (argc, argv, "c:a:i:mu")) != -1)
    switch (j)
      {
      case 'c':
        sprintf(buffer,"%s",optarg);
        c.controlSites = atoi(buffer);
        break;
        
      case 'a':
        sprintf(buffer,"%s",optarg);
        p.alpha = atof(buffer);
        break;

      case 'i':
        randomSeed = FALSE;
        sprintf(id,"%s",optarg);
        seed = atoi(id);
        sprintf(id,"_%s",optarg);
        break;

      case 'm':
        startM = TRUE;
        break;

      case 'u':
        startU = TRUE;
        break;
        
      default:
        usage();
      }
  
  /* Seed RNG */
  if (randomSeed == TRUE)
    rseed(&p);
  else
    setseed(&p,time(0) + seed);

  /* Memory allocation */
  c.K27 = i_vec_get( c.sites );
  g.methylate_index = i_vec_get( c.sites );
  g.transcribeDNA_index = i_vec_get( 1 );
  g.propensity = d_vec_get( c.sites + 1 );
  g.doReaction = malloc(g.propensity->len*sizeof( func_ptr_t ) );
  g.doReactionParam = i_vec_get( g.propensity->len );
  g.update = malloc(sizeof( flags ) );

  r.t = d_vec_get(p.maxReact + 1);
  r.firing = i_vec_get(p.maxReact + 1);
  r.t_out = d_vec_get(p.samples);
  r.K27 = i_mat_get(c.sites,p.samples);
  
  /* Initialisation */
  initialiseGillespieFunctions(&c,&g);

  FIRING = 0.0064;
  P_DEMETHYLATE = 0.1;
  P_METHYLATE = 0.00015;
        
  // Transcription
  // ------------------------------------------------------------
  p.firingThreshold = 1.0;
  p.firingRateMin = 0.0004; 
  p.firingRateMax = FIRING; // Optimise
  p.firingCap = 0.0166667; // Cap firing rate at ~ every minute.
  p.transcription_demethylate = P_DEMETHYLATE;
  p.transcription_turnover = 0.0;
  if (p.firingRateMax < p.firingRateMin) {
    fprintf(stderr,"Error: Max firing rate less than min firing \
rate. Setting k_min = k_max\n");
    p.firingRateMin = p.firingRateMax;
  }
  
  // Methylation/demethylation
  // ------------------------------------------------------------
  /* 5% noise. Represents basal activity of unstimulated PRC2 */
  p.noisy_me0_me1 = 9*P_METHYLATE/20.0;
  p.noisy_me1_me2 = 6*P_METHYLATE/20.0;
  p.noisy_me2_me3 = P_METHYLATE/20.0;
  
  /* ratio of 9:6:1 in "specificity constant" k_cat/K_M
     \cite{McCabe:2012kk} \cite{Sneeringer:2010dj} */
  p.me0_me1 = 9*P_METHYLATE; 
  p.me1_me2 = 6*P_METHYLATE; 
  p.me2_me3 = P_METHYLATE;
  
  /* 2 - 2.5 fold lower Kd for K27me3 than K27me2, together with
     4 - 5 fold allosteric activation give a me2/me3 "factor"
     of 10. That is, me3 is 10-times more likely to stimulate a
     methyl addition on a nearby nucleosome.
     \cite{Margueron:2009el} */
  p.me2factor = 0.1; 
  p.me3factor = 1;
  
  // Set results to zero for accumulation over each parameter set
  // ------------------------------------------------------------        
  gap = 0.0;
  Mavg = 0.0;
  probM = 0.0;
  probU = 0.0;
  fh = 0;
  tTot = tTotM = tTotU = 0.0;
  initM = initU = 0;
  firstPassageM = firstPassageU = 0.0;


  // Function tests
  // ------------------------------------------------------------        

  // runif
  setseed(&p,0);
  fptr = fopen("runif_SetSeed_TestRes.txt","w");
  for (i=0;i<1000;i++) {
    fprintf(fptr,"%0.10f\n",runif(p.gsl_r));
  }
  fclose(fptr);

  rseed(&p);
  fptr = fopen("runif_RandomSeed_TestRes.txt","w");
  for (i=0;i<1000;i++) {
    fprintf(fptr,"%0.10f\n",runif(p.gsl_r));
  }
  fclose(fptr);

  // methylate
  fptr = fopen("modifications_TestRes.txt","w");
  fprintf(fptr,"Initialise active:\n");
  initialiseActive(&c);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"\nMethylate position 0:\n");
  methylate(&c,&p,g.update,0);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"Methylate position 6:\n");
  methylate(&c,&p,g.update,6);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"Methylate position 6 again:\n");
  methylate(&c,&p,g.update,6);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"Methylate last position:\n");
  methylate(&c,&p,g.update,c.sites-1);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"\nDemethylate position 0:\n");
  demethylate(&c,&p,g.update,0);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"Demethylate position 6:\n");
  demethylate(&c,&p,g.update,6);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"Demethylate position 6 again:\n");
  demethylate(&c,&p,g.update,6);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"Demethylate last position:\n");
  demethylate(&c,&p,g.update,c.sites-1);
  fprintf_K27(fptr,&c);

  fprintf(fptr,"\nInitialise repressed:\n");
  initialiseRepressed(&c);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"\nReplicate DNA:\n");
  replicateDNA(&c,&p,g.update);
  fprintf_K27(fptr,&c);
  fprintf(fptr,"\nTranscribe (10x):\n");
  for (i=0;i<20;i++) {
    transcribeDNA(&c,&p,g.update,0);
    fprintf_K27(fptr,&c);    
  }
  fclose(fptr);  

  g.test = TRUE;
  g.test_fptr = fopen("gillespie_TestRes.txt","w");
  initialiseGillespieFunctions(&c,&g);

  fprintf(g.test_fptr,"\nChromatin state:\n");    
  initialiseRepressed(&c);
  for (i=0;i<10;i++) {
    transcribeDNA(&c,&p,g.update,0);
  }
  fprintf_K27(g.test_fptr,&c);

  p.firingThreshold = 1.0;
  p.firingRateMin = 0.0004; 
  p.firingRateMax = 0.04;
  p.firingFactor = 1.0;
  p.alpha = 0.0;
  
  fprintf(g.test_fptr,"\nFiringRate:\n");
  fprintf(g.test_fptr,"f_me2_me3 = 1.0, firingRate = %0.4f\n",firingRate(&p,1.0));
  fprintf(g.test_fptr,"f_me2_me3 = 0.5, firingRate = %0.4f\n",firingRate(&p,0.5));
  fprintf(g.test_fptr,"f_me2_me3 = 0.0, firingRate = %0.4f\n",firingRate(&p,0.0));
  
  c.controlSites = c.sites;
  fprintf(g.test_fptr,"\nControl Sites: %ld\n",c.controlSites);
  fprintf(g.test_fptr,"fracControlRegion_me2me3 = %0.3f\n",fracControlRegion_me2me3(&c));
  c.controlSites = 4;
  fprintf(g.test_fptr,"Control Sites: %ld\n",c.controlSites);
  fprintf(g.test_fptr,"fracControlRegion_me2me3 = %0.3f\n",fracControlRegion_me2me3(&c));

  fprintf(g.test_fptr,"\nChromatin state:\n");    
  initialiseRepressed(&c);
  for (i=0;i<10;i++) {
    transcribeDNA(&c,&p,g.update,0);
  }
  fprintf_K27(g.test_fptr,&c);

  fprintf(g.test_fptr,"\nNeighbours K27 factor:\n");
  for (i=0;i<c.sites;i++)
    fprintf(g.test_fptr,"pos = %ld, factor = %0.1f\n",i,neighboursK27factor(&c,&p,i));

  fprintf(g.test_fptr,"\nPropensities:\n");
  g.update->histone=TRUE;
  p.firingFactor = 1.0;
  updatePropensities(&c,&p,&g);
  for (i=0;i<g.propensity->len;i++)
    fprintf(g.test_fptr,"i = %ld , propensity = %0.10f\n",i,g.propensity->el[i]);
  
  fclose(g.test_fptr);  
  
  g.test = FALSE;
  // reset counters
  p.reactCount = 0;
  p.sampleCount = 0;
  p.cellCycleCount = 0;
  
  // Schedule first instance of the fixed time reactions
  g.t_nextRep = p.cellCycleDuration*3600;
  p.firingFactor = 1.0;
  
  /* Reaction loop */
  for (i=0;i<p.maxReact && p.cellCycleCount <= p.cellCycles;i++) {
    if (p.reactCount % p.sampleFreq == 0) {
      for (j=0;j<(c.sites);j++) {
        r.t_out->el[p.sampleCount] = r.t->el[p.reactCount];
        r.K27->el[j][p.sampleCount] = c.K27->el[j];
        r.t_outLastSample = p.sampleCount;
      }
      p.sampleCount++;
    }
    r.tMax = r.t->el[p.reactCount];
    p.reactCount++;
    gillespieStep(&c,&p,&g,&r);
  }
  
  if (p.resultsLastHourOnly == TRUE) {
    gap += tAverageGap_lastHour_nCycles(&c,&p,&r);
    Mavg += tAverage_me2_me3_lastHour_nCycles(&c,&p,&r);
    probM += prob_me2_me3_lastHour_nCycles(&c,&p,&r);
    probU += prob_me0_me1_lastHour_nCycles(&c,&p,&r);
  } else {
    gap += tAverageGap_nCycles(&c,&p,&r);
    Mavg += tAverage_me2_me3_nCycles(&c,&p,&r);
    probM += prob_me2_me3_nCycles(&c,&p,&r);
    probU += prob_me0_me1_nCycles(&c,&p,&r);
  }
  fh += numberHistoneStateFlips(&r);
  tTot += r.t->el[p.reactCount];
  
  
  firstPassage = firstPassageTime(&r,&initial);
  if (initial==-1) {
    firstPassageM += firstPassage;
    initM++;
    tTotM += r.t->el[p.reactCount];
  } else {
    firstPassageU += firstPassage;
    initU++;
    tTotU += r.t->el[p.reactCount];
  }
  
  if (fh != 0) lifetime = tTot/fh;
  else lifetime = -1.0;
  
  bistability = 4*probM*probU/(p.loci*p.loci);
  
  if (initM != 0) {
    fpM = firstPassageM/initM;
    tM = tTotM/initM;
  } else {
    fpM = -1.0;
    tM = -1.0;
  }
  
  if (initU != 0) {
    fpU= firstPassageU/initU;
    tU = tTotU/initU;
  } else {
    fpU = -1.0;
    tU = -1.0;
  }
  
  /* ---------------------------------------------------------------------------- */
  /* Tidy up and write results files */
  /* ---------------------------------------------------------------------------- */

  /* free all arrays */
  i_vec_free(c.K27);
  i_vec_free(g.methylate_index);
  i_vec_free(g.transcribeDNA_index);
  d_vec_free(g.propensity);
  free(g.doReaction);
  i_vec_free(g.doReactionParam);
  free(g.update);
  rfree(&p);

  if (p.resultsFinalLocus == TRUE) {
    /* print results for final locus */
    strcpy(fname,"t_\0"); strcat(fname,avgfile);
    fprint_t_out_nCycles(fname,&r);
    strcpy(fname,"me0_t_\0"); strcat(fname,avgfile);
    fprint_t_nCycles(fname,r.K27,me0,&r);
    strcpy(fname,"me1_t_\0"); strcat(fname,avgfile);
    fprint_t_nCycles(fname,r.K27,me1,&r);
    strcpy(fname,"me2_t_\0"); strcat(fname,avgfile);
    fprint_t_nCycles(fname,r.K27,me2,&r);
    strcpy(fname,"me3_t_\0"); strcat(fname,avgfile);
    fprint_t_nCycles(fname,r.K27,me3,&r);
    strcpy(fname,"Firing_t_\0"); strcat(fname,avgfile);
    fprint_firing_t_nCycles(fname,&r);
  }
  
  fptr = fopen("Log_TestRes.txt","w");
  writelog(fptr,&c,&p,&r);
  fclose(fptr);
  
  i_vec_free(r.firing);
  i_mat_free(r.K27);
  d_vec_free(r.t);
  d_vec_free(r.t_out);
  
  return(1);
}
