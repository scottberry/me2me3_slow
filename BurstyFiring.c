#include "definitions.h"

int main(int argc, char *argv[]) {
  FILE *fptr, *parFile, *burstFile;
  char *avgfile, fname[256]="", parameterSpace[256]="", burstyRes[256]="";
  chromatin c;
  parameters p;
  gillespie g;
  record r;
  quantification q;
  long i, j, locus;
  double FIRING, P_DEMETHYLATE, P_METHYLATE;
  double K_ON_MIN, K_ON_MAX, K_OFF;
  int p1, p2, p3;

  /* Code timing */
#ifdef __APPLE__
  mach_timebase_info_data_t info;
  uint64_t start, timeElapsed;
  mach_timebase_info(&info);
  start= mach_absolute_time();
#else
  clock_t start, end, elapsed;
  float timeElapsed;
  start = clock();
#endif

  /* ----------------- */
  /* Simulation setup  */
  /* ----------------- */
    
  /* Note: if sampling frequency is too low, data will not be
     collected in the last hour of each cell cycle, when parameter
     values are small. This will lead to program aborting due to
     Gap = NaN. For robustness in parameter searches use
     p.samples = p.maxReact. Also choose p.maxReact according to
     p.cellCycles so that sampling is frequent enough in relation to
     cell cycle. For 50 cell cycles, p.maxReact = 200000 is a good
     choice for a large parameter search. */

  c.sites = 60;
  p.loci = 10;
  p.maxReact = 200000;
  p.samples = 200000; 
  p.sampleFreq = p.maxReact/p.samples;

  /* Set program run parameters */
  p.cellCycles = 20;
  p.cellCycleDuration = 22.0; // (hours)
  p.optimSteps = 2; 

  /* SILAC specific parameters */
  p.silacExperiment = FALSE;
  p.silacLightCycles = 0;
  p.silacHeavyCycles = 0;
  
  /* Set program run type flags */
  p.DNAreplication = FALSE;
  p.resultsLastHourOnly = TRUE;
  p.resultsFinalLocus = TRUE;
  p.checkHistoneTurnover = FALSE;
  p.stochasticAlpha = FALSE;
  p.burstyFiring = TRUE; 
  g.test = FALSE;
  
  /* Parse command line */
  parseCommandLine(argc,argv,&c,&p);

  p.spatialResults = FALSE;

  if (p.stochasticAlpha == TRUE && p.burstyFiring == TRUE) {
    fprintf(stderr,"Error: burstyFiring and stochasticAlpha options are not compatible");
    exit(-1);
  }
  
  /* Seed RNG */
  if (p.randomSeed == TRUE)
    rseed(&p);
  else
    setseed(&p,time(0) + p.seed); 

  /* create base filename from specified run parameters */
  avgfile = parameterDependentBasename(&c,&p);

  /* open results file and write header */
  strcpy(parameterSpace,"ParamOptimRes_\0"); strcat(parameterSpace,avgfile); 
  parFile = fopen(parameterSpace,"w");
  fprintParameterSpaceHeader(parFile);

  if (p.burstyFiring==TRUE) {
    /* open bursty results file and write header */
    strcpy(burstyRes,"BurstyOptimRes_\0"); strcat(burstyRes,avgfile); 
    burstFile = fopen(burstyRes,"w");
    fprintBurstyResultsHeader(burstFile);
  }
  
  if (g.test==TRUE)
    g.test_fptr = fopen("TestGillespieFromMain.txt","w");

  /* allocate memory and initialise gillespie algorithm */
  allocateGillespieMemory(&c,&p,&g,&r);
  initialiseGillespieFunctions(&c,&g);
  if (p.burstyFiring == TRUE)
    initialiseGillespieFunctionsBurstyFiring(&c,&g);

  /* -------------------------- */
  /* Start loop over parameters */
  /* -------------------------- */
  
  for (p1=0;p1<p.optimSteps;p1++) { // 7
    for (p2=0;p2<p.optimSteps;p2++) {
      for (p3=0;p3<30;p3++) {
	  
        K_ON_MAX = pow(10,-0.1*(p1+30));
        K_OFF = pow(10,-0.1*(p2+20));
        P_DEMETHYLATE = pow(10,-0.1*(p3+4));
             
        // P_DEMETHYLATE = 0.004;
        // K_ON_MAX = 0.001;
        // K_OFF = 0.01;
        
        P_METHYLATE = 0.000008;
        K_ON_MIN = K_ON_MAX/40.0;
        // re-scale minimum firing rate to ensure same average levels
        // of transcription in the repressed state
        FIRING = 0.0001 * K_OFF / K_ON_MIN;

        
        // Transcription
        // -------------------------
        p.transcription_demethylate = P_DEMETHYLATE; 
        
        // Methylation/demethylation
        // -------------------------
        /* 5% noise. Represents basal activity of unstimulated PRC2 */
        p.noisy_me0_me1 = 9.0*P_METHYLATE/20.0;
        p.noisy_me1_me2 = 6.0*P_METHYLATE/20.0;
        p.noisy_me2_me3 = P_METHYLATE/20.0;
        
        /* ratio of 9:6:1 in "specificity constant" k_cat/K_M
           \cite{McCabe:2012kk} \cite{Sneeringer:2010dj} */
        p.me0_me1 = 9.0*P_METHYLATE; 
        p.me1_me2 = 6.0*P_METHYLATE; 
        p.me2_me3 = P_METHYLATE;
        
        /* 2 - 2.5 fold lower Kd for K27me3 than K27me2, together with
           4 - 5 fold allosteric activation give a me2/me3 "factor"
           of 10. That is, me3 is 10-times more likely to stimulate a
           methyl addition on a nearby nucleosome.
           \cite{Margueron:2009el} */
        p.me2factor = 0.1; 
        p.me3factor = 1.0;

        /* noisy demethylation independent of transcription */
        p.noisy_demethylate = P_DEMETHYLATE * FIRING * K_ON_MIN / K_OFF;
        
        // Bursty transcription
        // -------------------
        if (p.burstyFiring == TRUE) {
          p.k_onMin = K_ON_MIN;
          p.k_onMax = K_ON_MAX;
          p.k_off = K_OFF;
          p.constFiring = FIRING;
        }
            
        // Reset results to zero for each parameter set
        resetQuantification(&q);

        /* -------------- */
        /* loop over loci */
        /* -------------- */
        
        for (locus=0;locus<p.loci;locus++) {
          // fprintf(stderr,"locus %ld\n",locus);
          if (p.startM == TRUE) {
            initialiseRepressed(&c);
          } else if (p.startU == TRUE) {
            initialiseActive(&c);
          } else { 
            if (locus < floor(p.loci/2))
              initialiseRepressed(&c);
            else
              initialiseActive(&c);
          }

          if (p.burstyFiring == TRUE)
            initialisePromoterOFF(&c);
          
          /* reset counters */
          p.reactCount = 0;
          p.sampleCount = 0;
          p.cellCycleCount = 0;
          
          /* Schedule first instance of the fixed time reactions */
          g.t_nextRep = p.cellCycleDuration*3600;
          p.firingFactor = 1.0;
          
          /* Reaction loop */
          for (i=0;i<p.maxReact && p.cellCycleCount <= p.cellCycles;i++) {
            if (p.reactCount % p.sampleFreq == 0) {
              r.t_out->el[p.sampleCount] = r.t->el[p.reactCount];
              r.t_outLastSample = p.sampleCount;
              for (j=0;j<(c.sites);j++)
                r.K27->el[j][p.sampleCount] = c.K27->el[j];
              if (p.burstyFiring) {
                r.promoterON->el[p.sampleCount] = c.promoterON;
              }
              p.sampleCount++;
            }
            r.tMax = r.t->el[p.reactCount];
            p.reactCount++;
            gillespieStep(&c,&p,&g,&r);
          }
          accumulateQuantification(&c,&p,&r,&q);
        } /* end loop over loci */

        averageQuantification(&c,&p,&r,&q);
        fprintParameterSpaceResults(parFile,&p,&c,&q);

        if (p.burstyFiring == TRUE)
          fprintBurstyResults(burstFile,&p,&c,&q);
      }
    }
  }
  /* end loop over parameters */
  fclose(parFile);
  fclose(burstFile);
  
  /* print final results */
  if (p.resultsFinalLocus == TRUE) {
    fprintResultsFinalLocus(avgfile,&r);
    
    if (p.burstyFiring == TRUE ) {
      strcpy(fname,"promoter_\0"); strcat(fname,avgfile);
      fprint_promoterStatus_nCycles(fname,&r);
    }

    if (p.spatialResults == TRUE) {
      strcpy(fname,"spatial_\0"); strcat(fname,avgfile);
      fptr = fopen(fname,"w");
      fprintSpatialResults(fptr,&r);
      fclose(fptr);
    }
  }
  /* write log file */
  strcpy(fname,"Log_\0"); strcat(fname,avgfile);
  fptr = fopen(fname,"w");
  writelog(fptr,&c,&p,&r);

  /* free memory */
  freeGillespieMemory(&c,&p,&g,&r);
  
#ifdef __APPLE__
  timeElapsed = mach_absolute_time() - start;
  timeElapsed *= info.numer;
  timeElapsed /= info.denom;
  fprintf(fptr,"Simulation time: %f seconds\n", (float)(timeElapsed)/pow(10,9));
  fprintf(stdout,"Simulation time: %f seconds\n", (float)(timeElapsed)/pow(10,9));
#else
  end = clock();
  timeElapsed = (float)(end-start)/CLOCKS_PER_SEC;
  fprintf(fptr,"Simulation time: %f seconds\n", (float)(timeElapsed));
  fprintf(stdout,"Simulation time: %f seconds\n", (float)(timeElapsed));
#endif
  
  fclose(fptr);
  if (g.test==TRUE)
    fclose(g.test_fptr);
  return(1);
}
