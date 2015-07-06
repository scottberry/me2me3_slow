#include "definitions.h"

void usage(void)
{
	printf("Usage:\n");
	printf(" -c<control region>\n");
	printf(" -a<gene activation>\n");
        printf(" -i<identifier>\n");
        printf(" -m\n");
        printf(" -u\n");
	exit (8);
}

int main(int argc, char *argv[]) {
  FILE *fptr, *parFile, *silacAbsFile, *silacRelFile, *silacRelAvgFile;
  char avgfile[256]="", fname[256]="", tmp[256]="", buffer[256]="";
  char parameterSpace[256]="", ptmp[256]="", id[16]="";
  char silacAbs[256]="", silacRel[256]="", silacRelAvg[256]="";
  char *decimal = ".", *underscore = "_";
  chromatin c;
  parameters p;
  gillespie g;
  record r;
  signed char initial;
  double gap, Mavg, tTot, tTotM, tTotU, tM, tU, lifetime, me3_end;
  double firstPassage, firstPassageM, firstPassageU, fpU, fpM;
  long i, j, locus, fh, initM, initU, seed;
  double probM, probU, bistability;
  double FIRING, P_DEMETHYLATE, P_METHYLATE;
  int p1, p2, p3;
  logical startM = FALSE, startU = FALSE, randomSeed = TRUE;

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

  /* -------------------------------------------------------------------------------- */
  /* Simulation setup  */
  /* -------------------------------------------------------------------------------- */
  
  c.sites = 60;
  c.controlSites = c.sites; // can be replaced via command line
  
  /* Note: if sampling frequency is too low, data will not be
     collected in the last hour of each cell cycle, when parameter
     values are small. This will lead to program aborting due to
     Gap = NaN. For robustness in parameter searches use
     p.samples = p.maxReact. Also choose p.maxReact according to
     p.cellCycles so that sampling is frequent enough in relation to
     cell cycle. For 50 cell cycles, p.maxReact = 100000 is a good
     choice for a large parameter search. */
  
  p.loci = 1;
  p.maxReact = 30000;
  p.samples = 30000; 
  p.sampleFreq = p.maxReact/p.samples;

  p.cellCycles = 20;
  p.silacLightCycles = 5;
  p.silacHeavyCycles = 1;
  p.cellCycleDuration = 22.0; // (hours)
  p.G2duration = 4.0;
  p.activation = 1.0;
  p.firingThreshold = 1.0; // can be replaced via command line

  // fold-change from non-transcribing
  p.PRC2inhibition = 1.0; // can be replaced via command line
  
  p.DNAreplication = TRUE;
  p.resultsLastHourOnly = TRUE;
  p.silacExperiment = TRUE;
  p.resultsFinalLocus = TRUE;
  p.resultsSilacEachLocus = FALSE;

  p.optimSteps = 1; 
  
  // Test gillespie algorithm
  g.test = FALSE;
  
  /* Parse command line */
  opterr = 0;
  while ((j = getopt (argc, argv, "i:mut:")) != -1)
    switch (j)
      {
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

      case 't':
        sprintf(buffer,"%s",optarg);
        p.firingThreshold = atof(buffer);
        break;
        
      default:
        usage();
      }
  
  /* Seed RNG */
  if (randomSeed == TRUE)
    rseed(&p);
  else
    setseed(&p,time(0) + seed);

  /* Handle filename using command line args */
  sprintf(tmp,"s%ld",c.sites); strcat(avgfile,tmp); 
  sprintf(tmp,"ctrl%ld",c.controlSites); strcat(avgfile,tmp);
  sprintf(tmp,"cc%d",p.cellCycles); strcat(avgfile,tmp);
  sprintf(tmp,"%0.2f",p.activation);
  sprintf(ptmp,"a%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  sprintf(tmp,"%0.2f",p.firingThreshold);
  sprintf(ptmp,"fir%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  sprintf(tmp,"%0.2f",p.G2duration);
  sprintf(ptmp,"tau%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  sprintf(tmp,"st%ld",p.optimSteps); strcat(avgfile,tmp); 
  strcat(avgfile,id);
  strcat(avgfile,".txt\0");

  strcpy(parameterSpace,"ParamOptimRes_\0"); strcat(parameterSpace,avgfile); 

  parFile = fopen(parameterSpace,"w");
  fprintf(parFile,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\
\tFIRING_THRESHOLD\tP_DEMETHYLATE\tP_METHYLATE\tcontrolSites\tactivation\tgap\tMavg       \
\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU        \
\tavgInitU\ttTot\tprobM\tprobU\tbistability\tme3_end\n");
  fprintf(stderr,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\
\tFIRING_THRESHOLD\tP_DEMETHYLATE\tP_METHYLATE\tcontrolSites\tactivation\tgap\tMavg       \
\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU        \
\tavgInitU\ttTot\tprobM\tprobU\tbistability\tme3_end\n");

  if (g.test==TRUE)
    g.test_fptr = fopen("TestGillespieFromMain.txt","w");
  
  /* Memory allocation */
  c.K27 = i_vec_get( c.sites );
  c.silac = i_vec_get( c.sites );
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

  /* Silac initiation */
  if (p.silacExperiment == TRUE) {
    r.silac = i_mat_get(c.sites,p.samples);

    p.SILAC_0h = (double)3600*p.cellCycleDuration*(p.silacLightCycles+1);
    p.SILAC_10h = p.SILAC_0h + (double)3600*10;
    p.SILAC_24h = p.SILAC_0h + (double)3600*24;
    p.SILAC_48h = p.SILAC_0h + (double)3600*48;

    r.silacResultsLight_0h = d_vec_get(p.loci);
    r.silacResultsLight_10h = d_vec_get(p.loci);
    r.silacResultsLight_24h = d_vec_get(p.loci);
    r.silacResultsLight_48h = d_vec_get(p.loci);
    r.silacResultsHeavy_0h = d_vec_get(p.loci);
    r.silacResultsHeavy_10h = d_vec_get(p.loci);
    r.silacResultsHeavy_24h = d_vec_get(p.loci);
    r.silacResultsHeavy_48h = d_vec_get(p.loci);

    if (p.resultsSilacEachLocus == TRUE) {
      strcpy(silacAbs,"SilacAbs_\0"); strcat(silacAbs,avgfile); 
      silacAbsFile = fopen(silacAbs,"w");
      fprintf(silacAbsFile,"FIRING\tP_DEMETHYLATE\tP_METHYLATE\tFIRING_THRESHOLD\tlocus\ttime\tmod\tlabel\tlevel\n");
      strcpy(silacRel,"SilacRel_\0"); strcat(silacRel,avgfile); 
      silacRelFile = fopen(silacRel,"w");
      fprintf(silacRelFile,"FIRING\tP_DEMETHYLATE\tP_METHYLATE\tFIRING_THRESHOLD\tlocus\ttime\tmod\tlabel\tlevel\n");
    }
    
    strcpy(silacRelAvg,"SilacRelAverage_\0"); strcat(silacRelAvg,avgfile); 
    silacRelAvgFile = fopen(silacRelAvg,"w");
    fprintf(silacRelAvgFile,"FIRING\tP_DEMETHYLATE\tP_METHYLATE\tFIRING_THRESHOLD\ttime\tmod\tlabel\tlevel\n");
  }

  /* -------------------------------------------------------------------------------- */
  /* Start loop over parameters */
  /* -------------------------------------------------------------------------------- */
  for (p1=0;p1<1;p1++) { // 7
    for (p2=p.optimSteps-1;p2<p.optimSteps;p2++) {
      for (p3=p.optimSteps-1;p3<p.optimSteps;p3++) {
	  
        // !!! Set seed for debugging - remove for simulations
        // setseed(&p,0);
                      
        FIRING = 0.000277778*20;
        // P_DEMETHYLATE = pow(10,-0.15*(p2+4));
        // P_METHYLATE = pow(10,-0.12*(p3+26));
        
        // FIRING = 0.0256;
        P_DEMETHYLATE = 0.004;
        P_METHYLATE = 0.000008;
        
        // Transcription
        // ------------------------------------------------------------
        p.firingRateMin = 0.000277778; // Leave the repressed firing rate fixed at ~ every 60 min.
        p.firingRateMax = FIRING; // Optimise
        p.transcription_demethylate = P_DEMETHYLATE; // (rate per site per transcription event)
        if (p.firingRateMax < p.firingRateMin) {
          fprintf(stderr,"Error: Max firing rate less than min firing rate. Setting k_min = k_max\n");
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
        me3_end = 0.0;
        probM = 0.0;
        probU = 0.0;
        fh = 0;
        tTot = tTotM = tTotU = 0.0;
        initM = initU = 0;
        firstPassageM = firstPassageU = 0.0;
        
        /* -------------------------------------------------------------------------------- */
        /* loop over loci */
        /* -------------------------------------------------------------------------------- */
        
        for (locus=0;locus<p.loci;locus++) {
          // fprintf(stderr,"locus %ld\n",locus);
          if (argc > 1) {
            if (startM == TRUE) {
              initialiseRepressed(&c);
	      } else if (startU == TRUE) {
              initialiseActive(&c);
            } else { 
              if (locus < floor(p.loci/2))
                initialiseRepressed(&c);
              else
                initialiseActive(&c);
            }
          } else { 
            if (locus < floor(p.loci/2))
              initialiseRepressed(&c);
            else
              initialiseActive(&c);
          }

          // reset counters
          p.reactCount = 0;
          p.sampleCount = 0;
          p.cellCycleCount = 0;

          // reset silac parameters
          p.silacLabel = LIGHT;
          initialiseSilacLight(&c);
          p.SILAC_nextReport = p.SILAC_0h;
          p.SILAC_report = 1;
          
          // Schedule first instance of the fixed time reactions
          g.t_nextRep = p.cellCycleDuration*3600;
          g.t_nextEndG2 = (p.cellCycleDuration + p.G2duration)*3600;
          p.firingFactor = 1.0;
          
          /* Reaction loop */
          for (i=0;i<p.maxReact && p.cellCycleCount <= p.cellCycles;i++) {
            if (p.reactCount % p.sampleFreq == 0) {
              for (j=0;j<(c.sites);j++) {
                r.t_out->el[p.sampleCount] = r.t->el[p.reactCount];
                r.K27->el[j][p.sampleCount] = c.K27->el[j];
                r.silac->el[j][p.sampleCount] = c.silac->el[j];
                r.t_outLastSample = p.sampleCount;
              }
              p.sampleCount++;
            }
            r.tMax = r.t->el[p.reactCount];

            // SILAC report points
            if (r.t->el[p.reactCount] >= p.SILAC_nextReport && p.SILAC_report <= 4) {
              storeTripleSILAC_me3(locus,&p,&r);

              if (p.resultsSilacEachLocus == TRUE)
                fprintTripleSILAC_eachLocus(silacAbsFile,silacRelFile,locus,&p,&r);
              
              if (p.SILAC_report == 1) {
                p.SILAC_report = 2;
                p.SILAC_nextReport = p.SILAC_10h;
              } else if (p.SILAC_report == 2) {
                p.SILAC_report = 3;
                p.SILAC_nextReport = p.SILAC_24h;
              } else if (p.SILAC_report == 3) {
                p.SILAC_nextReport = p.SILAC_48h;
                p.SILAC_report = 4;
              } else if (p.SILAC_report == 4) {
                p.SILAC_report = 5;
              }
            }

            p.reactCount++;
            if (p.cellCycleCount >= p.silacLightCycles)
              p.silacLabel = HEAVY;
            if (p.cellCycleCount >= p.silacLightCycles + p.silacHeavyCycles)
              p.silacLabel = UNLABELLED;

            
            gillespieStep(&c,&p,&g,&r);
          }

          if (p.resultsLastHourOnly == TRUE) {
            gap += tAverageGap_lastHour_nCycles(&c,&p,&r);
            Mavg += tAverage_me2_me3_lastHour_nCycles(&c,&p,&r);
            me3_end += tAverage_me3_lastHour_nCycles(&c,&p,&r);
            probM += prob_lowExpression_lastHour_nCycles(&c,&p,&r);
            probU += prob_highExpression_lastHour_nCycles(&c,&p,&r);
          } else {
            gap += tAverageGap_nCycles(&c,&p,&r);
            Mavg += tAverage_me2_me3_nCycles(&c,&p,&r);
            probM += prob_me2_me3_nCycles(&c,&p,&r);
            probU += prob_me0_me1_nCycles(&c,&p,&r);
          }
          fh += numberHistoneStateFlips(&r);
          tTot += r.t->el[p.reactCount];

          if (isnan(gap)) {
            fprintf(stderr,"Error: gap is nan. Locus %ld\n",locus);
            fprintf(stderr,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\tP_DEMETHYLATE\tP_METHYLATE\n");
            fprintf(stderr,"%0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f\n",
                    p.me0_me1,p.me1_me2,p.me2_me3,p.me2factor,p.me3factor,FIRING,P_DEMETHYLATE,P_METHYLATE);
            exit(-1);
          }

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
        } /* end loop over loci */

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

        fprintf(parFile,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\
\t%0.10f\t%ld\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f\
\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.6f\n",
                p.me0_me1,p.me1_me2,p.me2_me3,p.me2factor,p.me3factor,
                FIRING,p.firingThreshold,
                P_DEMETHYLATE,P_METHYLATE,c.controlSites,p.activation,
                gap/p.loci,Mavg/p.loci,lifetime,initM,fpM,tM,initU,fpU,tU,tTot/p.loci,
                probM/p.loci,probU/p.loci,bistability,me3_end/p.loci);
        fprintf(stderr,"%0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  \
%0.10f  %ld  %0.4f  %0.4f  %0.4f  %0.4f  %ld  %0.4f  %0.4f  %ld  %0.4f  %0.4f \
%0.4f  %0.4f  %0.4f  %0.4f  %0.6f\n",
                p.me0_me1,p.me1_me2,p.me2_me3,p.me2factor,p.me3factor,
                FIRING,p.firingThreshold,
                P_DEMETHYLATE,P_METHYLATE,c.controlSites,p.activation,
                gap/p.loci,Mavg/p.loci,lifetime,initM,fpM,tM,initU,fpU,tU,tTot/p.loci,
                probM/p.loci,probU/p.loci,bistability,me3_end/p.loci);

        fprintTripleSILAC_average(silacRelAvgFile,&p,&r);
      }
    }
  }
 
  /* end loop over parameters */
  fclose(parFile);
  if (p.silacExperiment == TRUE) {
    if (p.resultsSilacEachLocus == TRUE) {
      fclose(silacAbsFile);
      fclose(silacRelFile);
    }
    fclose(silacRelAvgFile);
  }

  /* -------------------------------------------------------------------------------- */
  /* Tidy up and write results files */
  /* -------------------------------------------------------------------------------- */

  /* free all arrays */
  i_vec_free(c.K27);
  i_vec_free(c.silac);
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

    // light histones
    strcpy(fname,"LIGHT_me0_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me0,r.silac,LIGHT,&r);
    strcpy(fname,"LIGHT_me1_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me1,r.silac,LIGHT,&r);
    strcpy(fname,"LIGHT_me2_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me2,r.silac,LIGHT,&r);
    strcpy(fname,"LIGHT_me3_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me3,r.silac,LIGHT,&r);
    
    // heavy histones
    strcpy(fname,"HEAVY_me0_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me0,r.silac,HEAVY,&r);
    strcpy(fname,"HEAVY_me1_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me1,r.silac,HEAVY,&r);
    strcpy(fname,"HEAVY_me2_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me2,r.silac,HEAVY,&r);
    strcpy(fname,"HEAVY_me3_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me3,r.silac,HEAVY,&r);

    // unlabelled histones
    strcpy(fname,"UNLABELLED_me0_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me0,r.silac,UNLABELLED,&r);
    strcpy(fname,"UNLABELLED_me1_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me1,r.silac,UNLABELLED,&r);
    strcpy(fname,"UNLABELLED_me2_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me2,r.silac,UNLABELLED,&r);
    strcpy(fname,"UNLABELLED_me3_t_\0"); strcat(fname,avgfile);
    fprint_silac_t_nCycles(fname,r.K27,me3,r.silac,UNLABELLED,&r);    
  }
  
  strcpy(fname,"Log_\0"); strcat(fname,avgfile);
  fptr = fopen(fname,"w");
  writelog(fptr,&c,&p,&r);

  i_vec_free(r.firing);
  i_mat_free(r.K27);
  d_vec_free(r.t);
  d_vec_free(r.t_out);

  if (p.silacExperiment == TRUE) {
    i_mat_free(r.silac);
    d_vec_free(r.silacResultsLight_0h);
    d_vec_free(r.silacResultsLight_10h);
    d_vec_free(r.silacResultsLight_24h);
    d_vec_free(r.silacResultsLight_48h);
    d_vec_free(r.silacResultsHeavy_0h);
    d_vec_free(r.silacResultsHeavy_10h);
    d_vec_free(r.silacResultsHeavy_24h);
    d_vec_free(r.silacResultsHeavy_48h);
  }
  
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
