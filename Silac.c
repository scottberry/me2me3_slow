#include "definitions.h"

int main(int argc, char *argv[]) {
  FILE *fptr, *parFile, *silacAbsFile, *silacRelFile, *silacRelAvgFile;
  char *avgfile, fname[256]="", parameterSpace[256]="";
  char silacAbs[256]="", silacRel[256]="", silacRelAvg[256]="";
  chromatin c;
  parameters p;
  gillespie g;
  record r;
  quantification q;
  long i, j, locus;
  double FIRING, P_DEMETHYLATE, P_METHYLATE;
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

  /* -------------------------------------------------------------------------------- */
  /* Simulation setup  */
  /* -------------------------------------------------------------------------------- */
  
  /* Note: if sampling frequency is too low, data will not be
     collected in the last hour of each cell cycle, when parameter
     values are small. This will lead to program aborting due to
     Gap = NaN. For robustness in parameter searches use
     p.samples = p.maxReact. Also choose p.maxReact according to
     p.cellCycles so that sampling is frequent enough in relation to
     cell cycle. For 50 cell cycles, p.maxReact = 200000 is a good
     choice for a large parameter search. */

  c.sites = 60;
  p.loci = 20;
  p.maxReact = 30000;
  p.samples = 30000; 
  p.sampleFreq = p.maxReact/p.samples;

  /* Set program run parameters */
  p.cellCycles = 20;
  p.cellCycleDuration = 22.0; // (hours)
  p.optimSteps = 20;

  /* SILAC specific parameters */
  p.silacExperiment = TRUE;
  p.silacLightCycles = 5;
  p.silacHeavyCycles = 1;

  /* Set program run type flags */
  p.DNAreplication = FALSE;
  p.resultsLastHourOnly = TRUE;
  p.resultsFinalLocus = FALSE;
  p.checkHistoneTurnover = FALSE;
  p.resultsTranscribing = FALSE;
  g.test = FALSE;
  
  /* Parse command line */
  parseCommandLine(argc,argv,&c,&p);
    
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

  if (g.test==TRUE)
    g.test_fptr = fopen("TestGillespieFromMain.txt","w");
  
  /* allocate memory and initialise gillespie algorithm */
  allocateGillespieMemory(&c,&p,&g,&r);
  initialiseGillespieFunctions(&c,&g);
  
  /* allocate memory for SILAC results */  
  if (p.silacExperiment == TRUE) {
    allocateSilacRecordMemory(&c,&p,&r);

    strcpy(silacRelAvg,"SilacRelAverage_\0"); strcat(silacRelAvg,avgfile); 
    silacRelAvgFile = fopen(silacRelAvg,"w");
    fprintf(silacRelAvgFile,"FIRING\tP_DEMETHYLATE\tP_METHYLATE\tFIRING_THRESHOLD\ttime\tmod\tlabel\tlevel\n");
    
    if (p.resultsSilacEachLocus == TRUE) {
      strcpy(silacAbs,"SilacAbs_\0"); strcat(silacAbs,avgfile); 
      silacAbsFile = fopen(silacAbs,"w");
      fprintf(silacAbsFile,"FIRING\tP_DEMETHYLATE\tP_METHYLATE\tFIRING_THRESHOLD\tlocus\ttime\tmod\tlabel\tlevel\n");
      strcpy(silacRel,"SilacRel_\0"); strcat(silacRel,avgfile); 
      silacRelFile = fopen(silacRel,"w");
      fprintf(silacRelFile,"FIRING\tP_DEMETHYLATE\tP_METHYLATE\tFIRING_THRESHOLD\tlocus\ttime\tmod\tlabel\tlevel\n");
    }    
  }

  /* -------------------------------------------------------------------------------- */
  /* Start loop over parameters */
  /* -------------------------------------------------------------------------------- */
  for (p1=0;p1<1;p1++) { // 7
    for (p2=7;p2<p.optimSteps;p2++) { // min 7
      for (p3=12;p3<p.optimSteps;p3++) { // min 12
	  
        // setseed(&p,p.seed);
                      
        FIRING = 0.0001*40;
        P_DEMETHYLATE = pow(10,-0.15*(p2+4));
        P_METHYLATE = pow(10,-0.13*(p3+24));
        
        // FIRING = 0.0256;
        // P_DEMETHYLATE = 0.02;
        // P_METHYLATE = 0.00003;

        // Transcription
        // -------------
        /* Leave the repressed firing rate fixed at ~ every 2.8 hours. */
        p.firingRateMin = 0.0001; 
        p.firingRateMax = FIRING; // Optimise

        /* Cap firing rate at ~ every minute. */
        p.firingCap = 0.0166667;
        p.transcription_demethylate = P_DEMETHYLATE; 
        // p.transcription_turnover = 0.0; 
        p.transcriptionDelay = 0.0;
        
        if (p.firingRateMax < p.firingRateMin) {
          fprintf(stderr,"Error: Max firing rate less than min firing rate.");
          fprintf(stderr," Setting k_min = k_max\n");
          p.firingRateMin = p.firingRateMax;
        }
        
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
        p.noisy_demethylate = p.noisy_me2_me3;

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

            // Silac report points
            if (r.t->el[p.reactCount] >= p.SILAC_nextReport && p.SILAC_report <= 4) {
              storeTripleSILAC_me3(locus,&p,&r);
              if (p.resultsSilacEachLocus == TRUE)
                fprintTripleSILAC_eachLocus(silacAbsFile,silacRelFile,locus,&p,&r);
              incrementSilacReportPoint(&p);
            }

            p.reactCount++;
            gillespieStep(&c,&p,&g,&r);
          }

          accumulateQuantification(&c,&p,&r,&q);
        } /* end loop over loci */

        averageQuantification(&c,&p,&r,&q);
        fprintParameterSpaceResults(parFile,&p,&c,&q);
        fprintTripleSILAC_average(silacRelAvgFile,&p,&r);
      }
    }
  }
 
  /* end loop over parameters */
  fclose(parFile);
  if (p.silacExperiment == TRUE) {
    fclose(silacRelAvgFile);
    if (p.resultsSilacEachLocus == TRUE) {
      fclose(silacAbsFile);
      fclose(silacRelFile);
    }
  }

  /* print final results */
  if (p.resultsFinalLocus == TRUE) {
    fprintResultsFinalLocus(avgfile,&r);
    if (p.silacExperiment == TRUE) {
      fprintSilacResultsFinalLocus(avgfile,&r);
    }
  }

  /* write log file */
  strcpy(fname,"Log_\0"); strcat(fname,avgfile);
  fptr = fopen(fname,"w");
  writelog(fptr,&c,&p,&r);

  /* free memory */
  freeGillespieMemory(&c,&p,&g,&r);
  if (p.silacExperiment == TRUE) {
    freeSilacRecordMemory(&r);
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
