#include "definitions.h"

int main(int argc, char *argv[]) {
  FILE *fptr, *parFile;
  char *avgfile, fname[256]="", parameterSpace[256]="";
  chromatin c;
  parameters p;
  gillespie g;
  record r;
  quantification q;
  long i, j, locus;
  double FIRING, P_DEMETHYLATE, P_METHYLATE;
  double ALPHA = 1.0, ALPHA_INITIAL = 1.0, BETA = 1.0, BETA_INITIAL = 1.0;
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
  p.loci = 1;
  p.maxReact = 400000;
  p.samples = 400000; 
  p.sampleFreq = p.maxReact/p.samples;

  /* Set program run parameters */
  p.cellCycles = 16;
  p.initialCellCycles = 5;
  p.cellCycleDuration = 22.0; // (hours)
  p.optimSteps = 1; 

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
  p.countFiringEvents = TRUE;
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
  
  /* After initial cell cycles, alpha and beta are switched to the
     values specified on the command line. Store these command line
     alpha and beta for later dynamic changes, then reset p.alpha and
     p.beta to default values */
  ALPHA = p.alpha;
  BETA = p.beta;
  p.alpha = ALPHA_INITIAL;
  p.beta = BETA_INITIAL;

  /* open results file and write header */
  strcpy(parameterSpace,"ParamOptimRes_\0"); strcat(parameterSpace,avgfile); 
  parFile = fopen(parameterSpace,"w");
  fprintParameterSpaceHeader(parFile);

  if (g.test==TRUE)
    g.test_fptr = fopen("TestGillespieFromMain.txt","w");

  /* allocate memory and initialise gillespie algorithm */
  allocateGillespieMemory(&c,&p,&g,&r);
  initialiseGillespieFunctions(&c,&g);

  /* -------------------------------------------------------------------------------- */
  /* Start loop over parameters */
  /* -------------------------------------------------------------------------------- */
  for (p1=0;p1<1;p1++) { // 7
    for (p2=0;p2<p.optimSteps;p2++) { // 20
      for (p3=0;p3<p.optimSteps;p3++) { // 11
	  
        // !!! Set seed for debugging - remove for simulations
        //setseed(&p,p.seed);
        
        // FIRING = 0.000277778*pow(2,p1);
        // P_DEMETHYLATE = pow(10,-0.15*(p2+4));
        // P_METHYLATE = pow(10,-0.12*(p3+26));

        // ALPHA = 0.2*(p2+1);
        // BETA = 1.0 + 0.05*(p3-5); 
        
        FIRING = 0.0001*40;
        P_DEMETHYLATE = 0.004;
        P_METHYLATE = 0.000008;

        // Transcription
        // -------------
        /* Leave the repressed firing rate fixed at ~ every 60 min. */
        p.firingRateMin = 0.0001; 
        p.firingRateMax = FIRING; // Optimise

        /* Cap firing rate at ~ every minute. */
        p.firingCap = 0.0166667;
        p.transcription_demethylate = P_DEMETHYLATE; 
        // p.transcription_turnover = 0.0; (defined in parse.c)

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
        p.noisy_demethylate = P_DEMETHYLATE*p.firingRateMin;
                
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
          
          /* reset counters */
          p.reactCount = 0;
          p.sampleCount = 0;
          p.cellCycleCount = 0;
          resetFiringRecord(&r);
          
          // reset alpha/beta
          p.alpha = ALPHA_INITIAL;
          p.beta = BETA_INITIAL;

          // Schedule first instance of the fixed time reactions
          g.t_nextRep = p.cellCycleDuration*3600;
          p.firingFactor = 1.0;
          
          /* Reaction loop */
          for (i=0;i<p.maxReact && p.cellCycleCount <= p.cellCycles;i++) {
            if (p.cellCycleCount > p.initialCellCycles) {
              p.alpha = ALPHA;
              p.beta = BETA;
            }
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

          accumulateQuantification(&c,&p,&r,&q);
        } /* end loop over loci */

        averageQuantification(&c,&p,&r,&q);
        fprintParameterSpaceResults(parFile,&p,&c,&q);
      }
    }
  }
  
  /* end loop over parameters */
  fclose(parFile);

  /* print final results */
  if (p.resultsFinalLocus == TRUE)
    fprintResultsFinalLocus(avgfile,&r);

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
