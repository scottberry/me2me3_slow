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
  double FIRING, P_DEMETHYLATE, P_METHYLATE, BURST, MEAN;
  int p1, p2, p3;
  D_VEC *meth, *demeth;

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
  p.loci = 1;
  p.maxReact = 1000000;
  p.samples = 1000000; 
  p.sampleFreq = p.maxReact/p.samples;

  /* Set program run parameters */
  p.cellCycles = 20;
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
  p.stochasticAlpha = TRUE;
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
  if (p.stochasticAlpha == TRUE)
    initialiseGillespieFunctionsTransFactor(&c,&g);

  /* -------------------------- */
  /* Start loop over parameters */
  /* -------------------------- */

  meth = d_vec_get(6);
  demeth = d_vec_get(6);

  meth->el[0] = 0.000064;
  meth->el[1] = 0.000032;
  meth->el[2] = 0.000016;
  meth->el[3] = 0.000008;
  meth->el[4] = 0.0000056;
  meth->el[5] = 0.000004;

  demeth->el[0] = 0.1;
  demeth->el[1] = 0.07;
  demeth->el[2] = 0.03;
  demeth->el[3] = 0.004;
  demeth->el[4] = 0.001;
  demeth->el[5] = 0.001;
  
  for (p1=0;p1<1;p1++) { // 7
    for (p2=0;p2<p.optimSteps;p2++) {
      for (p3=0;p3<p.optimSteps;p3++) {
	  
        //setseed(&p,p.seed);
        
        // FIRING = 0.0001*pow(2,p1);
        // P_DEMETHYLATE = pow(10,-0.15*(p2+4));
        // P_METHYLATE = pow(10,-0.15*(p3+19));
             
        FIRING = 0.0001*40.0;
        // P_DEMETHYLATE = demeth->el[p1];
        // P_METHYLATE = meth->el[p1];
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
        // p.transcription_turnover = P_TURNOVER/2.0; (defined in parse.c)

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

        // Stochastic alpha
        // ----------------
        if (p.stochasticAlpha == TRUE) {
          BURST = p.stochasticTranslationEfficiency;
          MEAN = 1000.0;

          /* <p> = k_r * burst / gamma_p,
             where burst = k_p/gamma_r */
          p.gamma_r = 1.0/(2.0*3600.0);
          p.gamma_p = 1.0/(12.0*3600.0);
          p.k_r = p.gamma_p * MEAN / BURST;
          p.k_p = p.gamma_r * BURST;

          p.transFactorRNA = floor(p.k_r/p.gamma_r);
          p.transFactorProtein = floor(MEAN);
        }

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
              if (p.stochasticAlpha == TRUE) {
                r.transFactorProtein->el[p.sampleCount] = p.transFactorProtein;
                r.transFactorRNA->el[p.sampleCount] = p.transFactorRNA;
                r.alpha->el[p.sampleCount] = p.alpha;
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

  d_vec_free(meth);
  d_vec_free(demeth);
  
  /* print final results */
  if (p.resultsFinalLocus == TRUE) {
    fprintResultsFinalLocus(avgfile,&r);

    if (p.stochasticAlpha == TRUE) {
      strcpy(fname,"alpha_\0"); strcat(fname,avgfile);
      fprint_transFactorProtein_nCycles(fname,&r);
      strcpy(fname,"alphaOnly_\0"); strcat(fname,avgfile);
      fprint_alphaOnly_nCycles(fname,&r);
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
