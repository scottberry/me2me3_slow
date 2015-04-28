#include "definitions.h"
#include "model.c"

int main(int argc, char *argv[]) {
  FILE *fptr, *parFile;
  char avgfile[128]="", fname[128]="", tmp[128]="", buffer[128]="";
  char parameterSpace[128]="", ptmp[128]="";
  chromatin c;
  parameters p;
  gillespie g;
  record r;
  signed char initial;
  double new, old, gap, Mavg, tTot, tTotM, tTotU, tM, tU, lifetime;
  double firstPassage, firstPassageM, firstPassageU, fpU, fpM;
  long i, j, locus, fh, initM, initU;
  double probM, probU, bistability;
  double R_OFF, FIRING, P_OFF, P_DEMETHYLATE, P_METHYLATE, ENZYMATIC;
  int p1, p2, p3, p4;

  /* Code timing */
#ifdef __APPLE__
  mach_timebase_info_data_t info;
  uint64_t start, end, timeElapsed;
  mach_timebase_info(&info);
  start= mach_absolute_time();
#else
  clock_t start, end, elapsed;
  float timeElapsed;
  start = clock();
#endif

  /* NOTE: DNA REPLICATION OFF */

  /* -------------------------------------------------------------------------------- */
  /* Simulation setup  */
  /* -------------------------------------------------------------------------------- */
  c.sites = 60;

  p.loci = 50;
  p.maxReact = 500000;
  p.samples = 2000;
  p.sampleFreq = p.maxReact/p.samples;

  p.results = TRUE;

  /* ensure that firing_max does not fall below firing_min */
  p.optimSteps = 2; 

  if (argc > 1 && strcmp(argv[1],"P_OFF")==0)
    P_OFF = atof(argv[2]);
  else
    P_OFF = 0.0;	  

  /* Seed RNG */
  rseed(&p);
  // setseed(&p);

  /* Handle filename */
  sprintf(tmp,"s%ld",c.sites); strcat(avgfile,tmp); 
  sprintf(tmp,"r%ldtr0_",p.maxReact); strcat(avgfile,tmp);
  sprintf(ptmp,"%5.3f",P_OFF); 
  sprintf(tmp,"%*s",3,ptmp+2); strcat(avgfile,tmp);
  sprintf(tmp,"st%ld",p.optimSteps); strcat(avgfile,tmp); 
  strcat(avgfile,".txt\0");

  strcpy(parameterSpace,"ParamOptimRes_\0"); strcat(parameterSpace,avgfile); 

  parFile = fopen(parameterSpace,"w");
  fprintf(parFile,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\tP_DEMETHYLATE\tP_METHYLATE\tgap\tMavg\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU\tavgInitU\ttTot\tprobM\tprobU\tbistability\n");
  fprintf(stderr,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\tP_DEMETHYLATE\tP_METHYLATE\tgap\tMavg\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU\tavgInitU\ttTot\tprobM\tprobU\tbistability\n");

  /* Memory allocation */
  c.state = i_vec_get( c.sites );
  g.methylate_index = i_vec_get( c.sites );
  g.transcribeDNA_index = i_vec_get( 1 );
  g.propensity = d_vec_get( c.sites + 1 );
  g.doReaction = malloc(g.propensity->len*sizeof( func_ptr_t ) );
  g.doReactionParam = i_vec_get( g.propensity->len );
  g.update = malloc(sizeof( flags ) );

  if (p.results == TRUE) {
    r.t = d_vec_get(p.maxReact + 1);
    r.firing = i_vec_get(p.maxReact + 1);
    r.t_out = d_vec_get(p.samples);
    r.state = i_mat_get(c.sites,p.samples);
  }

  /* Initialisation */
  initialiseGillespieFunctions(&c,&g);

  /* -------------------------------------------------------------------------------- */
  /* Start loop over parameters */
  /* -------------------------------------------------------------------------------- */
  for (p1=0;p1<p.optimSteps;p1++) {
    for (p2=0;p2<p.optimSteps;p2++) {
      for (p3=0;p3<p.optimSteps;p3++) {
	  
        // !!! Set seed for debugging - remove for simulations
        // setseed(&p);
        
        FIRING = pow(10,-0.3*(p1+4));
        P_DEMETHYLATE = pow(10,-0.3*(p2+2));
        P_METHYLATE = pow(10,-0.3*(p3+5));
        
        /*
        FIRING = 0.63;
        P_DEMETHYLATE = 0.25;
        P_METHYLATE = 0.15;
        */
        
        // Transcription
        // ------------------------------------------------------------
        p.firingRateMin = 0.0005; // Leave the repressed firing rate fixed at ~ every 20 min.
        p.firingRateMax = FIRING; // Optimise
        p.transcription_demethylate = P_DEMETHYLATE; // (rate per site per transcription event)
        if (p.firingRateMax < p.firingRateMin) {
          fprintf(stderr,"Error: Max firing rate less than min firing rate. Setting k_min = k_max\n");
          p.firingRateMin = p.firingRateMax;
        }

        // Methylation/demethylation
        // ------------------------------------------------------------
        p.noisy_methylate = P_METHYLATE/20.0; // 5% noise
        p.me0_me1 = 10*P_METHYLATE;
        p.me1_me2 = P_METHYLATE;
        p.me2_me3 = P_METHYLATE;
        p.me2factor = 0.1;
        p.me3factor = 1;
  
        gap = 0.0;
        Mavg = 0.0;
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
            if (strcmp(argv[1],"M")==0) {
              initialiseRepressed(&c);
	      } else if (strcmp(argv[1],"U")==0) {
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
          p.reactCount = 0;
          p.sampleCount = 0;
          old = 0;
    
          // fprintf(stderr,"Starting reaction loop\n");
          for (i=0;i<p.maxReact;i++) {
            if (p.results == TRUE) {
              if (p.reactCount % p.sampleFreq == 0) {
                // fprintf(stderr,"Sample %ld\n",p.sampleCount);
                for (j=0;j<(c.sites);j++) {
                  r.t_out->el[p.sampleCount] = r.t->el[p.reactCount];
                  r.state->el[j][p.sampleCount] = c.state->el[j];
                }
                p.sampleCount++;
              }
            }
            // fprintf(stderr,"Reaction %ld\n",p.reactCount);
            p.reactCount++;
            gillespieStep(&c,&p,&g,&r);
      
            /* handle DNA replication deterministically, once per day */
            new = fmod(r.t->el[p.reactCount],86400);
      
            if (new < old) {
              // replicateDNA(&c,&p,g.update);
              //fprintf(stderr,"old = %0.4f, new = %0.4f\n",old,new);
            }
            old = new;
          }
          // fprintf(stderr,"Exiting reaction loop\n");
    
          /* calculate and accumulate results for this locus */
          gap += tAverageGap(&c,&p,&r);
          Mavg += tAverage_me2_me3(&c,&p,&r);
          probM += prob_me2_me3(&c,&p,&r);
          probU += prob_me0_me1(&c,&p,&r);
          fh += numberHistoneStateFlips(&r);
          tTot += r.t->el[p.reactCount];

          if (isnan(gap)) {
            fprintf(stderr,"Error: gap is nan. Locus %ld\n",locus);
            fprintf(stderr,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\tP_DEMETHYLATE\tP_METHYLATE\n");
            fprintf(stderr,"%0.6f  %0.6f  %0.6f  %0.6f  %0.6f  %0.6f  %0.6f  %0.6f\n",
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

        fprintf(parFile,"%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.4f\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n",
                p.me0_me1,p.me1_me2,p.me2_me3,p.me2factor,p.me3factor,FIRING,P_DEMETHYLATE,P_METHYLATE,
                gap/p.loci,Mavg/p.loci,lifetime,initM,fpM,tM,initU,fpU,tU,tTot/p.loci,probM/p.loci,probU/p.loci,bistability);
        fprintf(stderr,"%0.6f  %0.6f  %0.6f  %0.6f  %0.6f  %0.6f  %0.6f  %0.6f  %0.4f  %0.4f  %0.4f  %ld  %0.4f  %0.4f  %ld  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f\n",
                p.me0_me1,p.me1_me2,p.me2_me3,p.me2factor,p.me3factor,FIRING,P_DEMETHYLATE,P_METHYLATE,
                gap/p.loci,Mavg/p.loci,lifetime,initM,fpM,tM,initU,fpU,tU,tTot/p.loci,probM/p.loci,probU/p.loci,bistability);
      }
    }
  }
 
  /* end loop over parameters */
  fclose(parFile);

  /* -------------------------------------------------------------------------------- */
  /* Tidy up and write results files */
  /* -------------------------------------------------------------------------------- */

  /* free all arrays */
  i_vec_free(c.state);
  i_vec_free(g.methylate_index);
  i_vec_free(g.transcribeDNA_index);
  d_vec_free(g.propensity);
  free(g.doReaction);
  i_vec_free(g.doReactionParam);
  free(g.update);
  rfree(&p);

  /* print results for final locus */
  strcpy(fname,"t_\0"); strcat(fname,avgfile);
  fptr = fopen(fname,"w");
  d_vec_print(fptr,r.t_out);
  fclose(fptr);

  strcpy(fname,"me0_t_\0"); strcat(fname,avgfile);
  fprint_t(fname,r.state,me0);
  strcpy(fname,"me1_t_\0"); strcat(fname,avgfile);
  fprint_t(fname,r.state,me1);
  strcpy(fname,"me2_t_\0"); strcat(fname,avgfile);
  fprint_t(fname,r.state,me2);
  strcpy(fname,"me3_t_\0"); strcat(fname,avgfile);
  fprint_t(fname,r.state,me3);
  strcpy(fname,"Firing_t_\0"); strcat(fname,avgfile);
  fprint_firing_t(fname,&r);

  strcpy(fname,"Log_\0"); strcat(fname,avgfile);
  fptr = fopen(fname,"w");
  writelog(fptr,&c,&p,&r);

  if (p.results==TRUE) {
    i_vec_free(r.firing);
    i_mat_free(r.state);
    d_vec_free(r.t);
    d_vec_free(r.t_out);
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
  return(1);
}
