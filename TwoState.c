#include "definitions.h"
#include "model.c"

int main(int argc, char *argv[]) {
  FILE *fptr, *parFile;
  char avgfile[128]="", fname[128]="", tmp[128]="";
  char parameterSpace[128]="", ptmp[128]="";
  chromatin c;
  parameters p;
  gillespie g;
  record r;
  double new, old, gap, Mavg;
  long i, j, locus;
  double R_OFF, FIRING, P_OFF, P_DEMETHYLATE, ENZYMATIC;
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

  /* -------------------------------------------------------------------------------- */
  /* Simulation setup  */
  /* -------------------------------------------------------------------------------- */
  c.sites = 60;

  p.loci = 100;
  p.maxReact = 10000;
  p.samples = 2000;
  p.sampleFreq = p.maxReact/p.samples;

  p.results = TRUE;
  p.optimSteps = 2;

  P_OFF = 0.0;	  

  /* Seed RNG */
  rseed(&p);

  /* Handle filename */
  sprintf(tmp,"s%ld",c.sites); strcat(avgfile,tmp); 
  sprintf(tmp,"r%ldtr0_",p.maxReact); strcat(avgfile,tmp);
  sprintf(ptmp,"%5.3f",P_OFF); 
  sprintf(tmp,"%*s",3,ptmp+2); strcat(avgfile,tmp);
  sprintf(tmp,"st%ld",p.optimSteps); strcat(avgfile,tmp); 
  strcat(avgfile,".txt\0");

  strcpy(parameterSpace,"ParamOptimRes_\0"); strcat(parameterSpace,avgfile); 

  parFile = fopen(parameterSpace,"w");
  fprintf(parFile,"R_OFF\tENZYMATIC\tFIRING\tP_OFF\tP_DEMETHYLATE\tgap\tMavg\n");
  fprintf(stderr,"R_OFF\tENZYMATIC\tFIRING\tP_OFF\tP_DEMETHYLATE\tgap\tMavg\n");

  /* Memory allocation */
  c.state = i_vec_get( c.sites );
  g.bindRep_index = i_vec_get( c.sites );
  g.unbindRep_index = i_vec_get( c.sites );
  g.methylate_index = i_vec_get( c.sites );
  g.demethylate_index = i_vec_get( c.sites );
  g.transcribeDNA_index = i_vec_get( 1 );
  g.propensity = d_vec_get( 4*c.sites + 1 );
  g.doReaction = malloc(g.propensity->len*sizeof( func_ptr_t ) );
  g.doReactionParam = i_vec_get( g.propensity->len );
  g.update = malloc(sizeof( flags ) );

  if (p.results == TRUE) {
    r.t = d_vec_get(p.maxReact);
    r.firing = i_vec_get(p.maxReact);
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
	for (p4=0;p4<p.optimSteps;p4++) {
	  
	  R_OFF = pow(10,-0.25*p1); // log scaling
	  FIRING = pow(10,-0.25*p3); // log scaling
	  P_DEMETHYLATE = fabs(1.0-(double)(p4+1)/(p.optimSteps)); // between 0 and 1
	  ENZYMATIC = pow(10,-0.25*p2); // log scaling

	  // test parameters
	  /*
	  R_OFF = 1.0;
	  FIRING = 1.0;
	  P_DEMETHYLATE = 1.0;
	  P_OFF = 0.0;
	  ENZYMATIC = 0.316;
	  */
	  
	  // Protein binding 
	  // ------------------------------------------------------------
	  p.noisy_Rep_ON = 0.0333; // Leave this fixed at ~ PRC2 tries to bind each site every 30 seconds
  
	  // Stabilisation of Repressors by H3K27me3
	  // ------------------------------------------------------------
	  p.noisy_UR_Rep_OFF = R_OFF; // Optimise
	  p.noisy_MR_Rep_OFF = R_OFF/30; // Leave this fixed at a 30-fold stabilisation of proteins by M marks 
  
	  // Transcription
	  // ------------------------------------------------------------
	  p.firingRateMin = 0.000833; // Leave the repressed firing rate fixed at ~ every 20 min.
	  p.firingRateMax = FIRING; // Optimise
	  p.transcription_RepOFF = P_OFF; // (rate per site per transcription event)
	  p.transcription_demethylate = P_DEMETHYLATE; // (rate per site per transcription event)
  
	  // Methylation/demethylation
	  // ------------------------------------------------------------
	  p.noisy_demethylate = 0.0;
	  p.UR_methylate = ENZYMATIC; // Optimise
	  p.MR_methylate = 10*ENZYMATIC; // Leave this fixed at a 10-fold activation of proteins by M marks 
  
	  gap = 0;
	  Mavg = 0;

	  /* -------------------------------------------------------------------------------- */
	  /* loop over loci */
	  /* -------------------------------------------------------------------------------- */

	  for (locus=0;locus<p.loci;locus++) {
	    // fprintf(stderr,"locus %ld\n",locus);
	    if (argc > 1) {
	      if (strcmp(argv[1],"M")==0)
		initialiseRepressed(&c);
	      else if (strcmp(argv[1],"U")==0)
		initialiseActive(&c);
	    } else {
	      initialiseRandom(&c,&p);
	    }
	    p.reactCount = 0;
	    p.sampleCount = 0;
	    old = 0;
    
	    for (i=0;i<p.maxReact;i++) {
	      if (p.results == TRUE) {
		if (p.reactCount % p.sampleFreq == 0) {
	  
		  for (j=0;j<(c.sites);j++) {
		    r.t_out->el[p.sampleCount] = r.t->el[p.reactCount];
		    r.state->el[j][p.sampleCount] = c.state->el[j];
		  }
		  p.sampleCount++;
		}
	      }
	      p.reactCount++;
	      gillespieStep(&c,&p,&g,&r);
      
	      /* handle DNA replication deterministically, once per day */
	      new = fmod(r.t->el[p.reactCount],86400);
      
	      if (new < old) {
		replicateDNA(&c,&p,g.update);
		//fprintf(stderr,"old = %0.4f, new = %0.4f\n",old,new);
	      }
	      old = new;
	    }
    
	    /* calculate the "gap" for this locus */
	    gap += tAverageGap(&c,&p,&r);
	    Mavg += tAverageM(&c,&p,&r);

	  } /* end loop over loci */

	  fprintf(parFile,"%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n",R_OFF,
		  ENZYMATIC,FIRING,P_OFF,P_DEMETHYLATE,gap/p.loci,Mavg/p.loci);
	  fprintf(stderr,"%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n",R_OFF,
		  ENZYMATIC,FIRING,P_OFF,P_DEMETHYLATE,gap/p.loci,Mavg/p.loci);
	}
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
  free(g.bindRep_index);
  free(g.unbindRep_index);
  free(g.methylate_index);
  free(g.demethylate_index);
  free(g.transcribeDNA_index);
  free(g.propensity);
  free(g.doReaction);
  free(g.update);

  /* print results for final locus */
  strcpy(fname,"t_\0"); strcat(fname,avgfile);
  fptr = fopen(fname,"w");
  d_vec_print(fptr,r.t_out);
  fclose(fptr);

  strcpy(fname,"Meth_t_\0"); strcat(fname,avgfile);
  fprint_nMethylated_t(fname,r.state);

  strcpy(fname,"RepBound_t_\0"); strcat(fname,avgfile);
  fprint_nRepBound_t(fname,r.state);

  strcpy(fname,"Firing_t_\0"); strcat(fname,avgfile);
  fprint_firing_t(fname,&r);

  if (p.results==TRUE) {
    i_vec_free(r.firing);
    i_mat_free(r.state);
    d_vec_free(r.t);
    d_vec_free(r.t_out);
  }

  if (p.testProb == TRUE) {
    free(r.events);
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
  
  return(1);
}
