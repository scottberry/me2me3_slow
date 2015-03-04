#include "definitions.h"
#include "model.c"

int main(int argc, char *argv[]) {
  FILE *fptr;
  char avgfile[128]="", fname[128]="", tmp[128]="";
  chromatin c;
  parameters p;
  record r;
  double new, old;
  int i, j;

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

  /* parameters */

  c.sites = 60;

  // Fast dynamics for protein binding/unbinding 
  // ------------------------------------------------------------
  p.noisy_RepON = 0.015; // 10-fold slower binding than unbinding
  p.noisy_RepOFF = 0.15;

  // Fast dynamics for LHP1
  // ------------------------------------------------------------
  p.M_LHP1_ON = 0.4; 
  p.noisy_LHP1_OFF = 0.4;
  p.stabilised_LHP1_OFF = 0.04; // 10-fold (nearest-neighbour) stabilisation

  // Recruitment of Repressors by H3K27me3 and LHP1 (same value as noisy_OFF)
  // ------------------------------------------------------------
  p.M_bindRep = 0.4;
  p.LHP1_bindRep = 0.4;

  // Removal of repressors by transcription (U)
  // ------------------------------------------------------------
  p.U_unbindRep = 0.05;

  // Relatively slow dynamics for methylation/demethylation
  // ------------------------------------------------------------
  p.noisy_demethylate = 0.0005;
  p.U_demethylate = 0.005; // max. 10-fold change in nucleosome turnover between U and M states (global)

  p.UR_methylate = 0.000125;
  p.MR_methylate = 0.0025; // 20-fold allosteric activation of PHD-PRC2 (local)

  p.maxReact = 10000000;
  p.samples = 1000;
  p.sampleFreq = p.maxReact/p.samples;

  p.results = TRUE;

  /* Seed RNG */
  rseed(&p);

  /* Handle filename */
  sprintf(tmp,"s%ld",c.sites); strcat(avgfile,tmp); 
  sprintf(tmp,"r%ld",p.maxReact); strcat(avgfile,tmp);
  strcat(avgfile,".txt\0");

  /* Memory allocation */
  c.state = i_vec_get(c.sites);
  p.bindRep_index = i_vec_get(c.sites);
  p.unbindRep_index = i_vec_get(c.sites);
  p.bindLHP1_index = i_vec_get(c.sites);
  p.unbindLHP1_index = i_vec_get(c.sites);
  p.methylate_index = i_vec_get(c.sites);
  p.demethylate_index = i_vec_get(c.sites);
  p.propensity = d_vec_get(6*c.sites);
  p.doReaction = malloc(6*c.sites*sizeof( func_ptr_t ) );
  p.doReactionParam = i_vec_get(6*c.sites);
  p.update = malloc(sizeof( flags ) );

  if (p.results == TRUE) {
    r.t = d_vec_get(p.maxReact);
    r.t_out = d_vec_get(p.samples);
    r.state = i_mat_get(c.sites,p.samples);
  }

  /* Initialisation */
  initialiseGillespieFunctions(&c,&p);
  if (strcmp(argv[1],"M")==0)
    initialiseRepressed(&c);
  else if (strcmp(argv[1],"A")==0)
    initialiseActive(&c);
  else
    initialiseRandom(&c,&p);
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
    gillespieStep(&c,&p,&r);

    /* handle DNA replication deterministically, once per day */
    new = fmod(r.t->el[p.reactCount],86400);

    if (new < old) {
      replicateDNA(&c,p.update,&p);
      // fprintf(stderr,"old = %0.4f, new = %0.4f\n",old,new);
    }
    old = new;
  }

  /* free all arrays */
  i_vec_free(c.state);
  free(p.bindRep_index);
  free(p.unbindRep_index);
  free(p.bindLHP1_index);
  free(p.unbindLHP1_index);
  free(p.methylate_index);
  free(p.demethylate_index);
  free(p.propensity);
  free(p.doReaction);
  free(p.update);

  strcpy(fname,"t_\0"); strcat(fname,avgfile);
  fptr = fopen(fname,"w");
  d_vec_print(fptr,r.t_out);
  fclose(fptr);

  strcpy(fname,"Meth_t_\0"); strcat(fname,avgfile);
  fprint_nMethylated_t(fname,r.state);

  strcpy(fname,"RepBound_t_\0"); strcat(fname,avgfile);
  fprint_nRepBound_t(fname,r.state);

  strcpy(fname,"LHP1Bound_t_\0"); strcat(fname,avgfile);
  fprint_nLHP1Bound_t(fname,r.state);

  if (p.results==TRUE) {
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

/*
  c.sites = 60;
  p.noisy_RepON = 0.02;
  p.noisy_RepOFF = 0.12;
  p.noisy_demethylate = 0.0001;
  p.M_LHP1_ON = 0.1;
  p.noisy_LHP1_OFF = 0.2;
  p.stabilised_LHP1_OFF = 0.025;
  p.U_demethylate = 0.003;
  p.U_LHP1_OFF = 0.015;

  p.UR_methylate = 0.0005;
  p.MR_methylate = 0.002;
  p.M_bindRep = 0.2;
  p.LHP1_bindRep = 0.1;
*/

/*
  // Before adding U-mediated protein unbinding

  // Fast dynamics for protein binding/unbinding 
  // ------------------------------------------------------------
  p.noisy_RepON = 0.015; // 10-fold slower binding than unbinding
  p.noisy_RepOFF = 0.15;

  // Fast dynamics for LHP1
  // ------------------------------------------------------------
  p.M_LHP1_ON = 0.4; 
  p.noisy_LHP1_OFF = 0.4;
  p.stabilised_LHP1_OFF = 0.04; // 10-fold (nearest-neighbour) stabilisation

  // Recruitment of Repressors by H3K27me3 and LHP1 (same value as noisy_OFF)
  // ------------------------------------------------------------
  p.M_bindRep = 0.4;
  p.LHP1_bindRep = 0.4;

  // Relatively slow dynamics for methylation/demethylation
  // ------------------------------------------------------------
  p.noisy_demethylate = 0.0005;
  p.U_demethylate = 0.005; // max. 10-fold change in nucleosome turnover between U and M states (global)

  p.UR_methylate = 0.000125;
  p.MR_methylate = 0.0025; // 20-fold allosteric activation of PHD-PRC2 (local)
*/
