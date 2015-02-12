#include "definitions.h"
#include "model.c"

int main(int argc, char *argv[]) {
  FILE *fptr;
  char avgfile[128]="", fname[128]="", tmp[128]="";
  chromatin c;
  parameters p;
  record r;
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
  p.noisy_RepON = 20;
  p.noisy_RepOFF = 120;
  p.noisy_demethylate = 0.1;
  p.M_LHP1_ON = 100;
  p.noisy_LHP1_OFF = 200;
  p.stabilised_LHP1_OFF = 25;
  p.U_demethylate = 3;
  p.U_LHP1_OFF = 15;

  p.UR_methylate = 0.5;
  p.MR_methylate = 2.0;
  p.M_bindRep = 200;
  p.LHP1_bindRep = 100;

  p.maxReact = 1000000;
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
  initialiseRandom(&c,&p);
  p.reactCount = 0;
  p.sampleCount = 0;

  for (i=0;i<p.maxReact;i++) {
    if (p.results == TRUE) {
      if (p.reactCount % p.sampleFreq == 0) {
	// fprintf(stderr,"sampleFreq = %ld, reactCount = %ld, ",p.sampleFreq, p.reactCount);
	for (j=0;j<(c.sites);j++) {
	  r.t_out->el[p.sampleCount] = r.t->el[p.reactCount];
	  r.state->el[j][p.sampleCount] = c.state->el[j];
	}
	// fprintf(stderr,"sampleCount = %ld\n",p.sampleCount);
	p.sampleCount++;
      }
    }
    p.reactCount++;
    gillespieStep(&c,&p,&r);
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
