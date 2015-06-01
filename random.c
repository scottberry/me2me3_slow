/* 
   Functions to interface the Gnu Scientific Library pseudo-random
   number generator.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
 */

/* Choose a random number in the range [0,1) */
double runif(gsl_rng *r) {
  return((double)gsl_rng_uniform(r));
}

/* Seed the default GSL random number generator */
void rseed(parameters *p) {
  gsl_rng_env_setup();
  if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
  p->gsl_T = gsl_rng_default;
  p->gsl_r = gsl_rng_alloc(p->gsl_T);
  return;
}

/* Seed the default GSL random number generator reproducibly */
void setseed(parameters *p) {
  gsl_rng_env_setup();
  if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = 0;
  p->gsl_T = gsl_rng_default;
  p->gsl_r = gsl_rng_alloc(p->gsl_T);
  return;
}

/* Free the GSL random number generator */
void rfree(parameters *p) {
  gsl_rng_free(p->gsl_r);
  return;
}

