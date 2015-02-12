/* Choose a random number in the range [0,1) */
double runif(gsl_rng *r) {
  return((double)gsl_rng_uniform(r));
}
/* Choose a random number in the range [0,max) */
double runif_max(gsl_rng *r, double max) {
  return((double)gsl_rng_uniform(r)*max);
}

/* Choose a random number from a Poisson distribution of mean lambda */
long int rpois(gsl_rng *r, double lambda) {
  return((long int)gsl_ran_poisson(r,lambda));
}

/* Choose a random number from an exponential distribution of mean lambda */
long int rexp(gsl_rng *r, double lambda) {
  return((long int)gsl_ran_exponential(r,lambda));
}

/* Choose a random number from an exponential distribution with decay parameter lambda, truncated at max */
long int rexp_trunc(gsl_rng *r, double lambda, long int max) {
  long int l;
  l = (long int)gsl_ran_exponential(r,lambda);
  if (l > max)
    l = max;
  return(l);
}

/* Choose a random number between floor and ceiling inclusive */
long int runif_int(gsl_rng *r, long int floor, long int ceiling) {
  unsigned long int range, v;
  range = ceiling - floor + 1;
  v = gsl_rng_uniform_int (r,range);
  return((long int)(v-floor));
}

/* Seed the default GSL random number generator */
void rseed(parameters *p) {
  gsl_rng_env_setup();
  if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
  p->gsl_T = gsl_rng_default;
  p->gsl_r = gsl_rng_alloc(p->gsl_T);
  return;
}

/* Return a randomly chosen histone on the neighbouring nucleosome */
long int randomNeighbour(parameters *p, unsigned long int num1) {
  long int num2;
  double rnd;
  /* If even can interact with i-2, i-1, i+2, i+3 */
  /* If odd can interact with  i-3, i-2, i+1, i+2 */
  rnd = runif(p->gsl_r);
  if (rnd <= 0.25) {
    if (num1 % 2 == 0)
      num2 = num1 - 2;
    else
      num2 = num1 - 3;
  } else if (rnd <= 0.5) {
    if (num1 % 2 == 0)
      num2 = num1 - 1;
    else
      num2 = num1 - 2;
  } else if (rnd <= 0.75) {
    if (num1 % 2 == 0)
      num2 = num1 + 2;
    else
      num2 = num1 + 1;
  } else {
    if (num1 % 2 == 0)
      num2 = num1 + 3;
    else
      num2 = num1 + 2;
  } 
  return(num2);
}

