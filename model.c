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

/* initialise the protein binding randomly */
void initialiseRandom(chromatin *c, parameters *p) {
  int i;
  double rand;
  for (i=0;i<c->sites;i++) {
    rand = runif(p->gsl_r);
    if (rand < 1.0/2.0) 
      c->state->el[i] = U;
    else 
      c->state->el[i] = M;
  }

  return;
}

double d_vec_sum(D_VEC *d) {
  int i;
  double sum = 0;
  for (i=0;i<d->len;i++) {
    sum += d->el[i];
  }
  return(sum);
}

/* Individual reactions */
void bindRep(chromatin *c, flags *update, int pos) {
  if (c->state->el[pos] == U) {
    c->state->el[pos] = UR;
  } else if (c->state->el[pos] == M) {
    c->state->el[pos] = MR;
  }
  update->protein = TRUE;
  return;
}

void unbindRep(chromatin *c, flags *update, int pos) {
  if (c->state->el[pos] == UR) {
    c->state->el[pos] = U;
  } else if (c->state->el[pos] == MR) {
    c->state->el[pos] = M;
  }
  update->protein = TRUE;
  return;
}

void methylate(chromatin *c, flags *update, int pos) {
  if (c->state->el[pos] == U) {
    c->state->el[pos] = M;
  } else if (c->state->el[pos] == UR) {
    c->state->el[pos] = MR;
  }
  update->histone = TRUE;
  return;
}

void demethylate(chromatin *c, flags *update, int pos) {
  if (c->state->el[pos] == M) {
    c->state->el[pos] = U;
  } else if (c->state->el[pos] == MR) {
    c->state->el[pos] = UR;
  }
  update->histone = TRUE;
  return;
}

/* Called only once to initialise indices and function pointers */
void initialiseGillespieFunctions(chromatin *c, parameters *p) {
  int i;

  p->update->protein = TRUE;
  p->update->histone = TRUE;

  for (i=0;i<c->sites;i++) { // Rep binding
    p->doReaction[i] = bindRep; // point the function doReaction->el[i]
    p->doReactionParam->el[i] = i; // store the site number in a corresponding vector
    p->bindRep_index->el[i] = i; // index the address in the propensity/reaction array for this reaction
  }
  for (i=c->sites;i<2*c->sites;i++) { // Rep unbinding
    p->doReaction[i] = unbindRep;
    p->doReactionParam->el[i] = i-c->sites;
    p->unbindRep_index->el[i-c->sites] = i;
  }
  for (i=2*c->sites;i<3*c->sites;i++) { // methylate
    p->doReaction[i] = methylate;
    p->doReactionParam->el[i] = i-2*c->sites;
    p->methylate_index->el[i-2*c->sites] = i;
  }
  for (i=3*c->sites;i<4*c->sites;i++) { // demethylate
    p->doReaction[i] = demethylate;
    p->doReactionParam->el[i] = i-3*c->sites;
    p->demethylate_index->el[i-3*c->sites] = i;
  }
  return;
}

/* Calculate the fraction of "target" reactants */
double frac(I_VEC *vec, int target) {
  unsigned long int count = 0;
  double f;
  int i;
  for (i=0;i<vec->len;i++) {
    if (vec->el[i]==target) count++;
  }
  f = (double)count/vec->len;
  return(f);
}

/* Called after each reaction to update the propensities based on the state */
void updatePropensities(chromatin *c, parameters *p) {
  int i;
  double f_UR, f_MR, f_M, f_U;

  if (p->update->protein==TRUE | p->update->histone==TRUE) {

    f_UR = frac(c->state,UR);
    f_MR = frac(c->state,MR);
    f_M = frac(c->state,M);
    f_U = frac(c->state,U);

    for (i=0;i<c->sites;i++) { // bindRep
      if (c->state->el[i]==U || c->state->el[i]==M) {
	p->propensity->el[p->bindRep_index->el[i]] = p->noisy_RepON + (f_M + f_MR)*p->M_bindRep;
      } else {
	p->propensity->el[p->bindRep_index->el[i]] = 0.0;
      }
      if (c->state->el[i]==UR || c->state->el[i]==MR) { // unbindRep
	p->propensity->el[p->unbindRep_index->el[i]] = p->noisy_RepOFF;
      }	else {
	p->propensity->el[p->unbindRep_index->el[i]]= 0.0;
      }
      if (c->state->el[i]==UR || c->state->el[i]==U) { // methylate
	p->propensity->el[p->methylate_index->el[i]] = f_UR*p->UR_methylate + pow(f_MR,2)*p->MR_methylate/(pow(f_MR,2)+0.5);
      }	else {
	p->propensity->el[p->methylate_index->el[i]] = 0.0;
      }
      if (c->state->el[i]==M || c->state->el[i]==MR) { // demethylate
	p->propensity->el[p->demethylate_index->el[i]] = p->noisy_demethylate + (f_U + f_UR)*p->U_demethylate;
      } else {
	p->propensity->el[p->demethylate_index->el[i]] = 0.0;
      }
    }
    p->update->protein = FALSE; // reset the flag
    p->update->histone = FALSE; // reset the flag
  }
  return;
}

/* Single reaction for the Gillespie algorithm */
void gillespieStep(chromatin *c, parameters *p, record *r) {
  double delta_t, sum, p_s, r1, scaled_r2;
  long m, step, i;

  // update and sum propensities, call random numbers
  updatePropensities(c,p);

  p_s = d_vec_sum(p->propensity);
  r1 = runif(p->gsl_r);
  scaled_r2 = p_s*runif(p->gsl_r);

  // calculate time step
  delta_t = log(1.0/r1)/p_s;
  step = p->reactCount;
  r->t->el[step] = delta_t + r->t->el[step-1];

  // choose reaction m from propensities based on scaled_r2
  sum = 0;
  m = 0;
  /* 
  for (i=0;i<p->propensity->len;i++) {
    fprintf(stderr,"propensity = %0.4f\n",p->propensity->el[i]);  
  }
  */

  while (scaled_r2 > sum) {
    // fprintf(stderr,"m = %ld, scaled_r2 = %0.4f, sum = %0.4f, propensity = %0.4f\n", m, scaled_r2, sum, p->propensity->el[m]);  
    sum += p->propensity->el[m];
    m++;
  }
  m--;
  p->doReaction[m](c,p->update,p->doReactionParam->el[m]);
  return;
}

void fprint_nMethylated_t(char *fname, I_MAT *mat) {
  FILE *fptr;
  long unsigned methyl, i, j;
  
  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols;i++) {
    methyl = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == M || mat->el[j][i] == MR ) {
	methyl++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)methyl/mat->rows);
  }
  fclose(fptr);
  return;
}

void fprint_nRepBound_t(char *fname, I_MAT *mat) {
  FILE *fptr;
  long unsigned bound, i, j;

  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols;i++) {
    bound = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == UR || mat->el[j][i] == MR ) {
	bound++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)bound/mat->rows);
  }
  fclose(fptr);
  return;
}

