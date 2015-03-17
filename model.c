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

/* Free the GSL random number generator */
void rfree(parameters *p) {
  gsl_rng_free(p->gsl_r);
  return;
}

/* initialise the protein binding randomly */
void initialiseRepressed(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->state->el[i] = M;
  }
  return;
}

/* initialise the protein binding randomly */
void initialiseActive(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->state->el[i] = U;
  }
  return;
}

/* initialise the protein binding randomly */
void initialiseRandom(chromatin *c, parameters *p) {
  int i;
  double rand;
  rand = runif(p->gsl_r);
  if (rand <= 0.5)  {
    for (i=0;i<c->sites;i++) {
      c->state->el[i] = U;
    }
  } else {
    for (i=0;i<c->sites;i++) {
      c->state->el[i] = M;
    }
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
void bindRep(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == U) {
    c->state->el[pos] = UR;
  } else if (c->state->el[pos] == M) {
    c->state->el[pos] = MR;
  }
  update->protein = TRUE;
  return;
}

void unbindRep(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == UR) {
    c->state->el[pos] = U;
  } else if (c->state->el[pos] == MR) {
    c->state->el[pos] = M;
  }
  update->protein = TRUE;
  return;
}

void methylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == U) {
    c->state->el[pos] = M;
  } else if (c->state->el[pos] == UR) {
    c->state->el[pos] = MR;
  }
  update->histone = TRUE;
  return;
}

void demethylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == M) {
    c->state->el[pos] = U;
  } else if (c->state->el[pos] == MR) {
    c->state->el[pos] = UR;
  }
  update->histone = TRUE;
  return;
}

void transcribeDNA(chromatin *c, parameters *p, flags *update, int pos) {
  unsigned long i;
  for (i=0;i<c->sites;i++) {
    if(runif(p->gsl_r) <= p->transcription_RepOFF) {
      unbindRep(c,p,update,i);
    }
    if(runif(p->gsl_r) <= p->transcription_demethylate) {
      demethylate(c,p,update,i);
    }
    update->transcribed = TRUE;
  }
  return;
}

void replicateDNA(chromatin *c, parameters *p, flags *update) {
  unsigned long pos;
  for (pos=0;pos<c->sites;pos++) {
    if(runif(p->gsl_r)<=0.5) {
      c->state->el[pos] = U;
    }
  }
  update->protein = TRUE;
  update->histone = TRUE;
  return;
}

/* Called only once to initialise indices and function pointers */
void initialiseGillespieFunctions(chromatin *c, gillespie *g) {
  int i;
  
  g->update->protein = TRUE;
  g->update->histone = TRUE;
  g->update->transcribed = FALSE;

  for (i=0;i<c->sites;i++) { // Rep binding
    g->doReaction[i] = bindRep; // point the function doReaction->el[i]
    g->doReactionParam->el[i] = i; // store the site number in a corresponding vector
    g->bindRep_index->el[i] = i; // index the address in the propensity/reaction array for this reaction
  }
  for (i=c->sites;i<2*c->sites;i++) { // Rep unbinding
    g->doReaction[i] = unbindRep;
    g->doReactionParam->el[i] = i-c->sites;
    g->unbindRep_index->el[i-c->sites] = i;
  }
  for (i=2*c->sites;i<3*c->sites;i++) { // methylate
    g->doReaction[i] = methylate;
    g->doReactionParam->el[i] = i-2*c->sites;
    g->methylate_index->el[i-2*c->sites] = i;
  }
  for (i=3*c->sites;i<4*c->sites;i++) { // demethylate
    g->doReaction[i] = demethylate;
    g->doReactionParam->el[i] = i-3*c->sites;
    g->demethylate_index->el[i-3*c->sites] = i;
  }
  g->doReaction[4*c->sites] = transcribeDNA; // transcribeDNA
  g->doReactionParam->el[4*c->sites] = 0;
  g->transcribeDNA_index->el[0] = 4*c->sites;
  /*
  i_vec_print(stderr,g->bindRep_index);
  i_vec_print(stderr,g->unbindRep_index);
  i_vec_print(stderr,g->methylate_index);
  i_vec_print(stderr,g->demethylate_index); 
  i_vec_print(stderr,g->transcribeDNA_index); 
  */
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
void updatePropensities(chromatin *c, parameters *p, gillespie *g) {
   int i;
   double f_UR, f_MR, f_M, f_U;

  if (g->update->protein==TRUE || g->update->histone==TRUE) {

    f_UR = frac(c->state,UR);
    f_MR = frac(c->state,MR);
    f_M = frac(c->state,M);
    f_U = frac(c->state,U);

    for (i=0;i<c->sites;i++) { // bindRep
      if (c->state->el[i]==U || c->state->el[i]==M) {
	g->propensity->el[g->bindRep_index->el[i]] = p->noisy_Rep_ON;
      } else {
	g->propensity->el[g->bindRep_index->el[i]] = 0.0;
      }
      if (c->state->el[i]==UR) { // unbindRep
	g->propensity->el[g->unbindRep_index->el[i]] = p->noisy_UR_Rep_OFF;
      }	else if (c->state->el[i]==MR) {
	g->propensity->el[g->unbindRep_index->el[i]] = p->noisy_MR_Rep_OFF;
      }	else {
	g->propensity->el[g->unbindRep_index->el[i]]= 0.0;
      }
      if (c->state->el[i]==UR || c->state->el[i]==U) { // methylate
	g->propensity->el[g->methylate_index->el[i]] = f_UR*p->UR_methylate + f_MR*p->MR_methylate;
      }	else {
	g->propensity->el[g->methylate_index->el[i]] = 0.0;
      }
      if (c->state->el[i]==M || c->state->el[i]==MR) { // demethylate
	g->propensity->el[g->demethylate_index->el[i]] = p->noisy_demethylate;
      } else {
	g->propensity->el[g->demethylate_index->el[i]] = 0.0;
      } 
    }
    // transcribeDNA
    g->propensity->el[g->transcribeDNA_index->el[0]] = p->firingRateMax + (f_M + f_MR)*(p->firingRateMin - p->firingRateMax);
    //fprintf(stderr,"f_M + f_MR = %0.4f, firing rate = %0.4f\n",f_M+f_MR,g->propensity->el[g->transcribeDNA_index->el[0]]);

    g->update->protein = FALSE; // reset the flag
    g->update->histone = FALSE; // reset the flag

    /*i_vec_print(stderr,c->state);
    fprintf(stderr,"i = 0, propensity %0.4f\n",g->propensity->el[0]);
    fprintf(stderr,"i = 1, propensity %0.4f\n",g->propensity->el[1]);
    d_vec_print(stderr,g->propensity);*/
  }
  return;
}

/* Single reaction for the Gillespie algorithm */
void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r1, scaled_r2;
  long m, step, i;

  // update and sum propensities, call random numbers
  updatePropensities(c,p,g);

  p_s = d_vec_sum(g->propensity);
  r1 = runif(p->gsl_r);
  scaled_r2 = p_s*runif(p->gsl_r);

  // calculate time step
  delta_t = log(1.0/r1)/p_s;
  step = p->reactCount;
  //fprintf(stderr,"step = %ld\n",step);
  r->t->el[step] = delta_t + r->t->el[step-1];

  // choose reaction m from propensities based on scaled_r2
  sum = 0;
  m = 0;
  /* 
  for (i=0;i<g->propensity->len;i++) {
    fprintf(stderr,"propensity = %0.4f\n",g->propensity->el[i]);  
  }
  */

  while (scaled_r2 > sum) {
    // fprintf(stderr,"m = %ld, scaled_r2 = %0.4f, sum = %0.4f, propensity = %0.4f\n", m, scaled_r2, sum, g->propensity->el[m]);  
    sum += g->propensity->el[m];
    m++;
  }
  m--;
  g->doReaction[m](c,p,g->update,g->doReactionParam->el[m]);

  // record each firing event
  if (g->update->transcribed == TRUE) {
    r->firing->el[p->reactCount] = TRUE;
    g->update->transcribed = FALSE;
  } else {
    r->firing->el[p->reactCount] = FALSE;
  }

  return;
}

void fprint_nMethylated_t(char *fname, I_MAT *mat) {
  FILE *fptr;
  long unsigned methyl, i, j;
  
  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols;i++) {
    methyl = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == M || mat->el[j][i] == MR) {
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
      if (mat->el[j][i] == UR || mat->el[j][i] == MR) {
	bound++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)bound/mat->rows);
  }
  fclose(fptr);
  return;
}

void fprint_firing_t(char *fname, record *r) {
  FILE *fptr;
  long unsigned i;

  fptr = fopen(fname,"w");
  for (i=0;i<r->firing->len;i++) {
    if (r->firing->el[i]==TRUE)
      fprintf(fptr,"%0.4f\n",r->t->el[i]);
  }
  fclose(fptr);
  return;
}

double tAverageGap(chromatin *c, parameters *p, record *r) {
  unsigned long sumM, t, pos;
  double gapSum = 0;

  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==M || r->state->el[pos][t]==MR)
	sumM++;
    }
    gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(gapSum/r->t_out->el[p->samples-1]);
}

double tAverageM(chromatin *c, parameters *p, record *r) {
  unsigned long sumM, t, pos;
  double Mavg = 0;

  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==M || r->state->el[pos][t]==MR)
	sumM++;
    }
    Mavg += (double)sumM*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(Mavg/r->t_out->el[p->samples-1]);
}

/* Calculate the average lifetime of a state (max reported 1 flip per generation) */
unsigned long numberHistoneStateFlips(record *r) {
  signed char newState = 0, oldState = 0;
  unsigned long t, pos, flips=0, m=0, u=0;;
  
  for (t=0;t<r->state->cols;t++) {
    
    oldState = newState;
    
    m = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==M || r->state->el[pos][t]==MR) m++;
    }
    u = r->state->rows - m;

    if (oldState == 1) { // if previously U
      if (m >= 3*u) {
	newState = -1;
	flips++;
      }
    } else if (oldState == -1) { // if previously M
      if (u >= 3*m) {
	newState = 1;
	flips++;
      }
    } else if (oldState == 0) { // first time-point
      if (u >= 3*m)
	newState = 1;
      else
	newState = -1;
    }
  }
 
  return(flips);
}

double firstPassageTime(record *r, signed char *initial) {
  long unsigned m=0, u=0, pos,t=0;
  
  /* find initial state */
  for (pos=0;pos<r->state->rows;pos++) {
    if (r->state->el[pos][0]==M || r->state->el[pos][0]==MR) m++;
  }
  u = r->state->rows - m;
  
  if (m > u)
    *initial = -1;
  else
    *initial = 1;
    
  while ( t < r->state->cols &&
	  ((*initial == -1 && 3*m > u) || (*initial == 1 && 3*u > m))) {
    m = 0;
    u = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==M || r->state->el[pos][t]==MR) m++;
    }
    u = r->state->rows - m;
    //fprintf(stderr,"t = %0.2f, m %ld u %ld \n",r->t_out->el[t],m,u);
    t++;
  }
  return(r->t_out->el[t-1]);
}

/* write a log file */
int writelog(FILE *fptr, chromatin *c, parameters *p, record *r) {
  time_t curtime;
  struct tm *loctime;

  curtime = time (NULL);
  loctime = localtime (&curtime);

  fputs (asctime (loctime), fptr);
#ifdef __APPLE__
  fprintf(fptr,"Operating system: Mac OS\n");
#else
  fprintf(fptr,"Operating system: Unix\n");
#endif
  fprintf(fptr,"sites: %ld\n", c->sites);
  fprintf(fptr,"loci: %ld\n", p->loci);
  fprintf(fptr,"maxReact: %ld\n", p->maxReact);
  fprintf(fptr,"samples: %ld\n", p->samples);
  fprintf(fptr,"max. time (final locus): %0.2f seconds\n\n",r->t->el[r->t->len-1]);

  fprintf(fptr,"noisy_Rep_ON: %0.6f\n", p->noisy_Rep_ON);
  fprintf(fptr,"noisy_UR_Rep_OFF: %0.6f\n", p->noisy_UR_Rep_OFF);
  fprintf(fptr,"noisy_MR_Rep_OFF: %0.6f\n", p->noisy_MR_Rep_OFF);
  fprintf(fptr,"noisy_demethylate: %0.6f\n", p->noisy_demethylate);
  fprintf(fptr,"UR_methylate: %0.6f\n", p->UR_methylate);
  fprintf(fptr,"MR_methylate: %0.6f\n", p->MR_methylate);
  fprintf(fptr,"firingRateMax: %0.6f\n", p->firingRateMax);
  fprintf(fptr,"firingRateMin: %0.6f\n", p->firingRateMin);
  fprintf(fptr,"transcription_RepOFF: %0.6f\n", p->transcription_RepOFF);
  fprintf(fptr,"transcription_demethylate: %0.6f\n", p->transcription_demethylate);

  return(1);
}

/* 
long unsigned left(chromatin *c, long unsigned i) {
  long unsigned l;

  if (i > 0)
    l = i-1;
  else
    l = i;

  return(l);  
}

long unsigned right(chromatin *c, long unsigned i) {
  long unsigned r;

  if (i < c->sites-1)
    r = i+1;
  else
    r = i;

  return(r);  
}

void fprint_nLHP1Bound_t(char *fname, I_MAT *mat) {
  FILE *fptr;
  long unsigned bound, i, j;

  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols;i++) {
    bound = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == M_LHP1 || mat->el[j][i] == MR_LHP1) {
	bound++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)bound/mat->rows);
  }
  fclose(fptr);
  return;
}

*/
