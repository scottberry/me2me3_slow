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

/* initialise the protein binding randomly */
void initialiseRepressed(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->state->el[i] = me3;
  }
  return;
}

/* initialise the protein binding randomly */
void initialiseActive(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->state->el[i] = me0;
  }
  return;
}

/* initialise the protein binding randomly */
void initialiseRandom(chromatin *c, parameters *p) {
  int i;
  double rand;
  rand = runif(p->gsl_r);
  //fprintf(stderr,"rand = %0.4f\n",rand);

  if (rand <= 0.5)  {
    for (i=0;i<c->sites;i++) {
      c->state->el[i] = me0;
    }
  } else {
    for (i=0;i<c->sites;i++) {
      c->state->el[i] = me3;
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
void methylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == me0) {
    c->state->el[pos] = me1;
  } else if (c->state->el[pos] == me1) {
    c->state->el[pos] = me2;
  } else if (c->state->el[pos] == me2) {
    c->state->el[pos] = me3;
  }
  update->histone = TRUE;
  return;
}

void demethylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == me3) {
    c->state->el[pos] = me2;
  } else if (c->state->el[pos] == me2) {
    c->state->el[pos] = me1;
  } else if (c->state->el[pos] == me1) {
    c->state->el[pos] = me0;
  }
  update->histone = TRUE;
  return;
}

void transcribeDNA(chromatin *c, parameters *p, flags *update, int pos) {
  unsigned long i;
  for (i=0;i<c->sites;i++) {
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
      c->state->el[pos] = me0;
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

  for (i=0;i<c->sites;i++) { // methylate
    g->doReaction[i] = methylate;
    g->doReactionParam->el[i] = i;
    g->methylate_index->el[i] = i;
  }
  g->doReaction[c->sites] = transcribeDNA; // transcribeDNA
  g->doReactionParam->el[c->sites] = 0;
  g->transcribeDNA_index->el[0] = c->sites;

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

/* Calculate nearest neighbours */
double left(chromatin *c, parameters *p, int pos) {
  double s;
  if (pos>1) {
    if (c->state->el[pos-1] == me2)
      s = p->me2factor;
    else if (c->state->el[pos-1] == me3)
      s = p->me3factor;
  } else {
    s = 0;
  }
  return(s);
}

/* Calculate nearest neighbours */
double right(chromatin *c, parameters *p, int pos) {
  double s;
  if (pos<c->sites-1) {
    if (c->state->el[pos+1] == me2)
      s = p->me2factor;
    else if (c->state->el[pos+1] == me3)
      s = p->me3factor;
  } else {
    s = 0;
  }
  return(s);
}
  
/* Called after each reaction to update the propensities based on the state */
void updatePropensities(chromatin *c, parameters *p, gillespie *g) {
   int i;
   double f_me2_me3;
   
   if (g->update->histone==TRUE) {
     f_me2_me3 = frac(c->state,me2) + frac(c->state,me3);
     //fprintf(stderr,"updatePropensities: update->histone = TRUE, f_me2_me3 = %0.4f\n",f_me2_me3);
     //fprintf(stderr,"left = %0.4f, right = %0.4f\n",left(c,p,3),right(c,p,3));
     
     for (i=0;i<c->sites;i++) {
       if (c->state->el[i] == me0) { // methylate
	 g->propensity->el[g->methylate_index->el[i]] = p->noisy_methylate + p->me0_me1*(left(c,p,i)+right(c,p,i));
       } else if (c->state->el[i] == me1) {
	 g->propensity->el[g->methylate_index->el[i]] = p->noisy_methylate + p->me1_me2*(left(c,p,i)+right(c,p,i));
       } else if (c->state->el[i] == me2) {
	 g->propensity->el[g->methylate_index->el[i]] = p->noisy_methylate + p->me2_me3*(left(c,p,i)+right(c,p,i));
       } else {
         g->propensity->el[g->methylate_index->el[i]] = 0.0;
       }
     }
   
     // transcribeDNA
     g->propensity->el[g->transcribeDNA_index->el[0]] = p->firingFactor*(p->firingRateMax + f_me2_me3*(p->firingRateMin - p->firingRateMax));
     
     g->update->protein = FALSE; // reset the flag
     g->update->histone = FALSE; // reset the flag
     
   }
   /*
     for (i=0;i<g->propensity->len;i++) {
     fprintf(stderr,"i = %d, propensity %0.4f\n",i,g->propensity->el[i]);
     }
   */
   return;
}

/* Single reaction for the Gillespie algorithm */
void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r1=0.0, r2=0.0, scaled_r2=0.0;
  long m, step, i;

  // update and sum propensities, call random numbers
  updatePropensities(c,p,g);

  p_s = d_vec_sum(g->propensity);
  while (r1 == 0.0) { // ensure r1 > 0.0
    r1 = runif(p->gsl_r);
  }
  r2 = runif(p->gsl_r);
  scaled_r2 = p_s*r2;

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

  if (scaled_r2 > g->propensity->el[0]) { // need this "if" to protect against r2=0
    while (scaled_r2 > sum) {
      // fprintf(stderr,"m = %ld, scaled_r2 = %0.4f, sum = %0.4f, propensity = %0.4f\n", m, scaled_r2, sum, g->propensity->el[m]);  
      sum += g->propensity->el[m];
      m++;
    }
    m--;
  } else {
    m=0;
  }
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

void fprint_t(char *fname, I_MAT *mat, int target) {
  FILE *fptr;
  long unsigned count, i, j;
  
  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols;i++) {
    count = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == target) {
	count++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)count/mat->rows);
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
  long sumM, t, pos;
  double gapSum = 0;

  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
	sumM++;
    }
    gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(gapSum/r->t_out->el[p->samples-1]);
}

double prob_me2_me3(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0;
  
  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
	sumM++;
    }
    // fprintf(stderr,"sumM = %ld\t",sumM);
    if (4*sumM > 3*c->sites)
      time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
  }

  // fprintf(stderr,"final time_in_M = %0.4f, total time = %0.4f\n",time_in_M,r->t_out->el[p->samples-1]);
  return(time_in_M/r->t_out->el[p->samples-1]);
}

/* Calculate probability for the last hour of each cell cycle only */
double prob_me2_me3_lastHour(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0;
  double time_total = 0;
  
  for (t=1;t<r->state->cols;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour before replication
      sumM = 0;
      for (pos=0;pos<r->state->rows;pos++) {
        if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
          sumM++;
      }
      // fprintf(stderr,"sumM = %ld\t",sumM);
      if (4*sumM > 3*c->sites) {
        time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
        // fprintf(stderr,"M state: t = %0.2f, sumM = %ld \n",r->t_out->el[t],sumM);
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  // fprintf(stderr,"final time_in_M = %0.4f, total time = %0.4f\n",time_in_M,r->t_out->el[p->samples-1]);
  return(time_in_M/time_total);
}

double prob_me0_me1(chromatin *c, parameters *p, record *r) {
  long sumU, t, pos;
  double time_in_U = 0;
  
  for (t=1;t<r->state->cols;t++) {
    sumU = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me0 || r->state->el[pos][t]==me1)
	sumU++;
    }
    if (4*sumU > 3*c->sites)
      time_in_U += r->t_out->el[t] - r->t_out->el[t-1];
  }

  return(time_in_U/r->t_out->el[p->samples-1]);
}

/* Calculate probability for the last hour of each cell cycle only */
double prob_me0_me1_lastHour(chromatin *c, parameters *p, record *r) {
  long sumU, t, pos;
  double time_in_U = 0;
  double time_total = 0;
  
  for (t=1;t<r->state->cols;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour before replication
      sumU = 0;
      for (pos=0;pos<r->state->rows;pos++) {
        if (r->state->el[pos][t]==me0 || r->state->el[pos][t]==me1)
          sumU++;
      }
      if (4*sumU > 3*c->sites) {
        time_in_U += r->t_out->el[t] - r->t_out->el[t-1];
        // fprintf(stderr,"U state: t = %0.2f, sumU = %ld \n",r->t_out->el[t],sumU);
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  
  return(time_in_U/time_total);
}

double tAverage_me2_me3(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double Mavg = 0;

  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
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
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3) m++;
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
    if (r->state->el[pos][0]==me2 || r->state->el[pos][0]==me3) m++;
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
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3) m++;
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
  fprintf(fptr,"me0_me1: %0.6f\n", p->me0_me1);
  fprintf(fptr,"me1_me2: %0.6f\n", p->me1_me2);
  fprintf(fptr,"me2_me3: %0.6f\n", p->me2_me3);
  fprintf(fptr,"me2factor: %0.6f\n", p->me2factor);
  fprintf(fptr,"me3factor: %0.6f\n", p->me3factor);
  fprintf(fptr,"firingRateMax: %0.6f\n", p->firingRateMax);
  fprintf(fptr,"firingRateMin: %0.6f\n", p->firingRateMin);

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
