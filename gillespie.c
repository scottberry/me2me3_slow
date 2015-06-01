
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
