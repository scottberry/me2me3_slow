/* 
   Implementation of the direct Gillespie algorithm for simulation
   of models of chromatin-based epigenetics.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
 */

void initialiseRepressed(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->K27->el[i] = me3;
  }
  return;
}

void initialiseActive(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->K27->el[i] = me0;
  }
  return;
}

void initialiseRandom(chromatin *c, parameters *p) {
  int i;
  double rand;
  rand = runif(p->gsl_r);

  if (rand <= 0.5)  {
    for (i=0;i<c->sites;i++) {
      c->K27->el[i] = me0;
    }
  } else {
    for (i=0;i<c->sites;i++) {
      c->K27->el[i] = me3;
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

double enzymaticFactor(chromatin *c, int pos) {
  double s;
  
  if (c->K27->el[pos] == me2)
    s = p->me2factor;
  else if (c->K27->el[pos] == me3)
    s = p->me3factor;
  else
    s = 0.0;

  return(s);
}


/* Find histone mods on neighbouring nucleosome (left) */

double neighboursK27factor(chromatin *c, parameters *p, int pos) {
  double n = 0;

  if (pos % 2 == 0) { // (even) top histone tail
    n += enzymaticFactor(c,pos+1); // other tail on same nucleosome
    if (pos != 0) {
      n += enzymaticFactor(c,pos-1); // top tail on left neighbour nucleosome
      n += enzymaticFactor(c,pos-2); // bottom tail on left neighbour nucleosome
    }
    if (pos < c->sites - 3 ) {
      n += enzymaticFactor(c,pos+2); // top tail on right neighbour nucleosome
      n += enzymaticFactor(c,pos+3); // bottom tail on right neighbour nucleosome    
    }
  } else { // (odd) bottom histone tail
    n += enzymaticFactor(c,pos-1); // other tail on same nucleosome
    if (pos != 1) {
      n += enzymaticFactor(c,pos-2); // top tail on left neighbour nucleosome
      n += enzymaticFactor(c,pos-3); // bottom tail on left neighbour nucleosome
    }
    if (pos < c->sites - 2) {
      n += enzymaticFactor(c,pos+1); // top tail on right neighbour nucleosome
      n += enzymaticFactor(c,pos+2); // bottom tail on right neighbour nucleosome    
    }
  }
  return(n);
}

/* Update the propensities based on change in system K27. */

void updatePropensities(chromatin *c, parameters *p, gillespie *g) {
   int i;
   double f_me2_me3;
   
   if (g->update->histone==TRUE) {
     f_me2_me3 = frac(c->K27,me2) + frac(c->K27,me3);
     
     for (i=0;i<c->sites;i++) {
       if (c->K27->el[i] == me0) { // methylate
	 g->propensity->el[g->methylate_index->el[i]] = p->noisy_methylate + p->me0_me1*(neighboursK27factor(c,p,i));
       } else if (c->K27->el[i] == me1) {
	 g->propensity->el[g->methylate_index->el[i]] = p->noisy_methylate + p->me1_me2*(neighboursK27factor(c,p,i));
       } else if (c->K27->el[i] == me2) {
	 g->propensity->el[g->methylate_index->el[i]] = p->noisy_methylate + p->me2_me3*(neighboursK27factor(c,p,i));
       } else {
         g->propensity->el[g->methylate_index->el[i]] = 0.0;
       }
     }
   
     // transcribeDNA
     g->propensity->el[g->transcribeDNA_index->el[0]] =
       p->firingFactor*(p->firingRateMax + f_me2_me3*(p->firingRateMin - p->firingRateMax));
     
     g->update->protein = FALSE; // reset the flag
     g->update->histone = FALSE; // reset the flag
     
   }

   return;
}

/* Single iteration of the "direct" Gillespie algorithm. */

void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r1=0.0, r2=0.0, scaled_r2=0.0;
  long m, step, i;

  // update and sum propensities, call random numbers.
  updatePropensities(c,p,g);

  // Note: runif sometimes calls exactly 0.0, ensure that code is
  // robust to this value.
  
  p_s = d_vec_sum(g->propensity);
  while (r1 == 0.0) { // protect against r1==0.0
    r1 = runif(p->gsl_r);
  }
  r2 = runif(p->gsl_r);
  scaled_r2 = p_s*r2;

  // calculate time step  
  delta_t = log(1.0/r1)/p_s;
  step = p->reactCount;
  r->t->el[step] = delta_t + r->t->el[step-1];

  // choose reaction m from propensities based on scaled_r2
  sum = 0;
  m = 0;

  if (scaled_r2 > g->propensity->el[0]) { // protect against r2==0.0
    while (scaled_r2 > sum) {
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
