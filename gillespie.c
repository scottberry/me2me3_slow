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

double fracControlRegion_me2me3(chromatin *c) {
  unsigned long int count = 0;
  double f;
  int i;
  for (i=0;i<c->controlSites;i++) {
    if (c->K27->el[i]==me2 || c->K27->el[i]==me3) count++;
  }
  f = (double)count/c->controlSites;
  return(f);
}

double enzymaticFactor(chromatin *c, parameters *p, int pos) {
  double s;
  
  if (c->K27->el[pos] == me2)
    s = p->me2factor;
  else if (c->K27->el[pos] == me3)
    s = p->me3factor;
  else
    s = 0.0;

  return(s);
}


/* Find histone mods on neighbouring nucleosome */
double neighboursK27factor(chromatin *c, parameters *p, int pos) {
  double n = 0;

  if (pos % 2 == 0) { // (even) top histone tail
    n += enzymaticFactor(c,p,pos+1); // other tail on same nucleosome
    if (pos != 0) { // check not left-most nucleosome
      n += enzymaticFactor(c,p,pos-1); // top tail on left neighbour nucleosome
      n += enzymaticFactor(c,p,pos-2); // bottom tail on left neighbour nucleosome
    }
    if (pos < c->sites - 3 ) { // check not right-most nucleosome
      n += enzymaticFactor(c,p,pos+2); // top tail on right neighbour nucleosome
      n += enzymaticFactor(c,p,pos+3); // bottom tail on right neighbour nucleosome    
    }
  } else { // (odd) bottom histone tail
    n += enzymaticFactor(c,p,pos-1); // other tail on same nucleosome
    if (pos != 1) { // check not left-most nucleosome
      n += enzymaticFactor(c,p,pos-2); // top tail on left neighbour nucleosome
      n += enzymaticFactor(c,p,pos-3); // bottom tail on left neighbour nucleosome
    }
    if (pos < c->sites - 2) { // check not right-most nucleosome
      n += enzymaticFactor(c,p,pos+1); // top tail on right neighbour nucleosome
      n += enzymaticFactor(c,p,pos+2); // bottom tail on right neighbour nucleosome    
    }
  }
  return(n);
}

/* Update the propensities based on change in system K27. */

void updatePropensities(chromatin *c, parameters *p, gillespie *g) {
  int i;
  double f_me2_me3;
   
  if (g->update->histone==TRUE) {
    f_me2_me3 = fracControlRegion_me2me3(c);
    //f_me2_me3 = frac(c->K27,me2) + frac(c->K27,me3);
     
    for (i=0;i<c->sites;i++) {
      // fprintf(stderr,"i = %d, neighboursK27factor = %0.2f\n",i,neighboursK27factor(c,p,i));
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
      p->activation*p->firingFactor*(p->firingRateMax + f_me2_me3*(p->firingRateMin - p->firingRateMax));
     
    g->update->histone = FALSE; // reset the flag
     
  }

  return;
}

/* calculate the gillespie time step (and propensity sum in argument p_s) */
double gillespieTimeStep(parameters *p, gillespie *g, double *p_s) {
  double delta_t, r1=0.0;

  /* Note: runif sometimes calls exactly 0.0, ensure that code is
     robust to this value. */
  
  (*p_s) = d_vec_sum(g->propensity);
  while (r1 == 0.0) { // protect against r1==0.0
    r1 = runif(p->gsl_r);
  }

  // calculate time step  
  delta_t = log(1.0/r1)/(*p_s);

  /* Need some safeguard to determine if delta_t more than a cell
     cycle, in which case two DNA replications will need to have taken
     place and program should exit. */
  
  if (delta_t > 3600*(p->cellCycleDuration-p->G2duration)) {
    fprintf(stderr,"Error (gillespieStep): delta_t > cell cycle.");
    exit(-1);
  }
  
  return(delta_t); 
}
  
/* Single iteration of the "direct" Gillespie algorithm. Note that
   this varies somewhat from a pure stochastic simulation with
   exponentially distributed waiting times because there are events
   which occur at precise time points (DNA replication and release of
   G2 transcriptional inhibition). Some reactions are discarded as the
   system state has changed between the two simulation steps.
   Consequently the algorithm is non-Markovian. See
   \cite{Barrio:2006kb} for a description of the problems with delays
   in SSA simulations. */

void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r2=0.0, scaled_r2=0.0, next_t;
  long m, step, i;
  logical gillespieReaction = TRUE;

  step = p->reactCount;
  
  // update propensities
  updatePropensities(c,p,g);
  
  // calculate time step (p_s also stores propensity sum).
  delta_t = gillespieTimeStep(p,g,&p_s);
  next_t = delta_t + r->t->el[step-1];

  /* Determine if the new time is after DNA replication. If so,
     interrupt Gillespie algorithm to replicate DNA (at a precise
     time). If not, check if G2 transcription inhibition needs to be
     released and then proceed with a Gillespie algorithm choice for
     next system update. */

  //fprintf(stderr,"next_t = %0.0f ",next_t);
  
  r->new = fmod(next_t,3600*p->cellCycleDuration); 
  //fprintf(stderr,"r->new = %0.0f, r->old %0.0f, cell cycle %d\n",r->new, r->old,p->cellCycleCount);
  if (r->new < r->old) {
    // fprintf(stderr,"replicate\n");
    // DNA replication
    // ------------------------------

    // increment counters and record new (deterministic) time
    p->cellCycleCount++;    
    r->t->el[step] = next_t - r->new;

    // replicate DNA and instigate G2 transcriptional inhibition
    if (p->DNAreplication == TRUE)  {
      replicateDNA(c,p,g->update);
      if (p->G2duration > 0.0) p->firingFactor = 0.5;
    }
    // record time of last replication
    r->t_lastRep = r->t->el[step];

    // ensure that DNA replication does not occur next time
    r->old = 0.0;

    // do not proceed with gillespie reaction this time
    gillespieReaction = FALSE;

  } else if (p->G2duration > 0.0 &&
             (next_t - r->t_lastRep) > (3600*p->G2duration) &&
             p->firingFactor == 0.5) {
    
    // Release G2 firing inhibition
    // ------------------------------
    /* Note that this solution is better but not perfect since the
       next delta_t chosen may actually be less than the current
       delta_t, which will lead effectively to premature release of
       the G2 transcriptional inhibition. This is likely to be a
       small effect as G2duration >> delta_t. */
    
    // update firing factor
    p->firingFactor = 1.0;

    // reupdate propensities based on new firing factor
    g->update->histone = TRUE;
    updatePropensities(c,p,g);
    
    // recalculate time step (p_s also recalculated).
    delta_t = gillespieTimeStep(p,g,&p_s);
    next_t = delta_t + r->t->el[step-1];
    r->new = fmod(next_t,3600*p->cellCycleDuration); 
    
    // proceed with gillespie reaction selection
    gillespieReaction = TRUE;
    
  }

  if (gillespieReaction == TRUE) {

    // Gillespie reaction.
    // ------------------------------

    // store new (stochastic) time
    r->t->el[step] = next_t;
    
    // choose reaction m from propensities based on scaled_r2
    r2 = runif(p->gsl_r);
    scaled_r2 = p_s*r2;

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

    // record time of each firing event
    if (g->update->transcribed == TRUE) {
      r->firing->el[p->reactCount] = TRUE;
      g->update->transcribed = FALSE;
    } else {
      r->firing->el[p->reactCount] = FALSE;
    }

    // update DNA replcation counter
    r->old = r->new;    
  }
  
  
return;
}
