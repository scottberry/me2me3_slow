#include "definitions.h"
/* 
   Implementation of the direct Gillespie algorithm for simulation
   of models of chromatin-based epigenetics.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

/* Allocate memory for chromatin, records, gillespie objects */
void allocateGillespieMemory(chromatin *c, parameters *p, gillespie *g, record *r) {

  // Chromatin
  c->K27 = i_vec_get( c->sites );

  // Records
  r->t = d_vec_get(p->maxReact + 1);
  r->firing = i_vec_get(p->maxReact + 1);
  r->t_out = d_vec_get(p->samples);
  r->K27 = i_mat_get(c->sites,p->samples);

  // Gillespie objects
  g->methylate_index = i_vec_get( c->sites );
  g->transcribeDNA_index = i_vec_get( 1 );
  g->propensity = d_vec_get( c->sites + 1 );

  g->doReaction = malloc(g->propensity->len*sizeof( func_ptr_t ) );
  g->doReactionParam = i_vec_get( g->propensity->len );
  g->update = malloc(sizeof( flags ) );
  
  return;
}

/* Free memory associated with Gillespie simulation */
void freeGillespieMemory(chromatin *c, parameters *p, gillespie *g, record *r) {

  // Chromatin
  i_vec_free(c->K27);

  // Records
  i_vec_free(r->firing);
  i_mat_free(r->K27);
  d_vec_free(r->t);
  d_vec_free(r->t_out);
  
  // Gillespie objects
  i_vec_free(g->methylate_index);
  i_vec_free(g->transcribeDNA_index);

  d_vec_free(g->propensity);
  free(g->doReaction);
  i_vec_free(g->doReactionParam);
  free(g->update);
  rfree(p);
    
  return;
}

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

void initialiseMixed(chromatin *c, parameters *p) {
  int i;

  for (i=0;i<c->sites;i++) {
    if (runif(p->gsl_r) <= 0.5)  {
      c->K27->el[i] = me0;
    } else {
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

/* Initialise Gillespie functions such as propensity vectors and
   indices. Called only once to initialise indices and function pointers */
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
  
  if (g->test == TRUE) {
    fprintf(g->test_fptr,"Methylate index:\n");
    i_vec_print(g->test_fptr,g->methylate_index);
    fprintf(g->test_fptr,"Transcribe index:\n");
    i_vec_print(g->test_fptr,g->transcribeDNA_index); 
    fprintf(g->test_fptr,"doReactionParam:\n");
    i_vec_print(g->test_fptr,g->doReactionParam); 
  }
  return;
}

/* Calculate the fraction of 'target' reactants */
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

/* Calculate the fraction of 'target' reactants for a control
   region of size specified in c->controlSites */
double fracControlRegion_me2me3(chromatin *c) {
  unsigned long int count = 0;
  double f;
  int i;
  for (i=0;i<c->controlSites;i++) {
    if (c->K27->el[i]==me3) count++;
  }
  f = (double)count/c->controlSites;
  return(f);
}

/* Transcriptional firing-rate function: linear between f_max and
   f_min, which occurs at f_me2_me3 = firingThreshold. Above
   firingThreshold f = f_min */
double firingRate(parameters *p, double f_me2_me3) {
  double f;
  if (f_me2_me3 < p->firingThreshold) {
    f = p->firingRateMax - (f_me2_me3 * (p->firingRateMax - p->firingRateMin) / p->firingThreshold);
  } else {
    f = p->firingRateMin;
  }
  return(f);
}

/* Update the propensities based on change in system */
void updatePropensities(chromatin *c, parameters *p, gillespie *g) {
  int i;
  double f_me2_me3;
   
  if (g->update->histone==TRUE) {
    f_me2_me3 = fracControlRegion_me2me3(c);
    
    for (i=0;i<c->sites;i++) {
      if (c->K27->el[i] == me0) { // methylate
        g->propensity->el[g->methylate_index->el[i]] = p->noisy_me2_me3 + p->me2_me3 * f_me2_me3;
      } else {
        g->propensity->el[g->methylate_index->el[i]] = 0.0;
      }
    }
   
    // transcribeDNA
    g->propensity->el[g->transcribeDNA_index->el[0]] = firingRate(p,f_me2_me3);

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

  return(delta_t); 
}
  
/* Single iteration of the "direct" Gillespie algorithm */
void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r2=0.0, scaled_r2=0.0, next_t;
  long m, step;

  // local variable
  step = p->reactCount;
  
  // update propensities
  updatePropensities(c,p,g);

  // Calculate time step
  // (propensity sum returned through argument for use in Gillespie
  // reaction selection).
  delta_t = gillespieTimeStep(p,g,&p_s);
  next_t = delta_t + r->t->el[step-1];

  
  /* Determine if the new time is after DNA replication. If so,
     interrupt Gillespie algorithm to replicate DNA (at a precise
     time). If not, proceed with Gillespie reaction selection. */
  
  if (next_t > g->t_nextRep) {

    // DNA replication.
    // ------------------------------

    // increment counters and record new system time
    p->cellCycleCount++;
    r->t->el[step] = g->t_nextRep;

    // Update Silac Label
    if (p->silacExperiment == TRUE) {
      if (p->cellCycleCount > p->silacLightCycles)
        p->silacLabel = HEAVY;
      if (p->cellCycleCount > p->silacLightCycles + p->silacHeavyCycles)
        p->silacLabel = UNLABELLED;
    }
    
    // replicate DNA
    if (p->DNAreplication == TRUE)  {
      replicateDNA(c,p,g->update);
    }
    
    // schedule next DNA replication
    g->t_nextRep += p->cellCycleDuration*3600;
    
    // fprintf(stderr,"replicate: t = %0.2f hours\n",r->t->el[step]/3600.0);
    
  } else {
    
    // Gillespie reaction selection.
    // ------------------------------
    
    // store new (stochastically chosen) time
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
  }
  
  return;
}

