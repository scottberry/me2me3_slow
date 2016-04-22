#include "definitions.h"
/* 
   Implementation of the direct Gillespie algorithm for simulation
   of models of chromatin-based epigenetics.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

/* Allocate memory for chromatin, records, gillespie objects, 
   silac states and histone turnover records */
void allocateGillespieMemory(chromatin *c, parameters *p, gillespie *g, record *r) {

  // Chromatin
  c->K27 = i_vec_get( c->sites );

  // Records
  r->t = d_vec_get(p->maxReact + 1);
  r->firing = i_vec_get(p->maxReact + 1);
  r->t_out = d_vec_get(p->samples);
  r->K27 = i_mat_get(c->sites,p->samples);

  // Silac
  if (p->silacExperiment == TRUE)
    c->silac = i_vec_get( c->sites );

  // Gillespie objects
  g->methylate_index = i_vec_get( c->sites );
  g->demethylate_index = i_vec_get( c->sites );
  g->transcribeDNA_index = i_vec_get( 1 );
  if (p->stochasticAlpha == TRUE) {
    g->propensity = d_vec_get( 2*c->sites + 5 );
    g->increaseProtein_index = i_vec_get( 1 );
    g->decreaseProtein_index = i_vec_get( 1 );
    g->increaseRNA_index = i_vec_get( 1 );
    g->decreaseRNA_index = i_vec_get( 1 );
    r->transFactorProtein = i_vec_get( p->samples );
    r->transFactorRNA = i_vec_get( p->samples );
    r->alpha = d_vec_get( p->samples );
  } else if (p->burstyFiring == TRUE) {
    g->propensity = d_vec_get( 2*c->sites + 3 );    
    g->activatePromoter_index = i_vec_get( 1 );
    g->deactivatePromoter_index = i_vec_get( 1 );
    r->promoterON = i_vec_get(p->samples);
  } else {
    g->propensity = d_vec_get( 2*c->sites + 1 );
  }

  g->doReaction = malloc(g->propensity->len*sizeof( func_ptr_t ) );
  g->doReactionParam = i_vec_get( g->propensity->len );
  g->update = malloc(sizeof( flags ) );

  // Histone turnover
  if (p->checkHistoneTurnover == TRUE) {
    c->turnover = i_vec_get(4);
    r->turnover = d_mat_get(p->loci,4);
    c->variant = i_vec_get(c->sites);
    r->variant = i_mat_get(c->sites,p->samples);
  }
  
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
  
  // Silac
  if (p->silacExperiment == TRUE)
    i_vec_free(c->silac);

  // Gillespie objects
  i_vec_free(g->methylate_index);
  i_vec_free(g->demethylate_index);
  i_vec_free(g->transcribeDNA_index);
  if (p->stochasticAlpha == TRUE) {
    i_vec_free(g->increaseProtein_index);
    i_vec_free(g->decreaseProtein_index);
    i_vec_free(g->increaseRNA_index);
    i_vec_free(g->decreaseRNA_index);
    i_vec_free(r->transFactorProtein);
    i_vec_free(r->transFactorRNA);
    d_vec_free(r->alpha);
  } else if (p->burstyFiring) {
    i_vec_free(g->activatePromoter_index);
    i_vec_free(g->deactivatePromoter_index);
    i_vec_free(r->promoterON);
  }

  d_vec_free(g->propensity);
  free(g->doReaction);
  i_vec_free(g->doReactionParam);
  free(g->update);
  rfree(p);

  // Histone turnover
  if (p->checkHistoneTurnover == TRUE) {
    i_vec_free(c->turnover);  
    d_mat_free(r->turnover);
    i_vec_free(c->variant);
    i_mat_free(r->variant);
  }
    
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

void initialiseSilacLight(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->silac->el[i] = LIGHT;
  }
  return;
}

void initialiseH3_1(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->variant->el[i] = H3_1;
  }
  return;
}

void initialisePromoterOFF(chromatin *c) {
  c->promoterON = FALSE;
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
  for (i=c->sites;i<2*c->sites;i++) { // demethylate
    g->doReaction[i] = demethylate;
    g->doReactionParam->el[i] = i-c->sites;
    g->demethylate_index->el[i-c->sites] = i;
  }
  g->doReaction[2*c->sites] = transcribeDNA; // transcribeDNA
  g->doReactionParam->el[2*c->sites] = 0;
  g->transcribeDNA_index->el[0] = 2*c->sites;
  
  if (g->test == TRUE) {
    fprintf(g->test_fptr,"Methylate index:\n");
    i_vec_print(g->test_fptr,g->methylate_index);
    fprintf(g->test_fptr,"Demethylate index:\n");
    i_vec_print(g->test_fptr,g->demethylate_index);
    fprintf(g->test_fptr,"Transcribe index:\n");
    i_vec_print(g->test_fptr,g->transcribeDNA_index); 
    fprintf(g->test_fptr,"doReactionParam:\n");
    i_vec_print(g->test_fptr,g->doReactionParam); 
  }
  return;
}

/* Additional function to add trans-factor simulations to propensity
   matrix */
void initialiseGillespieFunctionsTransFactor(chromatin *c, gillespie *g) {

  g->doReaction[2*c->sites+1] = decreaseProtein; // decreaseP
  g->doReactionParam->el[2*c->sites+1] = 0;
  g->decreaseProtein_index->el[0] = 2*c->sites+1;

  g->doReaction[2*c->sites+2] = increaseProtein; // increaseP
  g->doReactionParam->el[2*c->sites+2] = 0;
  g->increaseProtein_index->el[0] = 2*c->sites+2;

  g->doReaction[2*c->sites+3] = decreaseRNA; // decreaseR
  g->doReactionParam->el[2*c->sites+3] = 0;
  g->decreaseRNA_index->el[0] = 2*c->sites+3;

  g->doReaction[2*c->sites+4] = increaseRNA; // increaseR
  g->doReactionParam->el[2*c->sites+4] = 0;
  g->increaseRNA_index->el[0] = 2*c->sites+4;

  return;
}

/* Additional function to add promoter status */
void initialiseGillespieFunctionsBurstyFiring(chromatin *c, gillespie *g) {

  g->doReaction[2*c->sites+1] = activatePromoter;
  g->doReactionParam->el[2*c->sites+1] = 0;
  g->activatePromoter_index->el[0] = 2*c->sites+1;

  g->doReaction[2*c->sites+2] = deactivatePromoter;
  g->doReactionParam->el[2*c->sites+2] = 0;
  g->deactivatePromoter_index->el[0] = 2*c->sites+2;

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
    if (c->K27->el[i]==me2 || c->K27->el[i]==me3) count++;
  }
  f = (double)count/c->controlSites;
  return(f);
}

/* Return the enzymatic activity factor for PRC2 associated with the
   recruitment potential of an me2 or me3 histone mark */ 
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

/* Find histone mods on neighbouring nucleosome and return the
   relative activity of PRC2 for a particular histone position 'pos' */
double neighboursK27factor(chromatin *c, parameters *p, int pos) {
  double n = 0;

  if (pos % 2 == 0) { // (even) left histone tail
    n += enzymaticFactor(c,p,pos+1); // other tail on same nucleosome
    if (pos != 0) { // check not left-most nucleosome in domain
      n += enzymaticFactor(c,p,pos-1); // left tail on left neighbour nucleosome
      n += enzymaticFactor(c,p,pos-2); // right tail on left neighbour nucleosome
    }
    if (pos < c->sites - 3 ) { // check not right-most nucleosome in domain
      n += enzymaticFactor(c,p,pos+2); // left tail on right neighbour nucleosome
      n += enzymaticFactor(c,p,pos+3); // right tail on right neighbour nucleosome    
    }
  } else { // (odd) right histone tail
    n += enzymaticFactor(c,p,pos-1); // other tail on same nucleosome
    if (pos != 1) { // check not left-most nucleosome in domain
      n += enzymaticFactor(c,p,pos-2); // left tail on left neighbour nucleosome
      n += enzymaticFactor(c,p,pos-3); // right tail on left neighbour nucleosome
    }
    if (pos < c->sites - 2) { // check not right-most nucleosome in domain
      n += enzymaticFactor(c,p,pos+1); // left tail on right neighbour nucleosome
      n += enzymaticFactor(c,p,pos+2); // right tail on right neighbour nucleosome    
    }
  }
  return(n);
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

double k_on(parameters *p, double f_me2_me3) {
  double f;
  if (f_me2_me3 < p->firingThreshold) {
    f = p->k_onMax - (f_me2_me3 * (p->k_onMax - p->k_onMin) / p->firingThreshold);
  } else {
    f = p->k_onMin;
  }
  return(f);
}

/* Update the propensities based on change in system */
void updatePropensities(chromatin *c, parameters *p, gillespie *g) {
  int i;
  double f_me2_me3;
   
  if (g->update->histone==TRUE) {
    f_me2_me3 = fracControlRegion_me2me3(c);

    if (p->stochasticAlpha==TRUE)
      p->alpha = (double)p->transFactorProtein/1000.0;

    // chromatin
    for (i=0;i<c->sites;i++) {
      // fprintf(stderr,"i = %d, neighboursK27factor = %0.2f\n",i,neighboursK27factor(c,p,i));
      if (c->K27->el[i] == me0) { // methylate
        g->propensity->el[g->methylate_index->el[i]] = p->beta*(p->noisy_me0_me1 + p->me0_me1*(neighboursK27factor(c,p,i)));
        g->propensity->el[g->demethylate_index->el[i]] = 0.0;
      } else if (c->K27->el[i] == me1) {
        g->propensity->el[g->methylate_index->el[i]] = p->beta*(p->noisy_me1_me2 + p->me1_me2*(neighboursK27factor(c,p,i)));
        g->propensity->el[g->demethylate_index->el[i]] = p->noisy_demethylate;
      } else if (c->K27->el[i] == me2) {
        g->propensity->el[g->methylate_index->el[i]] = p->beta*(p->noisy_me2_me3 + p->me2_me3*(neighboursK27factor(c,p,i)));
        g->propensity->el[g->demethylate_index->el[i]] = p->noisy_demethylate;
      } else {
        g->propensity->el[g->methylate_index->el[i]] = 0.0;
        g->propensity->el[g->demethylate_index->el[i]] = p->noisy_demethylate;
      }
    }

    // transcription
    if (p->burstyFiring == FALSE ) { // if not bursty
      g->propensity->el[g->transcribeDNA_index->el[0]] = p->firingFactor * p->alpha * firingRate(p,f_me2_me3);
            
    } else { // if bursty
      if (c->promoterON == TRUE) {
        g->propensity->el[g->deactivatePromoter_index->el[0]] = p->k_off;
        g->propensity->el[g->transcribeDNA_index->el[0]] = p->constFiring;
      } else {
        g->propensity->el[g->activatePromoter_index->el[0]] = p->alpha * k_on(p,f_me2_me3);
        g->propensity->el[g->transcribeDNA_index->el[0]] = 0.0;
      }
    }
    if (p->capFiring == TRUE) {
      // cap the firing rate 
      if (g->propensity->el[g->transcribeDNA_index->el[0]] > p->firingCap)
        g->propensity->el[g->transcribeDNA_index->el[0]] = p->firingCap;
    }
    g->update->histone = FALSE; // reset the flag
  }
  return;
}

/* Update the propensities based on change in system, when including
   variable activation through trans-factor simulations */
void updatePropensitiesTransFactor(parameters *p, gillespie *g) {

  // update trans-regulator concentration
  g->propensity->el[g->decreaseProtein_index->el[0]] = p->gamma_p * (double)p->transFactorProtein;
  g->propensity->el[g->increaseProtein_index->el[0]] = p->k_p * (double)p->transFactorRNA;

  g->propensity->el[g->decreaseRNA_index->el[0]] = p->gamma_r * (double)p->transFactorRNA;
  g->propensity->el[g->increaseRNA_index->el[0]] = p->k_r;
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
  
/* Single iteration of the "direct" Gillespie algorithm. Note that
   this varies somewhat from a pure stochastic simulation with
   exponentially distributed waiting times because there are events
   which occur at precise time points (i.e. DNA replication). 
   Some time steps proposed by the simulation algorithm are discarded
   as the system state has changed between the two simulation steps.
   See \cite{Barrio:2006kb} for a description of the problems with
   delays in SSA simulations. The following implementation has been
   referred to as "delays as duration" approach in
   \cite{Barbuti:2009ih} or altenatively \cite{Bratsun:2005cs}. */

void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r2=0.0, scaled_r2=0.0, next_t;
  long m, step;

  // local variable
  step = p->reactCount;
  
  // update propensities
  updatePropensities(c,p,g);
  if (p->stochasticAlpha == TRUE)
    updatePropensitiesTransFactor(p,g);

  /*  for (i=0;i<g->propensity->len;i++)
    fprintf(stderr,"%ld,%0.10f\n",i,g->propensity->el[i]);
  */
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

