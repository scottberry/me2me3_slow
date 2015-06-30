#include "definitions.h"
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

void initialiseSilacLight(chromatin *c) {
  int i;
  for (i=0;i<c->sites;i++) {
    c->silac->el[i] = LIGHT;
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
        g->propensity->el[g->methylate_index->el[i]] = p->noisy_me0_me1 + p->me0_me1*(neighboursK27factor(c,p,i));
      } else if (c->K27->el[i] == me1) {
        g->propensity->el[g->methylate_index->el[i]] = p->noisy_me1_me2 + p->me1_me2*(neighboursK27factor(c,p,i));
      } else if (c->K27->el[i] == me2) {
        g->propensity->el[g->methylate_index->el[i]] = p->noisy_me2_me3 + p->me2_me3*(neighboursK27factor(c,p,i));
      } else {
        g->propensity->el[g->methylate_index->el[i]] = 0.0;
      }
    }
   
    // transcribeDNA
    g->propensity->el[g->transcribeDNA_index->el[0]] =
      p->activation*p->firingFactor*(p->firingRateMax + f_me2_me3*(p->firingRateMin - p->firingRateMax));
    
    /*
    g->propensity->el[g->transcribeDNA_index->el[0]] =
      p->activation * p->firingFactor *
      (p->firingRateMax - (p->firingRateMax - p->firingRateMin) * pow(f_me2_me3,p->firingHill) /
       (pow(f_me2_me3,p->firingHill) + pow(p->firingK,p->firingHill)));
    */
    g->update->histone = FALSE; // reset the flag
     
  }

  return;
}

/* Update the propensities based on change in system K27. */

void updatePropensitiesTranscriptionInhibit(chromatin *c, parameters *p, gillespie *g) {
  int i;
  double f_me2_me3, PRC2activity;
   
  if (g->update->histone==TRUE) {
    f_me2_me3 = fracControlRegion_me2me3(c);
    //f_me2_me3 = frac(c->K27,me2) + frac(c->K27,me3);

    if (c->transcribing == TRUE)
      PRC2activity = 1.0/p->PRC2inhibition;
    else
      PRC2activity = 1.0;
    
    for (i=0;i<c->sites;i++) {
      // fprintf(stderr,"i = %d, neighboursK27factor = %0.2f\n",i,neighboursK27factor(c,p,i));
      if (c->K27->el[i] == me0) { // methylate
        g->propensity->el[g->methylate_index->el[i]] = p->noisy_me0_me1 + p->me0_me1*(neighboursK27factor(c,p,i))*PRC2activity;
      } else if (c->K27->el[i] == me1) {
        g->propensity->el[g->methylate_index->el[i]] = p->noisy_me1_me2 + p->me1_me2*(neighboursK27factor(c,p,i))*PRC2activity;
      } else if (c->K27->el[i] == me2) {
        g->propensity->el[g->methylate_index->el[i]] = p->noisy_me2_me3 + p->me2_me3*(neighboursK27factor(c,p,i))*PRC2activity;
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

  return(delta_t); 
}
  
/* Single iteration of the "direct" Gillespie algorithm. Note that
   this varies somewhat from a pure stochastic simulation with
   exponentially distributed waiting times because there are events
   which occur at precise time points (DNA replication and release of
   G2 transcriptional inhibition). Some time steps proposed by the
   simulation algorithm are discarded as the system state has changed
   between the two simulation steps. Consequently the algorithm is
   non-Markovian. See \cite{Barrio:2006kb} for a description of the
   problems with delays in SSA simulations. The following
   implementation has been referred to as "delays as
   duration" approach in \cite{Barbuti:2009ih} */

void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r2=0.0, scaled_r2=0.0, next_t;
  long m, step;

  // local variable
  step = p->reactCount;
  
  // update propensities
  updatePropensities(c,p,g);
  
  // calculate time step
  // (p_s also returns propensity sum for use in Gillespie reaction selection).
  delta_t = gillespieTimeStep(p,g,&p_s);
  next_t = delta_t + r->t->el[step-1];

  /* Determine if the new time is after DNA replication or release of
     G2 inhibition. If so, interrupt Gillespie algorithm to replicate
     DNA or release G2 transcription inhibition (at a precise time).
     If not, proceed with Gillespie reaction selection. */

  if (next_t > g->t_nextRep) {

    // DNA replication.
    // ------------------------------

    // increment counters and record new system time
    p->cellCycleCount++;    
    r->t->el[step] = g->t_nextRep;

    // replicate DNA and instigate G2 transcriptional inhibition
    if (p->DNAreplication == TRUE)  {
      replicateDNA(c,p,g->update);
      if (p->G2duration > 0.0) p->firingFactor = 0.5;
    }
    
    // schedule next DNA replication
    g->t_nextRep += p->cellCycleDuration*3600;
    
    // fprintf(stderr,"replicate: t = %0.2f hours\n",r->t->el[step]/3600.0);
    
  } else if (next_t > g->t_nextEndG2 && p->G2duration > 0.0) {

    // Release G2 firing inhibition.
    // ------------------------------
    
    // update firing factor
    p->firingFactor = 1.0;

    // force propensity recalculation next step
    g->update->histone = TRUE;

    // increment system time
    r->t->el[step] = g->t_nextEndG2;

    // schedule next G2 inhibition end
    g->t_nextEndG2 += p->cellCycleDuration*3600;

    // fprintf(stderr,"release G2: t = %0.2f hours\n",r->t->el[step]/3600.0);
    
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

/* Gillespie algorithm with delays for each transcription event 
   Add variables:
   logical c.transcribing 
   double p.transcriptionDelay (duration of transcription event)
   double g.t_nextEndTranscription
*/

void gillespieStepTranscriptionDelays(chromatin *c, parameters *p, gillespie *g, record *r) {
  double delta_t, sum, p_s, r2=0.0, scaled_r2=0.0, next_t;
  long m, step;
  logical chooseGillespieReaction = TRUE;

  step = p->reactCount;
  
  // update propensities
  updatePropensitiesTranscriptionInhibit(c,p,g);
  
  // calculate time step
  // (p_s also returns propensity sum for use in Gillespie reaction selection).
  delta_t = gillespieTimeStep(p,g,&p_s);
  next_t = delta_t + r->t->el[step-1];

  /* Check delayed reactions.
     (priority in r->t update given to DNA replication if two delayed
     reactions occur simultaneously) */
  /* ---------------------------------------- */

  /* Terminate transcription */
  if (next_t > g->t_nextEndTranscription && c->transcribing == TRUE) {
    c->transcribing = FALSE; // update transcription status
    g->update->histone = TRUE; // force propensity recalculation
    r->t->el[step] = g->t_nextEndTranscription; // update system time
    chooseGillespieReaction = FALSE; // discard gillespie reaction
    //fprintf(stderr,"terminate: t = %0.2f sec\n",r->t->el[step]);
  }
  
  /*  Release G2 firing inhibition. */
  if (next_t > g->t_nextEndG2 && p->G2duration > 0.0) {
    p->firingFactor = 1.0; // update firing factor
    g->update->histone = TRUE; // force propensity recalculation
    r->t->el[step] = g->t_nextEndG2; // update system time
    g->t_nextEndG2 += p->cellCycleDuration*3600; // schedule next 
    chooseGillespieReaction = FALSE; // discard gillespie reaction
    // fprintf(stderr,"release G2: t = %0.2f hours\n",r->t->el[step]/3600.0);
  }

  /* DNA replication */
  if (next_t > g->t_nextRep) {
    p->cellCycleCount++; // increment counters 
    r->t->el[step] = g->t_nextRep; // update system time
    if (p->DNAreplication == TRUE)  {
      replicateDNA(c,p,g->update); // replicate DNA
      if (p->G2duration > 0.0) p->firingFactor = 0.5; // inhibit transcription
    }
    g->t_nextRep += p->cellCycleDuration*3600; // schedule next
    chooseGillespieReaction = FALSE; // discard gillespie reaction
    // fprintf(stderr,"replicate: t = %0.2f hours\n",r->t->el[step]/3600.0);
  }
  
  /* End delayed reactions */
  /* ---------------------------------------- */
  
  if (chooseGillespieReaction == TRUE) {

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

    // record each firing event
    if (g->update->transcribed == TRUE) {
      r->firing->el[step] = TRUE;
      // fprintf(stderr,"fire: t = %0.2f sec\n",r->t->el[step]);
      c->transcribing = TRUE;
      g->t_nextEndTranscription = r->t->el[step] + p->transcriptionDelay;
      g->update->transcribed = FALSE;
    } else {
      r->firing->el[step] = FALSE;
    }
    if (p->resultsTranscribing == TRUE) {
      r->transcribing->el[step] = c->transcribing; // update record
    }
  }
  
  return;
}
