#include "definitions.h"
/* 
   Individual state modifications for Gillsepie algorithm simulations.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

void transcribeDNA(chromatin *c, parameters *p, flags *update, int pos) {
  unsigned long i;
  double rand;

  for (i=0;i<c->sites;i++) {
    rand = runif(p->gsl_r);
    if (rand < p->transcription_demethylate) {
      demethylate(c,p,update,i);
    }
  }
  update->transcribed = TRUE;
  return;
}

void methylate(chromatin *c, parameters *p, flags *update, int pos) {
  c->K27->el[pos] = me3;
  update->histone = TRUE;
  return;
}

void demethylate(chromatin *c, parameters *p, flags *update, int pos) {
  c->K27->el[pos] = me0;
  update->histone = TRUE;
  return;
}


