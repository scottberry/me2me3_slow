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

void replicateDNA(chromatin *c, parameters *p, flags *update) {
  unsigned long pos;
  for (pos=0;pos<c->sites;pos+=2) {
    // note +=2 in loop implies only check once per nucleosome
    if(runif(p->gsl_r)<=0.5) {
      // replace with unmodified nucleosome
      c->K27->el[pos] = me0;
      c->K27->el[pos+1] = me0;
    }
  }
  update->histone = TRUE;
  return;
}

