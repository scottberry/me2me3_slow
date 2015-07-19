#include "definitions.h"
/* 
   Individual state modifications for Gillsepie algorithm simulations.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

void methylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->K27->el[pos] == me0) {
    c->K27->el[pos] = me1;
  } else if (c->K27->el[pos] == me1) {
    c->K27->el[pos] = me2;
  } else if (c->K27->el[pos] == me2) {
    c->K27->el[pos] = me3;
  }
  update->histone = TRUE;
  return;
}

void demethylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->K27->el[pos] == me3) {
    c->K27->el[pos] = me2;
  } else if (c->K27->el[pos] == me2) {
    c->K27->el[pos] = me1;
  } else if (c->K27->el[pos] == me1) {
    c->K27->el[pos] = me0;
  }
  update->histone = TRUE;
  return;
}
/*
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
*/

void transcribeDNA(chromatin *c, parameters *p, flags *update, int pos) {
  unsigned long i;
  double rand;

  for (i=0;i<c->sites;i++) {
    rand = runif(p->gsl_r);
    if (rand < p->transcription_demethylate) {
      demethylate(c,p,update,i);
    } else if (rand < p->transcription_demethylate + p->transcription_turnover) {
      if (i % 2 == 0) { // even histone
        if (p->checkHistoneTurnover==TRUE) {
          if (c->K27->el[i]==me0) c->turnover->el[0]++;
          if (c->K27->el[i+1]==me0) c->turnover->el[0]++;
          if (c->K27->el[i]==me1) c->turnover->el[1]++;
          if (c->K27->el[i+1]==me1) c->turnover->el[1]++;
          if (c->K27->el[i]==me2) c->turnover->el[2]++;
          if (c->K27->el[i+1]==me2) c->turnover->el[2]++;
          if (c->K27->el[i]==me3) c->turnover->el[3]++;
          if (c->K27->el[i+1]==me3) c->turnover->el[3]++;
        }
        c->K27->el[i] = me0;
        c->K27->el[i+1] = me0;
        update->histone = TRUE;
      } else { // odd histone
        if (p->checkHistoneTurnover==TRUE) {
          if (c->K27->el[i]==me0) c->turnover->el[0]++;
          if (c->K27->el[i-1]==me0) c->turnover->el[0]++;
          if (c->K27->el[i]==me1) c->turnover->el[1]++;
          if (c->K27->el[i-1]==me1) c->turnover->el[1]++;
          if (c->K27->el[i]==me2) c->turnover->el[2]++;
          if (c->K27->el[i-1]==me2) c->turnover->el[2]++;
          if (c->K27->el[i]==me3) c->turnover->el[3]++;
          if (c->K27->el[i-1]==me3) c->turnover->el[3]++;
        }
        c->K27->el[i] = me0;
        c->K27->el[i-1] = me0;
        update->histone = TRUE;
      }
    }
  }
  update->transcribed = TRUE;
  return;
}

void replicateDNA(chromatin *c, parameters *p, flags *update) {
  unsigned long pos;
  for (pos=0;pos<c->sites;pos+=2) {
    if(runif(p->gsl_r)<=0.5) {
      c->K27->el[pos] = me0;
      c->K27->el[pos+1] = me0;
      if (p->silacExperiment == TRUE) {
        if (p->silacLabel == LIGHT) {
          c->silac->el[pos] = LIGHT;
          c->silac->el[pos+1] = LIGHT;
        } else if (p->silacLabel == HEAVY) {
          c->silac->el[pos] = HEAVY;
          c->silac->el[pos+1] = HEAVY;
        } else {
          c->silac->el[pos] = UNLABELLED;
          c->silac->el[pos+1] = UNLABELLED;
        }
      }
    }
  }
  update->histone = TRUE;
  return;
}

void decreaseRNA(chromatin *c, parameters *p, flags *update, int pos) {
  p->transFactorRNA--;
  update->histone = TRUE;
  return;  
}

void increaseRNA(chromatin *c, parameters *p, flags *update, int pos) {
  p->transFactorRNA++;
  update->histone = TRUE;
  return;  
}

void decreaseProtein(chromatin *c, parameters *p, flags *update, int pos) {
  p->transFactorProtein--;
  update->histone = TRUE;
  return;  
}

void increaseProtein(chromatin *c, parameters *p, flags *update, int pos) {
  p->transFactorProtein++;
  update->histone = TRUE;
  return;  
}
