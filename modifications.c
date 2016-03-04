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
    } else if (rand < p->transcription_demethylate + p->transcription_turnover) {
      if (i % 2 == 0) {
        /*---EVEN---*/
        if (p->checkHistoneTurnover==TRUE) {
          // record identity of mark removed by turnover event
          if (c->K27->el[i]==me0) c->turnover->el[0]++;
          if (c->K27->el[i+1]==me0) c->turnover->el[0]++;
          if (c->K27->el[i]==me1) c->turnover->el[1]++;
          if (c->K27->el[i+1]==me1) c->turnover->el[1]++;
          if (c->K27->el[i]==me2) c->turnover->el[2]++;
          if (c->K27->el[i+1]==me2) c->turnover->el[2]++;
          if (c->K27->el[i]==me3) c->turnover->el[3]++;
          if (c->K27->el[i+1]==me3) c->turnover->el[3]++;
          // replace with H3.3 rather than H3.1 for
          // transcription-coupled histone exchange
          c->variant->el[i] = H3_3;
          c->variant->el[i+1] = H3_3;
        }
        // remove histone modifications
        c->K27->el[i] = me0;
        c->K27->el[i+1] = me0;
        
        if (p->silacExperiment == TRUE) {
          // ensure new histone has correct silac label
          if (p->silacLabel == LIGHT) {
            c->silac->el[i] = LIGHT;
            c->silac->el[i+1] = LIGHT;
          } else if (p->silacLabel == HEAVY) {
            c->silac->el[i] = HEAVY;
            c->silac->el[i+1] = HEAVY;
          } else {
            c->silac->el[i] = UNLABELLED;
            c->silac->el[i+1] = UNLABELLED;
          }
        }
        update->histone = TRUE;
      } else {
        /*---ODD---*/
        if (p->checkHistoneTurnover==TRUE) {
          // record identity of mark removed by turnover event
          if (c->K27->el[i]==me0) c->turnover->el[0]++;
          if (c->K27->el[i-1]==me0) c->turnover->el[0]++;
          if (c->K27->el[i]==me1) c->turnover->el[1]++;
          if (c->K27->el[i-1]==me1) c->turnover->el[1]++;
          if (c->K27->el[i]==me2) c->turnover->el[2]++;
          if (c->K27->el[i-1]==me2) c->turnover->el[2]++;
          if (c->K27->el[i]==me3) c->turnover->el[3]++;
          if (c->K27->el[i-1]==me3) c->turnover->el[3]++;
          // replace with H3.3 rather than H3.1 for
          // transcription-coupled histone exchange
          c->variant->el[i] = H3_3;
          c->variant->el[i-1] = H3_3;
        }
        // remove histone modifications
        c->K27->el[i] = me0;
        c->K27->el[i-1] = me0;
        
        if (p->silacExperiment == TRUE) {
          // ensure new histone has current silac label
          if (p->silacLabel == LIGHT) {
            c->silac->el[i] = LIGHT;
            c->silac->el[i-1] = LIGHT;
          } else if (p->silacLabel == HEAVY) {
            c->silac->el[i] = HEAVY;
            c->silac->el[i-1] = HEAVY;
          } else {
            c->silac->el[i] = UNLABELLED;
            c->silac->el[i-1] = UNLABELLED;
          }
        }
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
    // note +=2 in loop implies only check once per nucleosome
    if(runif(p->gsl_r)<=0.5) {
      // replace with unmodified nucleosome
      c->K27->el[pos] = me0;
      c->K27->el[pos+1] = me0;
      if (p->checkHistoneTurnover==TRUE) {
        // replace with H3.1 rather than H3.3 for
        // replication-specific incorporation
        c->variant->el[pos] = H3_1;
        c->variant->el[pos+1] = H3_1;
      }
      if (p->silacExperiment == TRUE) {
        // ensure new histone has current silac label
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

void activatePromoter(chromatin *c, parameters *p, flags *update, int pos) {
  c->promoterON = TRUE;
  update->histone = TRUE;
  return;  
}

void deactivatePromoter(chromatin *c, parameters *p, flags *update, int pos) {
  c->promoterON = FALSE;
  update->histone = TRUE;
  return;  
}
