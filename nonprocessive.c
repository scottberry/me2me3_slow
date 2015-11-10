#include "definitions.h"
/* 
   Non-processive methylation and demethylation reactions
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

/* Non-processive methylation */
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

/* Non-processive demethylation */
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

