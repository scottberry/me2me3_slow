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

void replicateDNA(chromatin *c, parameters *p, flags *update) {
  unsigned long pos;
  for (pos=0;pos<c->sites;pos+=2) {
    if(runif(p->gsl_r)<=0.5) {
      c->K27->el[pos] = me0;
      c->K27->el[pos+1] = me0;
    }
  }
  update->histone = TRUE;
  return;
}
