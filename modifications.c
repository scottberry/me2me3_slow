/* Individual reactions */
void methylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == me0) {
    c->state->el[pos] = me1;
  } else if (c->state->el[pos] == me1) {
    c->state->el[pos] = me2;
  } else if (c->state->el[pos] == me2) {
    c->state->el[pos] = me3;
  }
  update->histone = TRUE;
  return;
}

void demethylate(chromatin *c, parameters *p, flags *update, int pos) {
  if (c->state->el[pos] == me3) {
    c->state->el[pos] = me2;
  } else if (c->state->el[pos] == me2) {
    c->state->el[pos] = me1;
  } else if (c->state->el[pos] == me1) {
    c->state->el[pos] = me0;
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
  for (pos=0;pos<c->sites;pos++) {
    if(runif(p->gsl_r)<=0.5) {
      c->state->el[pos] = me0;
    }
  }
  update->protein = TRUE;
  update->histone = TRUE;
  return;
}
