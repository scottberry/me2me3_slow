
void fprint_t(char *fname, I_MAT *mat, int target) {
  FILE *fptr;
  long unsigned count, i, j;
  
  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols;i++) {
    count = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == target) {
	count++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)count/mat->rows);
  }
  fclose(fptr);
  return;
}

void fprint_t_nCycles(char *fname, I_MAT *mat, int target, record *r) {
  FILE *fptr;
  long unsigned count, i, j;
  
  // fprintf(stderr,"rows %ld, cols %ld\n",mat->rows,mat->cols);
  
  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols && r->t->el[i] < r->tMax;i++) {
    count = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == target) {
	count++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)count/(double)i);
  }
  fclose(fptr);
  return;
}

void fprint_firing_t(char *fname, record *r) {
  FILE *fptr;
  long unsigned i;

  fptr = fopen(fname,"w");
  for (i=0;i<r->firing->len;i++) {
    if (r->firing->el[i]==TRUE)
      fprintf(fptr,"%0.4f\n",r->t->el[i]);
  }
  fclose(fptr);
  return;
}

void fprint_firing_t_nCycles(char *fname, record *r) {
  FILE *fptr;
  long unsigned i;

  fptr = fopen(fname,"w");
  for (i=0;i<r->firing->len && r->t->el[i]<r->tMax;i++) {
    if (r->firing->el[i]==TRUE)
      fprintf(fptr,"%0.4f\n",r->t->el[i]);
  }
  fclose(fptr);
  return;
}

double tAverageGap(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double gapSum = 0;

  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
	sumM++;
    }
    gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(gapSum/r->t_out->el[p->samples-1]);
}

// average Gap only over the number of cell cycles simulated
double tAverageGap_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double gapSum = 0;

  for (t=1;t<r->state->cols && t<r->t_outLastSample;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
	sumM++;
    }
    gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(gapSum/r->t_out->el[r->t_outLastSample]);
}

double prob_me2_me3(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0;
  
  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
	sumM++;
    }
    // fprintf(stderr,"sumM = %ld\t",sumM);
    if (4*sumM > 3*c->sites)
      time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
  }

  // fprintf(stderr,"final time_in_M = %0.4f, total time = %0.4f\n",time_in_M,r->t_out->el[p->samples-1]);
  return(time_in_M/r->t_out->el[p->samples-1]);
}

/* Calculate probability for the last hour of each cell cycle only */
double prob_me2_me3_lastHour(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0;
  double time_total = 0;
  
  for (t=1;t<r->state->cols;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour before replication
      sumM = 0;
      for (pos=0;pos<r->state->rows;pos++) {
        if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
          sumM++;
      }
      // fprintf(stderr,"sumM = %ld\t",sumM);
      if (4*sumM > 3*c->sites) {
        time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
        // fprintf(stderr,"M state: t = %0.2f, sumM = %ld \n",r->t_out->el[t],sumM);
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  // fprintf(stderr,"final time_in_M = %0.4f, total time = %0.4f\n",time_in_M,r->t_out->el[p->samples-1]);
  return(time_in_M/time_total);
}

/* Calculate probability for the last hour of each cell cycle only
   (with cycles capped) */
double prob_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0;
  double time_total = 0;

  for (t=1;t<r->state->cols && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour before replication
      sumM = 0;
      for (pos=0;pos<r->state->rows;pos++) {
        if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
          sumM++;
      }
      // fprintf(stderr,"sumM = %ld\t",sumM);
      if (4*sumM > 3*c->sites) {
        time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
        // fprintf(stderr,"M state: t = %0.2f, sumM = %ld \n",r->t_out->el[t],sumM);
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  // fprintf(stderr,"final time_in_M = %0.4f, total time = %0.4f\n",time_in_M,r->t_out->el[p->samples-1]);
  return(time_in_M/time_total);
}

double prob_me0_me1(chromatin *c, parameters *p, record *r) {
  long sumU, t, pos;
  double time_in_U = 0;
  
  for (t=1;t<r->state->cols;t++) {
    sumU = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me0 || r->state->el[pos][t]==me1)
	sumU++;
    }
    if (4*sumU > 3*c->sites)
      time_in_U += r->t_out->el[t] - r->t_out->el[t-1];
  }

  return(time_in_U/r->t_out->el[p->samples-1]);
}

/* Calculate probability for the last hour of each cell cycle only */
double prob_me0_me1_lastHour(chromatin *c, parameters *p, record *r) {
  long sumU, t, pos;
  double time_in_U = 0;
  double time_total = 0;
  
  for (t=1;t<r->state->cols;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour before replication
      sumU = 0;
      for (pos=0;pos<r->state->rows;pos++) {
        if (r->state->el[pos][t]==me0 || r->state->el[pos][t]==me1)
          sumU++;
      }
      if (4*sumU > 3*c->sites) {
        time_in_U += r->t_out->el[t] - r->t_out->el[t-1];
        // fprintf(stderr,"U state: t = %0.2f, sumU = %ld \n",r->t_out->el[t],sumU);
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  
  return(time_in_U/time_total);
}

/* Calculate probability for the last hour of each cell cycle only
   (with cycles capped) */
double prob_me0_me1_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumU, t, pos;
  double time_in_U = 0;
  double time_total = 0;

  for (t=1;t<r->state->cols && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour before replication
      sumU = 0;
      for (pos=0;pos<r->state->rows;pos++) {
        if (r->state->el[pos][t]==me0 || r->state->el[pos][t]==me1)
          sumU++;
      }
      if (4*sumU > 3*c->sites) {
        time_in_U += r->t_out->el[t] - r->t_out->el[t-1];
        // fprintf(stderr,"U state: t = %0.2f, sumU = %ld \n",r->t_out->el[t],sumU);
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  
  return(time_in_U/time_total);
}

double tAverage_me2_me3(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double Mavg = 0;

  for (t=1;t<r->state->cols;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
	sumM++;
    }
    Mavg += (double)sumM*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(Mavg/r->t_out->el[p->samples-1]);
}

double tAverage_me2_me3_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double Mavg = 0;

  for (t=1;t<r->state->cols && t<r->t_outLastSample;t++) {
    sumM = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3)
	sumM++;
    }
    Mavg += (double)sumM*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(Mavg/r->t_out->el[r->t_outLastSample-1]);
}

/* Calculate the average lifetime of a state */
unsigned long numberHistoneStateFlips(record *r) {
  signed char newState = 0, oldState = 0;
  unsigned long t, pos, flips=0, m=0, u=0;;
  
  for (t=0;t<r->state->cols && t<r->t_outLastSample;t++) {
    
    oldState = newState;
    
    m = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3) m++;
    }
    u = r->state->rows - m;

    if (oldState == 1) { // if previously U
      if (m >= 3*u) {
	newState = -1;
	flips++;
      }
    } else if (oldState == -1) { // if previously M
      if (u >= 3*m) {
	newState = 1;
	flips++;
      }
    } else if (oldState == 0) { // first time-point
      if (u >= 3*m)
	newState = 1;
      else
	newState = -1;
    }
  }
 
  return(flips);
}

double firstPassageTime(record *r, signed char *initial) {
  long unsigned m=0, u=0, pos,t=0;
  
  /* find initial state */
  for (pos=0;pos<r->state->rows;pos++) {
    if (r->state->el[pos][0]==me2 || r->state->el[pos][0]==me3) m++;
  }
  u = r->state->rows - m;
  
  if (m > u)
    *initial = -1;
  else
    *initial = 1;
    
  while ( t < r->state->cols && t < r->t_outLastSample &&
	  ((*initial == -1 && 3*m > u) || (*initial == 1 && 3*u > m))) {
    m = 0;
    u = 0;
    for (pos=0;pos<r->state->rows;pos++) {
      if (r->state->el[pos][t]==me2 || r->state->el[pos][t]==me3) m++;
    }
    u = r->state->rows - m;
    //fprintf(stderr,"t = %0.2f, m %ld u %ld \n",r->t_out->el[t],m,u);
    t++;
  }
  return(r->t_out->el[t-1]);
}

/* write a log file */
int writelog(FILE *fptr, chromatin *c, parameters *p, record *r) {
  time_t curtime;
  struct tm *loctime;

  curtime = time (NULL);
  loctime = localtime (&curtime);

  fputs (asctime (loctime), fptr);
#ifdef __APPLE__
  fprintf(fptr,"Operating system: Mac OS\n");
#else
  fprintf(fptr,"Operating system: Unix\n");
#endif
  fprintf(fptr,"sites: %ld\n", c->sites);
  fprintf(fptr,"loci: %ld\n", p->loci);
  fprintf(fptr,"maxReact: %ld\n", p->maxReact);
  fprintf(fptr,"samples: %ld\n", p->samples);
  fprintf(fptr,"cellCycleDuration: %0.2f hours\n", p->cellCycleDuration);
  fprintf(fptr,"G2duration: %0.2f hours\n", p->G2duration);
  fprintf(fptr,"cellCycles: %d\n", p->cellCycles);
  fprintf(fptr,"me0_me1: %0.6f\n", p->me0_me1);
  fprintf(fptr,"me1_me2: %0.6f\n", p->me1_me2);
  fprintf(fptr,"me2_me3: %0.6f\n", p->me2_me3);
  fprintf(fptr,"me2factor: %0.6f\n", p->me2factor);
  fprintf(fptr,"me3factor: %0.6f\n", p->me3factor);
  fprintf(fptr,"firingRateMax: %0.6f\n", p->firingRateMax);
  fprintf(fptr,"firingRateMin: %0.6f\n", p->firingRateMin);

  fprintf(fptr,"transcription_demethylate: %0.6f\n", p->transcription_demethylate);

  return(1);
}
