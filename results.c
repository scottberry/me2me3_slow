/* 
   Output functions for calulating results and writing to file.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
 */

/* Print simulation time at sample points. Length depends directly on
   r.tMax, which is determined by p.cellCycles. */

void fprint_t_out_nCycles(char *fname, record *r) {
  FILE *fptr;
  long unsigned count, i, j;
  
  fptr = fopen(fname,"w");
  for (i = 0;i < r->t_out->len && i < r->t_outLastSample;i++) {
    fprintf(fptr,"%0.4f\n",r->t_out->el[i]);
    // fprintf(stderr,"%0.4f %0.4f\n",r->tMax, r->t_out->el[i]);
  }
  fclose(fptr);
  return;
}

/* Average results over time and print a time-dependent results 
   vector. Length depends directly on r.tMax, which is 
   determined by p.cellCycles. */

void fprint_t_nCycles(char *fname, I_MAT *mat, int target, record *r) {
  FILE *fptr;
  long unsigned count, i, j;
  
  // fprintf(stderr,"rows %ld, cols %ld\n",mat->rows,mat->cols);
  
  fptr = fopen(fname,"w");
  for (i = 0;i < mat->cols && i < r->t_outLastSample;i++) {
    count = 0;
    for (j=0;j < mat->rows;j++) {
      if (mat->el[j][i] == target) {
	count++;
      }
    }
    fprintf(fptr,"%0.4f\n",(double)count/(double)j);
  }
  fclose(fptr);
  return;
}

/* Print absolute time of each firing event.
   Length depends directly on r.tMax, which is 
   determined by p.cellCycles. */

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

/* Calculate the "Gap" parameter, as defined in Dodd et al. 2007:
   |M-A|/(M+A). Average over time for a single locus. Note that
   function evaluation accounts for non-constant time-step using 
   a discrete "rectangular" integration approach. Average between
   samples 1 and r.t_outLastSample. Determined by whether 
   p.maxReact or steps required to simulate p.cellCycles is 
   greater. */

double tAverageGap_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double gapSum = 0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
	sumM++;
    }
    gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(gapSum/r->t_out->el[r->t_outLastSample]);
}

/* See above - but for last hour of each cell cycle only */

double tAverageGap_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double gapSum = 0.0, time_total = 0.0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
          sumM++;
      }
      gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }

  if (time_total == 0.0)
    fprintf(stderr,"Error: tAverageGap_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(gapSum/time_total);
}

/* Calculate the probability over time of being in the me2/me3 
   K27. */

double prob_me2_me3_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
	sumM++;
    }
    if (4*sumM > 3*c->sites)
      time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
  }

  return(time_in_M/r->t_out->el[r->t_outLastSample]);
}

/* Calculate the probability over time of being in the me2/me3 
   K27 (for the last hour of each cell cycle only). Do not 
   try to include results beyond p.cellCycles. */

double prob_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double time_in_M = 0.0;
  double time_total = 0.0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
          sumM++;
      }
      
      if (4*sumM > 3*c->sites) {
        time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  
  if (time_total == 0.0)
    fprintf(stderr,"Error: prob_me2_me3_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(time_in_M/time_total);
}

/* Calculate the probability over time of being in the me0/me1 
   K27. */

double prob_me0_me1_nCycles(chromatin *c, parameters *p, record *r) {
  long sumU, t, pos;
  double time_in_U = 0;
  
  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    sumU = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me0 || r->K27->el[pos][t]==me1)
	sumU++;
    }
    if (4*sumU > 3*c->sites)
      time_in_U += r->t_out->el[t] - r->t_out->el[t-1];
  }

  return(time_in_U/r->t_out->el[r->t_outLastSample]);
}

/* Calculate the probability over time of being in the me2/me3 
   K27 (for the last hour of each cell cycle only). Do not 
   try to include results beyond p.cellCycles. */

double prob_me0_me1_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumU = 0, t, pos;
  double time_in_U = 0.0;
  double time_total = 0.0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumU = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me0 || r->K27->el[pos][t]==me1)
          sumU++;
      }
      if (4*sumU > 3*c->sites) {
        time_in_U += r->t_out->el[t] - r->t_out->el[t-1];
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  
  if (time_total == 0.0)
    fprintf(stderr,"Error: prob_me0_me1_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(time_in_U/time_total);
}

/* Calculate the average number of histones in me2/me3 over time do
   not try to include results beyond p.cellCycles. */

double tAverage_me2_me3_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double Mavg = 0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
	sumM++;
    }
    Mavg += (double)sumM*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
  }
  return(Mavg/r->t_out->el[r->t_outLastSample-1]);
}

/* See above - but for last hour of each cell cycle only */

double tAverage_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double Mavg = 0.0, time_total = 0.0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
          sumM++;
      }
      Mavg += (double)sumM*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
    }
    time_total += r->t_out->el[t] - r->t_out->el[t-1];
  }

  if (time_total == 0.0)
    fprintf(stderr,"Error: tAverage_me2_me3_lastHour_nCycles. No samples for last hour of cell cycle.\n");
  
  return(Mavg/time_total);
}

/* Calculate the number of times the K27 switches from high me0/me1
   to high me2/me3. */

unsigned long numberHistoneStateFlips(record *r) {
  signed char newState = 0, oldState = 0;
  unsigned long t, pos, flips=0, m=0, u=0;;
  
  for (t=0;t<r->K27->cols && t<r->t_outLastSample;t++) {
    
    oldState = newState;
    
    m = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3) m++;
    }
    u = r->K27->rows - m;

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

/* Calculate the first passage time (taking care not to exceed
   p.t_outLastSample. */

double firstPassageTime(record *r, signed char *initial) {
  long unsigned m=0, u=0, pos,t=0;
  
  /* find initial state */
  for (pos=0;pos<r->K27->rows;pos++) {
    if (r->K27->el[pos][0]==me2 || r->K27->el[pos][0]==me3) m++;
  }
  u = r->K27->rows - m;
  
  if (m > u)
    *initial = -1;
  else
    *initial = 1;
    
  while ( t < r->K27->cols && t < r->t_outLastSample &&
	  ((*initial == -1 && 3*m > u) || (*initial == 1 && 3*u > m))) {
    m = 0;
    u = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3) m++;
    }
    u = r->K27->rows - m;
    t++;
  }
  return(r->t_out->el[t-1]);
}

/* Write the log file. */

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
  fprintf(fptr,"controlSites: %ld\n", c->controlSites);
  fprintf(fptr,"loci: %ld\n", p->loci);
  fprintf(fptr,"maxReact: %ld\n", p->maxReact);
  fprintf(fptr,"samples: %ld\n", p->samples);
  fprintf(fptr,"DNAreplication:");
  if (p->DNAreplication == TRUE)
    fprintf(fptr," TRUE\n");
  else
    fprintf(fptr," FALSE\n");
  fprintf(fptr,"resultsLastHour:");
  if (p->resultsLastHourOnly == TRUE)
    fprintf(fptr," TRUE\n");
  else
    fprintf(fptr," FALSE\n");
  fprintf(fptr,"SILAC:");
  if (p->SILAC == TRUE)
    fprintf(fptr," TRUE\n");
  else
    fprintf(fptr," FALSE\n");
  fprintf(fptr,"cellCycleDuration: %0.2f hours\n", p->cellCycleDuration);
  fprintf(fptr,"G2duration: %0.2f hours\n", p->G2duration);
  fprintf(fptr,"cellCycles: %d\n", p->cellCycles);
  fprintf(fptr,"me0_me1: %0.10f\n", p->me0_me1);
  fprintf(fptr,"me1_me2: %0.10f\n", p->me1_me2);
  fprintf(fptr,"me2_me3: %0.10f\n", p->me2_me3);
  fprintf(fptr,"me2factor: %0.4f\n", p->me2factor);
  fprintf(fptr,"me3factor: %0.4f\n", p->me3factor);
  fprintf(fptr,"firingRateMax: %0.10f\n", p->firingRateMax);
  fprintf(fptr,"firingRateMin: %0.10f\n", p->firingRateMin);
  fprintf(fptr,"transcription_demethylate: %0.6f\n", p->transcription_demethylate);

  return(1);
}
