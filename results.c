#include "definitions.h"
/* 
   Output functions for calulating results and writing to file.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
 */

/* Replace a character of a string */
char *str_replace(char *orig, char *rep, char *with) {
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep
    int len_with; // length of with
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    if (!orig)
        return NULL;
    if (!rep)
        rep = "";
    len_rep = strlen(rep);
    if (!with)
        with = "";
    len_with = strlen(with);

    ins = orig;
    for (count = 0; (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }

    // first time through the loop, all the variable are set correctly
    // from here on,
    //    tmp points to the end of the result string
    //    ins points to the next occurrence of rep in orig
    //    orig points to the remainder of orig after "end of rep"
    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}

/* Print simulation time at sample points. Length depends directly on
   r.tMax, which is determined by p.cellCycles. */

void fprint_t_out_nCycles(char *fname, record *r) {
  FILE *fptr;
  long unsigned i;
  
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

/* Average results over time and print a time-dependent results 
   vector. Length depends directly on r.tMax, which is 
   determined by p.cellCycles. */

void fprint_silac_t_nCycles(char *fname, I_MAT *mat, int target, I_MAT *silac, int silac_target, record *r) {
  FILE *fptr;
  long unsigned count, i, j;
  
  fptr = fopen(fname,"w");
  for (i=0;i<mat->cols && i < r->t_outLastSample;i++) {
    count = 0;
    for (j=0;j<mat->rows;j++) {
      if (mat->el[j][i] == target && silac->el[j][i] == silac_target) {
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

/* Print transcription status for each time point
   Length depends directly on r.tMax, which is 
   determined by p.cellCycles */

void fprint_transcribing_t_nCycles(char *fname, record *r) {
  FILE *fptr;
  long unsigned i;

  fptr = fopen(fname,"w");
  fprintf(fptr,"t\ttranscribing\n");
  for (i=0;i<r->firing->len && r->t->el[i]<r->tMax;i++) {
    fprintf(fptr,"%0.4f\t%ld\n",r->t->el[i],r->transcribing->el[i]);
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
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }

  if (time_total == 0.0)
    fprintf(stderr,"Error: tAverage_me2_me3_lastHour_nCycles. No samples for last hour of cell cycle.\n");
  
  return(Mavg/time_total);
}

double tAverage_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double Mavg = 0.0, time_total = 0.0;

  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me3)
          sumM++;
      }
      Mavg += (double)sumM*(r->t_out->el[t]-r->t_out->el[t-1])/(c->sites);
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
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
  double fpt = 0.0;
  
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
  if (t==r->t_outLastSample)
    fpt = r->tMax;
  else 
    fpt = r->t_out->el[t-1];
  return(fpt);
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
  if (p->silacExperiment == TRUE) {
    fprintf(fptr," TRUE\n");
    fprintf(fptr," silacLightCycles = %ld \n",p->silacLightCycles);
  } else
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
  fprintf(fptr,"activation: %0.4f\n", p->activation);
  fprintf(fptr,"transcription_demethylate: %0.6f\n", p->transcription_demethylate);

  return(1);
}

void fprintTripleSILAC_eachLocus(FILE *fptrAbs, FILE *fptrRel, long locus, parameters *p, record *r) {
  int count_me0_LIGHT = 0;
  int count_me0_HEAVY = 0;
  int count_me1_LIGHT = 0;
  int count_me1_HEAVY = 0;
  int count_me2_LIGHT = 0;
  int count_me2_HEAVY = 0;
  int count_me3_LIGHT = 0;
  int count_me3_HEAVY = 0;
  int j, t, tot_HEAVY, tot_LIGHT;
  double time;
  
  t = p->reactCount-1;
  time = p->SILAC_nextReport;
  
  for (j=0;j<r->K27->rows;j++) {
    if (r->K27->el[j][t] == me0 && r->silac->el[j][t] == LIGHT) count_me0_LIGHT++;
    if (r->K27->el[j][t] == me0 && r->silac->el[j][t] == HEAVY) count_me0_HEAVY++;
    if (r->K27->el[j][t] == me1 && r->silac->el[j][t] == LIGHT) count_me1_LIGHT++;
    if (r->K27->el[j][t] == me1 && r->silac->el[j][t] == HEAVY) count_me1_HEAVY++;
    if (r->K27->el[j][t] == me2 && r->silac->el[j][t] == LIGHT) count_me2_LIGHT++;
    if (r->K27->el[j][t] == me2 && r->silac->el[j][t] == HEAVY) count_me2_HEAVY++;
    if (r->K27->el[j][t] == me3 && r->silac->el[j][t] == LIGHT) count_me3_LIGHT++;
    if (r->K27->el[j][t] == me3 && r->silac->el[j][t] == HEAVY) count_me3_HEAVY++;
  }

  /* Print absolute values */
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme0\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me0_LIGHT/(double)r->K27->rows);
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme0\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me0_HEAVY/(double)r->K27->rows);
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme1\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me1_LIGHT/(double)r->K27->rows);
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme1\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me1_HEAVY/(double)r->K27->rows);
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme2\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me2_LIGHT/(double)r->K27->rows);
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme2\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me2_HEAVY/(double)r->K27->rows);
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme3\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me3_LIGHT/(double)r->K27->rows);
  fprintf(fptrAbs,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme3\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me3_HEAVY/(double)r->K27->rows);

  /* Print relative values */
  tot_LIGHT = count_me0_LIGHT + count_me1_LIGHT + count_me2_LIGHT + count_me3_LIGHT;
  tot_HEAVY = count_me0_HEAVY + count_me1_HEAVY + count_me2_HEAVY + count_me3_HEAVY;
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\ttot\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)tot_LIGHT/r->K27->rows);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\ttot\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)tot_HEAVY/r->K27->rows);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme0\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me0_LIGHT/(double)tot_LIGHT);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme0\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me0_HEAVY/(double)tot_HEAVY);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme1\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me1_LIGHT/(double)tot_LIGHT);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme1\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me1_HEAVY/(double)tot_HEAVY);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme2\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me2_LIGHT/(double)tot_LIGHT);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme2\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me2_HEAVY/(double)tot_HEAVY);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme3\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me3_LIGHT/(double)tot_LIGHT);
  fprintf(fptrRel,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\tme3\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,locus,time,(double)count_me3_HEAVY/(double)tot_HEAVY);

  return;
} 

void storeTripleSILAC_me3(long locus, parameters *p, record *r) {
  int count_me3_LIGHT = 0;
  int count_me3_HEAVY = 0;
  int j, t, tot_HEAVY = 0, tot_LIGHT = 0;
  double time;
  
  t = p->reactCount-1;
  time = p->SILAC_nextReport;
  
  for (j=0;j<r->K27->rows;j++) {
    if (r->K27->el[j][t] == me3 && r->silac->el[j][t] == LIGHT) count_me3_LIGHT++;
    if (r->K27->el[j][t] == me3 && r->silac->el[j][t] == HEAVY) count_me3_HEAVY++;
    if (r->silac->el[j][t] == LIGHT) tot_LIGHT++;
    if (r->silac->el[j][t] == HEAVY) tot_HEAVY++;
  }

  /* Store relative values */
  if (p->SILAC_report == 1) {
    if (tot_LIGHT == 0) r->silacResultsLight_0h->el[locus] = 0.0;
    else r->silacResultsLight_0h->el[locus] = (double)count_me3_LIGHT/(double)tot_LIGHT;
    if (tot_HEAVY == 0) r->silacResultsHeavy_0h->el[locus] = 0.0;
    else r->silacResultsHeavy_0h->el[locus] = (double)count_me3_HEAVY/(double)tot_HEAVY;
  }
  if (p->SILAC_report == 2) {
    if (tot_LIGHT == 0) r->silacResultsLight_10h->el[locus] = 0.0;
    else r->silacResultsLight_10h->el[locus] = (double)count_me3_LIGHT/(double)tot_LIGHT;
    if (tot_HEAVY == 0) r->silacResultsLight_10h->el[locus] = 0.0;
    else r->silacResultsHeavy_10h->el[locus] = (double)count_me3_HEAVY/(double)tot_HEAVY;
  }
  if (p->SILAC_report == 3) {
    if (tot_LIGHT == 0) r->silacResultsLight_24h->el[locus] = 0.0;
    else r->silacResultsLight_24h->el[locus] = (double)count_me3_LIGHT/(double)tot_LIGHT;
    if (tot_HEAVY == 0) r->silacResultsLight_24h->el[locus] = 0.0;
    else r->silacResultsHeavy_24h->el[locus] = (double)count_me3_HEAVY/(double)tot_HEAVY;
  }
  if (p->SILAC_report == 4) {
    if (tot_LIGHT == 0) r->silacResultsLight_48h->el[locus] = 0.0;
    else r->silacResultsLight_48h->el[locus] = (double)count_me3_LIGHT/(double)tot_LIGHT;
    if (tot_HEAVY == 0) r->silacResultsLight_48h->el[locus] = 0.0;
    else r->silacResultsHeavy_48h->el[locus] = (double)count_me3_HEAVY/(double)tot_HEAVY;
  }
  return;
} 


void fprintTripleSILAC_average(FILE *fptr, parameters *p, record *r) {
  double L0, L10, L24, L48;
  double H0, H10, H24, H48;
  int locus;
  double time;

  L0 = L10 = L24 = L48 = 0.0;
  H0 = H10 = H24 = H48 = 0.0;
  
  for (locus=0;locus<p->loci;locus++) {
    L0 += r->silacResultsLight_0h->el[locus];
    H0 += r->silacResultsHeavy_0h->el[locus];
    L10 += r->silacResultsLight_10h->el[locus];
    H10 += r->silacResultsHeavy_10h->el[locus];
    L24 += r->silacResultsLight_24h->el[locus];
    H24 += r->silacResultsHeavy_24h->el[locus];
    L48 += r->silacResultsLight_48h->el[locus];
    H48 += r->silacResultsHeavy_48h->el[locus];
    //fprintf(stderr,"L48 = %0.10f\n",L48);
  }

  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t0.0\tme3\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,L0/(double)p->loci);
  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t0.0\tme3\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,H0/(double)p->loci);
  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t10.0\tme3\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,L10/(double)p->loci);
  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t10.0\tme3\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,H10/(double)p->loci);
  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t24.0\tme3\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,L24/(double)p->loci);
  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t24.0\tme3\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,H24/(double)p->loci);
  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t48.0\tme3\tLIGHT\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,L48/(double)p->loci);
  fprintf(fptr,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t48.0\tme3\tHEAVY\t%0.10f\n",p->firingRateMax,p->transcription_demethylate,p->me2_me3,p->firingThreshold,H48/(double)p->loci);
  
  return;
}
