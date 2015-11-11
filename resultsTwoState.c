#include "definitions.h"
/* 
   Memory and output functions for calulating results and writing files.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

/* Reset all quantifications for different parameter sets */
void resetQuantification(quantification *q) {
  q->gap = 0.0;
  q->Mavg = 0.0;
  q->me3_end = 0.0;
  q->probM = 0.0;
  q->probU = 0.0;
  q->fh = 0;
  q->tTot = q->tTotM = q->tTotU = 0.0;
  q->initM = q->initU = 0;
  q->firstPassageM = q->firstPassageU = 0.0;
  q->firingEvents = 0;
  return;
}

/* Accumulate quantified results for each locus, so specific
   trajectories can be discarded and over-written in memory */
void accumulateQuantification(chromatin *c, parameters *p, record *r, quantification *q) {  

  if (p->resultsLastHourOnly == TRUE) {
    q->gap += tAverageGap_lastHour_nCycles(c,p,r);
    q->Mavg += tAverage_me2_me3_lastHour_nCycles(c,p,r);
    q->probM += prob_lowExpression_lastHour_nCycles(c,p,r);
    q->probU += prob_highExpression_lastHour_nCycles(c,p,r);
    // q->probM += prob_me2_me3_lastHour_nCycles(c,p,r);
    // q->probU += prob_me0_me1_lastHour_nCycles(c,p,r);
  } else {
    q->gap += tAverageGap_nCycles(c,p,r);
    q->Mavg += tAverage_me2_me3_nCycles(c,p,r);
    q->probM += prob_me2_me3_nCycles(c,p,r);
    q->probU += prob_me0_me1_nCycles(c,p,r);
  }
  q->me3_end += tAverage_me3_lastHour_nCycles(c,p,r);
  q->fh += numberHistoneStateFlips(r);
  q->tTot += r->t->el[p->reactCount];
  
  // Note: this metric works best when threshold = 1.0
  // (based on histone state)
  // q->firstPassage = firstPassageTime(r,&q->initial);

  // Note: this metric requires firing rate changes
  // (based on expression quartiles)
  q->firstPassage = firstPassageTimeExpression(r,p,&q->initial);

  if (q->initial==-1) {
    q->firstPassageM += q->firstPassage;
    q->initM++;
    q->tTotM += r->t->el[p->reactCount];
  } else {
    q->firstPassageU += q->firstPassage;
    q->initU++;
    q->tTotU += r->t->el[p->reactCount];
  }
  
  return;
}

void averageQuantification(chromatin *c, parameters *p, record *r, quantification *q) {  

  if (q->fh != 0) q->lifetime = q->tTot/q->fh;
  else q->lifetime = -1.0;

  q->bistability = 4*q->probM*q->probU/(p->loci*p->loci);

  if (q->initM != 0) {
    q->fpM = q->firstPassageM/q->initM;
    q->tM = q->tTotM/q->initM;
  } else {
    q->fpM = -1.0;
    q->tM = -1.0;
  }

  if (q->initU != 0) {
    q->fpU= q->firstPassageU/q->initU;
    q->tU = q->tTotU/q->initU;
  } else {
    q->fpU = -1.0;
    q->tU = -1.0;
  }

  return;
}

/* Create a parameter-dependent filename */
char *parameterDependentBasename(chromatin *c, parameters *p) {
  char *avgfile;
  char *decimal = ".", *underscore = "_";
  char ptmp[256]="", tmp[256]="";

  avgfile = malloc(256);
  sprintf(tmp,"s%ld",c->sites); strcat(avgfile,tmp); 
  sprintf(tmp,"ctrl%ld",c->controlSites); strcat(avgfile,tmp);
  sprintf(tmp,"cc%d",p->cellCycles); strcat(avgfile,tmp);
  sprintf(tmp,"%0.2f",p->firingThreshold);
  sprintf(ptmp,"thresh%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  if (p->DNAreplication == TRUE) {
    sprintf(tmp,"Rep"); strcat(avgfile,tmp);
  } else {
    sprintf(tmp,"NoRep"); strcat(avgfile,tmp); 
  }
  sprintf(tmp,"_st%ld",p->optimSteps); strcat(avgfile,tmp); 

  strcat(avgfile,p->id);
  strcat(avgfile,".txt\0");

  return(avgfile);
}

/* Print header for the results files */
void fprintParameterSpaceHeader(FILE *parFile) {
  fprintf(parFile,"FIRING\
\tFIRING_THRESHOLD\tP_DEMETHYLATE\tP_METHYLATE\
\tcontrolSites\tgap\tMavg\
\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU\
\tavgInitU\ttTot\tprobM\tprobU\tbistability\tme3_end\n");
  fprintf(stderr,"FIRING\
\tFIRING_THRESHOLD\tP_DEMETHYLATE\tP_METHYLATE\
\tcontrolSites\tgap\tMavg\
\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU\
\tavgInitU\ttTot\tprobM\tprobU\tbistability\tme3_end\n");
  return;
}

/* Print accumulated results for each parameter set */
void fprintParameterSpaceResults(FILE *parFile, parameters *p, chromatin *c, quantification *q) {        
  fprintf(parFile,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%ld\t%0.4f\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f\
\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.6f\n",
          p->firingRateMax,
          p->firingThreshold,
          p->transcription_demethylate,
          p->me2_me3,
          c->controlSites,
          q->gap/p->loci,
          q->Mavg/p->loci,
          q->lifetime,
          q->initM,
          q->fpM,
          q->tM,
          q->initU,
          q->fpU,
          q->tU,
          q->tTot/p->loci,
          q->probM/p->loci,
          q->probU/p->loci,
          q->bistability,
          q->me3_end/p->loci);
  fprintf(stderr,"%0.10f  %0.10f  %0.10f  %0.10f  %ld  %0.4f  %0.4f  %0.4f  %ld  %0.4f  %0.4f  %ld  %0.4f  %0.4f \
%0.4f  %0.4f  %0.4f  %0.4f  %0.6f\n",
          p->firingRateMax,
          p->firingThreshold,
          p->transcription_demethylate,
          p->me2_me3,
          c->controlSites,
          q->gap/p->loci,
          q->Mavg/p->loci,
          q->lifetime,
          q->initM,
          q->fpM,
          q->tM,
          q->initU,
          q->fpU,
          q->tU,
          q->tTot/p->loci,
          q->probM/p->loci,
          q->probU/p->loci,
          q->bistability,
          q->me3_end/p->loci);
  return;
}

/* Replace a character of a string 
   Courtesy of jmucchiello 
   http://stackoverflow.com/questions/779875/what-is-the-function-to-replace-string-in-c */
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

/* Average results over locus and print a time-dependent results 
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
   function evaluation accounts for non-constant time-step. Average
   between samples 1 and r.t_outLastSample. Determined by whether 
   p.maxReact or steps required to simulate p.cellCycles is 
   greater. */

double tAverageGap_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double gapSum = 0;

  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me3)
        sumM++;
    }
    gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t+1]-r->t_out->el[t])/(c->sites);
  }
  return(gapSum/r->t_out->el[r->t_outLastSample]);
}

/* See above - but for last hour of each cell cycle only */

double tAverageGap_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double gapSum = 0.0, time_total = 0.0;

  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me3)
          sumM++;
      }
      gapSum += (double)labs(2*sumM-c->sites)*(r->t_out->el[t+1]-r->t_out->el[t])/(c->sites);
      time_total += r->t_out->el[t+1] - r->t_out->el[t];
    }
  }

  if (time_total == 0.0)
    fprintf(stderr,"Error: tAverageGap_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(gapSum/time_total);
}

/* Calculate the probability over time of being in the me2/me3 
   K27 state */
double prob_me2_me3_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0;

  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me3)
        sumM++;
    }
    if (4*sumM > 3*c->sites)
      time_in_M += r->t_out->el[t+1] - r->t_out->el[t];
  }

  return(time_in_M/r->t_out->el[r->t_outLastSample]);
}

/* Calculate the probability over time of being in the lower
   expression quartile */
double prob_lowExpression_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0, f_me2_me3, lower_quartile;

  lower_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)/4.0;
  
  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me3)
        sumM++;
    }
    f_me2_me3 = (double)sumM/(double)c->sites;
    if (firingRate(p,f_me2_me3) <= lower_quartile)
      time_in_M += r->t_out->el[t+1] - r->t_out->el[t];
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

  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me3)
          sumM++;
      }
      
      if (4*sumM > 3*c->sites) {
        time_in_M += r->t_out->el[t+1] - r->t_out->el[t];
      }
      time_total += r->t_out->el[t+1] - r->t_out->el[t];
    }
  }
  
  if (time_total == 0.0)
    fprintf(stderr,"Error: prob_me2_me3_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(time_in_M/time_total);
}

/* Calculate the probability over time of being in low expression
   state. (for the last hour of each cell cycle only) */
double prob_lowExpression_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double time_in_M = 0.0;
  double time_total = 0.0;
  double f_me2_me3, lower_quartile;

  lower_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)/4.0;
    
  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me3)
          sumM++;
      }
      f_me2_me3 = (double)sumM/(double)c->sites;
      if (firingRate(p,f_me2_me3) <= lower_quartile) {
        time_in_M += r->t_out->el[t+1] - r->t_out->el[t];
      }
      time_total += r->t_out->el[t+1] - r->t_out->el[t];
    }
  }
  
  if (time_total == 0.0)
    fprintf(stderr,"Error: prob_lowExpression_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(time_in_M/time_total);
}

/* Calculate the probability over time of being in the me0/me1 
   K27 state */
double prob_me0_me1_nCycles(chromatin *c, parameters *p, record *r) {
  long sumU, t, pos;
  double time_in_U = 0;
  
  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    sumU = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me0)
        sumU++;
    }
    if (4*sumU > 3*c->sites)
      time_in_U += r->t_out->el[t+1] - r->t_out->el[t];
  }

  return(time_in_U/r->t_out->el[r->t_outLastSample]);
}

/* Calculate the probability over time of being in the high expression
   quartile */
double prob_highExpression_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0, f_me2_me3, upper_quartile;

  upper_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)*3.0/4.0;
  
  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me3)
        sumM++;
    }
    f_me2_me3 = (double)sumM/(double)c->sites;
    if (firingRate(p,f_me2_me3) > upper_quartile)
      time_in_M += r->t_out->el[t+1] - r->t_out->el[t];
  }

  return(time_in_M/r->t_out->el[r->t_outLastSample]);
}

/* Calculate the probability over time of being in the me2/me3 
   K27 (for the last hour of each cell cycle only). Do not 
   try to include results beyond p.cellCycles. */
double prob_me0_me1_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumU = 0, t, pos;
  double time_in_U = 0.0;
  double time_total = 0.0;

  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumU = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me0)
          sumU++;
      }
      if (4*sumU > 3*c->sites) {
        time_in_U += r->t_out->el[t+1] - r->t_out->el[t];
      }
      time_total += r->t_out->el[t+1] - r->t_out->el[t];
    }
  }
  
  if (time_total == 0.0)
    fprintf(stderr,"Error: prob_me0_me1_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(time_in_U/time_total);
}

/* Calculate the probability over time of being in high expression
   state. (for the last hour of each cell cycle only) */
double prob_highExpression_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double time_in_M = 0.0;
  double time_total = 0.0;
  double f_me2_me3, upper_quartile;

  upper_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)*3.0/4.0;
  
  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me3)
          sumM++;
      }
      f_me2_me3 = (double)sumM/(double)c->sites;
      if (firingRate(p,f_me2_me3) >= upper_quartile) {
        time_in_M += r->t_out->el[t+1] - r->t_out->el[t];
      }
      time_total += r->t_out->el[t+1] - r->t_out->el[t];
    }
  }
  
  if (time_total == 0.0)
    fprintf(stderr,"Error: prob_highExpression_lastHour_nCycles. No samples for last hour of cell cycle.\n");

  return(time_in_M/time_total);
}

/* Calculate the average number of histones in me2/me3 over time do
   not try to include results beyond p.cellCycles. */
double tAverage_me2_me3_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double Mavg = 0;

  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me3)
        sumM++;
    }
    Mavg += (double)sumM*(r->t_out->el[t+1]-r->t_out->el[t])/(c->sites);
  }
  return(Mavg/r->t_out->el[r->t_outLastSample-1]);
}

/* See above - but for last hour of each cell cycle only */
double tAverage_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double Mavg = 0.0, time_total = 0.0;

  for (t=0;t<r->K27->cols-1 && t<r->t_outLastSample-1;t++) {
    if (fmod(r->t_out->el[t],p->cellCycleDuration*3600) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me3)
          sumM++;
      }
      Mavg += (double)sumM*(r->t_out->el[t+1]-r->t_out->el[t])/(c->sites);
      time_total += r->t_out->el[t+1] - r->t_out->el[t];
    }
  }

  if (time_total == 0.0)
    fprintf(stderr,"Error: tAverage_me2_me3_lastHour_nCycles. No samples for last hour of cell cycle.\n");
  
  return(Mavg/time_total);
}

/* Calculate average cell-cycle end value of me3 */
double tAverage_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, sample, lastSample, pos, cycle;
  double Mavg = 0.0, t_end;

  cycle = 1;
  t_end = (double)cycle * p->cellCycleDuration * 3600;
  while (t_end < r->tMax) { // check t_end in range
    // find last sample in each cell cycle
    sample = lastSample = 0;
    while (r->t_out->el[sample] < t_end) {
      lastSample = sample;
      sample++;
    }
    // calculate me3 level for the last sample
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][lastSample]==me3)
        sumM++;
    }
    Mavg += (double)sumM/c->sites;
    
    // increment cycle and t_end
    cycle++;
    t_end = (double)cycle * p->cellCycleDuration * 3600;
  }
  
  return(Mavg/(double)(cycle-1));
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
      if (r->K27->el[pos][t]==me3) m++;
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
    if (r->K27->el[pos][0]==me3) m++;
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
      if (r->K27->el[pos][t]==me3) m++;
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

/* Calculate the first passage time (taking care not to exceed
   p.t_outLastSample. */

double firstPassageTimeExpression(record *r, parameters *p, signed char *initial) {
  long unsigned m=0, pos,t=0;
  double lower_quartile, median, upper_quartile, f;
  double fpt = 0.0;

  lower_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)/4.0;
  median = p->firingRateMin + (p->firingRateMax - p->firingRateMin)/2.0;
  upper_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)*3.0/4.0;

  // fprintf(stderr,"upr = %0.8f,med = %0.8f,lwr = %0.8f,",upper_quartile,median,lower_quartile);
  
  /* find initial state */
  for (pos=0;pos<r->K27->rows;pos++) {
    if (r->K27->el[pos][0]==me3) m++;
  }

  f = firingRate(p,(double)m/r->K27->rows);
  //   fprintf(stderr,"m = %ld,f = %0.8f\n",m,f);
  
  /* find initial state */
  if (f < median)
    *initial = -1;
  else
    *initial = 1;
  
  while ( t < r->K27->cols && t < r->t_outLastSample &&
          ((*initial == -1 && f < upper_quartile) || (*initial == 1 && f > lower_quartile))) {
    m = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me3) m++;
    }
    f = firingRate(p,(double)m/r->K27->rows);
    t++;
  }
  if (t==r->t_outLastSample)
    fpt = r->tMax;
  else 
    fpt = r->t_out->el[t-1];
  return(fpt);
}

/* Print log file */
int writelog(FILE *fptr, chromatin *c, parameters *p, record *r) {
  time_t curtime;
  struct tm *loctime;

  curtime = time (NULL);
  loctime = localtime (&curtime);
  fputs (asctime (loctime), fptr);

  fprintf(fptr,"Program name: %s\n",p->executable);

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
  fprintf(fptr,"optimSteps: %ld\n", p->optimSteps);
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
  fprintf(fptr,"cellCycleDuration: %0.2f hours\n", p->cellCycleDuration);
  fprintf(fptr,"cellCycles: %d\n", p->cellCycles);
  fprintf(fptr,"me2_me3: %0.10f\n", p->me2_me3);
  fprintf(fptr,"firingRateMax: %0.10f\n", p->firingRateMax);
  fprintf(fptr,"firingRateMin: %0.10f\n", p->firingRateMin);
  fprintf(fptr,"firingThreshold: %0.10f\n", p->firingThreshold);
  fprintf(fptr,"transcription_demethylate: %0.10f\n", p->transcription_demethylate);
  
  return(1);
}

/* Interface for printing all time-dependent results files for the
   final locus */
void fprintResultsFinalLocus(char *avgfile, record *r) {
  char fname[256]="";
  
  /* print results for final locus */
  strcpy(fname,"t_\0"); strcat(fname,avgfile);
  fprint_t_out_nCycles(fname,r);
  strcpy(fname,"me0_t_\0"); strcat(fname,avgfile);
  fprint_t_nCycles(fname,r->K27,me0,r);
  strcpy(fname,"me3_t_\0"); strcat(fname,avgfile);
  fprint_t_nCycles(fname,r->K27,me3,r);
  strcpy(fname,"Firing_t_\0"); strcat(fname,avgfile);
  fprint_firing_t_nCycles(fname,r);

  return;
}

/* Calculate the number of firing events in the last cell cycle */
long countFiringEventsLastCellCycle(parameters *p, record *r) {
  long i, count = 0;
  double start, end;

  // find start time for last cell cycle
  start = r->tMax - 3600*p->cellCycleDuration;
  end = r->tMax;

  for(i=0;i<r->t->len;i++) {
    if(r->t->el[i] > start && r->t->el[i] < end && r->firing->el[i]==TRUE) {
      count++;
    }
  }
  return(count);
}

/* Reset the stored record of transcriptional firing times */
void resetFiringRecord(record *r) {
  long i;
  for (i=0;i<r->firing->len;i++)
    r->firing->el[i]=FALSE;
  return;
}
