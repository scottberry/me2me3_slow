#include "definitions.h"
/* 
   Output functions for calulating results and writing files.
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

void allocateSilacRecordMemory(chromatin *c, parameters *p, record *r) {

  r->silac = i_mat_get(c->sites,p->samples);
  p->SILAC_0h = (double)3600*p->cellCycleDuration*(p->silacLightCycles+1);
  p->SILAC_10h = p->SILAC_0h + (double)3600*10;
  p->SILAC_24h = p->SILAC_0h + (double)3600*24;
  p->SILAC_48h = p->SILAC_0h + (double)3600*48;
  
  r->silacResultsLight_0h = d_vec_get(p->loci);
  r->silacResultsLight_10h = d_vec_get(p->loci);
  r->silacResultsLight_24h = d_vec_get(p->loci);
  r->silacResultsLight_48h = d_vec_get(p->loci);
  r->silacResultsHeavy_0h = d_vec_get(p->loci);
  r->silacResultsHeavy_10h = d_vec_get(p->loci);
  r->silacResultsHeavy_24h = d_vec_get(p->loci);
  r->silacResultsHeavy_48h = d_vec_get(p->loci);

  return;
}

void freeSilacRecordMemory(record *r) {

  i_mat_free(r->silac);
  d_vec_free(r->silacResultsLight_0h);
  d_vec_free(r->silacResultsLight_10h);
  d_vec_free(r->silacResultsLight_24h);
  d_vec_free(r->silacResultsLight_48h);
  d_vec_free(r->silacResultsHeavy_0h);
  d_vec_free(r->silacResultsHeavy_10h);
  d_vec_free(r->silacResultsHeavy_24h);
  d_vec_free(r->silacResultsHeavy_48h);

  return;
}

void incrementSilacReportPoint(parameters *p) {
  if (p->SILAC_report == 1) {
    p->SILAC_report = 2;
    p->SILAC_nextReport = p->SILAC_10h;
  } else if (p->SILAC_report == 2) {
    p->SILAC_report = 3;
    p->SILAC_nextReport = p->SILAC_24h;
  } else if (p->SILAC_report == 3) {
    p->SILAC_nextReport = p->SILAC_48h;
    p->SILAC_report = 4;
  } else if (p->SILAC_report == 4) {
    p->SILAC_report = 5;
  }
  return;
}

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
  q->totalHistoneTurnover = 0.0;
  q->alphaSD = 0.0;
  q->alphaMean = 0.0;
  return;
}

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
  if (p->stochasticAlpha == TRUE) {
    q->alphaMean += tAverageAlpha(r)/p->loci;
    q->alphaSD += tAverageAlphaSD(r,q->alphaMean)/p->loci;
  }
  // Note: this metric works best when threshold = 1.0
  // q->firstPassage = firstPassageTime(r,&q->initial);

  // Note: this metric requires firing rate changes to work
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
  sprintf(tmp,"%0.2f",p->alpha);
  sprintf(ptmp,"a%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  sprintf(tmp,"%0.2f",p->beta);
  sprintf(ptmp,"b%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  sprintf(tmp,"%0.2f",p->firingThreshold);
  sprintf(ptmp,"thresh%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  sprintf(tmp,"%0.8f",p->transcription_turnover);
  sprintf(ptmp,"turn%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
  sprintf(tmp,"%0.2f",p->PRC2inhibition);
  sprintf(ptmp,"p%s",str_replace(tmp,decimal,underscore)); strcat(avgfile,ptmp);
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

void fprintParameterSpaceHeader(FILE *parFile) {
  fprintf(parFile,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\
\tFIRING_THRESHOLD\tP_DEMETHYLATE\tP_METHYLATE\tP_TURNOVER\
\tcontrolSites\talpha\talphaMean\
\talphaSD\tbeta\ttau\tgap\tMavg\
\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU\
\tavgInitU\ttTot\tprobM\tprobU\tbistability\tme3_end\ttotTurnover\n");
  fprintf(stderr,"me0_me1\tme1_me2\tme2_me3\tme2factor\tme3factor\tFIRING\
\tFIRING_THRESHOLD\tP_DEMETHYLATE\tP_METHYLATE\tP_TURNOVER\
\tcontrolSites\talpha\talphaMean\
\talphaSD\tbeta\ttau\tgap\tMavg\
\tlifetime\tinitM\tfirstPassageM\tavgInitM\tinitU\tfirstPassageU\
\tavgInitU\ttTot\tprobM\tprobU\tbistability\tme3_end\ttotTurnover\n");
  return;
}

void fprintParameterSpaceResults(FILE *parFile, parameters *p, chromatin *c, quantification *q) {        
  fprintf(parFile,"%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\
\t%0.10f\t%0.10f\t%ld\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f\t%ld\t%0.4f\t%0.4f \
\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.6f\t%0.10f\n",
          p->me0_me1,p->me1_me2,p->me2_me3,p->me2factor,p->me3factor,
          p->firingRateMax,p->firingThreshold,
          p->transcription_demethylate,p->me2_me3,p->transcription_turnover,c->controlSites,p->alpha,q->alphaMean,q->alphaSD,
          p->beta,p->G2duration,
          q->gap/p->loci,q->Mavg/p->loci,q->lifetime,q->initM,q->fpM,q->tM,q->initU,q->fpU,q->tU,q->tTot/p->loci,
          q->probM/p->loci,q->probU/p->loci,q->bistability,q->me3_end/p->loci,q->totalHistoneTurnover/p->loci);
  fprintf(stderr,"%0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  %0.10f  \
%0.10f %0.10f  %ld  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f  %ld  %0.4f  %0.4f  %ld  %0.4f  %0.4f \
%0.4f  %0.4f  %0.4f  %0.4f  %0.6f  %0.10f\n",
          p->me0_me1,p->me1_me2,p->me2_me3,p->me2factor,p->me3factor,
          p->firingRateMax,p->firingThreshold,
          p->transcription_demethylate,p->me2_me3,p->transcription_turnover,c->controlSites,p->alpha,q->alphaMean,q->alphaSD,
          p->beta,p->G2duration,
          q->gap/p->loci,q->Mavg/p->loci,q->lifetime,q->initM,q->fpM,q->tM,q->initU,q->fpU,q->tU,q->tTot/p->loci,
          q->probM/p->loci,q->probU/p->loci,q->bistability,q->me3_end/p->loci,q->totalHistoneTurnover/p->loci);
  return;
}

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

/* Calculate the probability over time of being in the low expression
   state. */

double prob_lowExpression_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0, f_me2_me3, lower_quartile;

  lower_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)/4.0;
  
  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
        sumM++;
    }
    f_me2_me3 = (double)sumM/(double)c->sites;
    if (firingRate(p,f_me2_me3) <= lower_quartile)
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

/* Calculate the probability over time of being in low expression
   state. (for the last hour of each cell cycle only) */

double prob_lowExpression_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double time_in_M = 0.0;
  double time_total = 0.0;
  double f_me2_me3, lower_quartile;

  lower_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)/4.0;
    
  // fprintf(stderr,"lower_quartile = %0.10f \n",lower_quartile);
  
  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
          sumM++;
      }
      f_me2_me3 = (double)sumM/(double)c->sites;
      // fprintf(stderr,"f_me2_me3 = %0.4f, firing %0.4f\n",f_me2_me3,firingRate(p,f_me2_me3));
      if (firingRate(p,f_me2_me3) <= lower_quartile) {
        time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
    }
  }
  
  if (time_total == 0.0)
    fprintf(stderr,"Error: prob_lowExpression_lastHour_nCycles. No samples for last hour of cell cycle.\n");

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

/* Calculate the probability over time of being in the high expression
   state. */

double prob_highExpression_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM, t, pos;
  double time_in_M = 0, f_me2_me3, upper_quartile;

  upper_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)*3.0/4.0;
  
  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    sumM = 0;
    for (pos=0;pos<r->K27->rows;pos++) {
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
        sumM++;
    }
    f_me2_me3 = (double)sumM/(double)c->sites;
    if (firingRate(p,f_me2_me3) > upper_quartile)
      time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
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

/* Calculate the probability over time of being in high expression
   state. (for the last hour of each cell cycle only) */

double prob_highExpression_lastHour_nCycles(chromatin *c, parameters *p, record *r) {
  long sumM = 0, t, pos;
  double time_in_M = 0.0;
  double time_total = 0.0;
  double f_me2_me3, upper_quartile;

  upper_quartile = p->firingRateMin + (p->firingRateMax - p->firingRateMin)*3.0/4.0;
  
  for (t=1;t<r->K27->cols && t<r->t_outLastSample;t++) {
    if (fmod(r->t_out->el[t],3600*p->cellCycleDuration) >= 3600*(p->cellCycleDuration - 1)) { // if within last hour
      sumM = 0;
      for (pos=0;pos<r->K27->rows;pos++) {
        if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3)
          sumM++;
      }
      f_me2_me3 = (double)sumM/(double)c->sites;
      if (firingRate(p,f_me2_me3) >= upper_quartile) {
        time_in_M += r->t_out->el[t] - r->t_out->el[t-1];
      }
      time_total += r->t_out->el[t] - r->t_out->el[t-1];
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
    if (r->K27->el[pos][0]==me2 || r->K27->el[pos][0]==me3) m++;
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
      if (r->K27->el[pos][t]==me2 || r->K27->el[pos][t]==me3) m++;
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

/* Write the log file. */

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
  fprintf(fptr,"SILAC:");
  if (p->silacExperiment == TRUE) {
    fprintf(fptr," TRUE\n");
    fprintf(fptr,"SilacLightCycles = %ld \n",p->silacLightCycles);
    fprintf(fptr,"SilacHeavyCycles = %ld \n",p->silacHeavyCycles);
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
  fprintf(fptr,"firingThreshold: %0.10f\n", p->firingThreshold);
  fprintf(fptr,"stochasticAlpha:");
  if (p->stochasticAlpha == TRUE) {
    fprintf(fptr," TRUE\n");
    fprintf(fptr,"stochasticTranslationEfficiency: %ld\n",p->stochasticTranslationEfficiency);    
  } else {
    fprintf(fptr," FALSE\n");
    fprintf(fptr,"alpha: %0.4f\n", p->alpha);
  }
  fprintf(fptr,"beta: %0.4f\n", p->beta);
  fprintf(fptr,"transcription_demethylate: %0.10f\n", p->transcription_demethylate);
  fprintf(fptr,"transcription_turnover: %0.10f\n", p->transcription_turnover);
  fprintf(fptr,"noisy_demethylate: %0.10f\n", p->noisy_demethylate);
  
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

void fprintHistoneTurnover(FILE *fptr, parameters *p, record *r) {
  double turn_me0 = 0.0, turn_me1 = 0.0, turn_me2 = 0.0, turn_me3 = 0.0;
  long locus;
  
  fprintf(fptr,"species\tturnover\n");
  
  for (locus=0;locus<p->loci;locus++) {
    turn_me0 += r->turnover->el[locus][0];
    turn_me1 += r->turnover->el[locus][1];
    turn_me2 += r->turnover->el[locus][2];
    turn_me3 += r->turnover->el[locus][3];
  }

  fprintf(fptr,"me0\t%0.12f\n",turn_me0/p->loci);
  fprintf(fptr,"me1\t%0.12f\n",turn_me1/p->loci);
  fprintf(fptr,"me2\t%0.12f\n",turn_me2/p->loci);
  fprintf(fptr,"me3\t%0.12f\n",turn_me3/p->loci);
  fprintf(fptr,"total\t%0.12f\n",(turn_me0 + turn_me1 + turn_me2 + turn_me3)/p->loci);
  return; 
}

void fprintResultsFinalLocus(char *avgfile, record *r) {
  char fname[256]="";
  
  /* print results for final locus */
  strcpy(fname,"t_\0"); strcat(fname,avgfile);
  fprint_t_out_nCycles(fname,r);
  strcpy(fname,"me0_t_\0"); strcat(fname,avgfile);
  fprint_t_nCycles(fname,r->K27,me0,r);
  strcpy(fname,"me1_t_\0"); strcat(fname,avgfile);
  fprint_t_nCycles(fname,r->K27,me1,r);
  strcpy(fname,"me2_t_\0"); strcat(fname,avgfile);
  fprint_t_nCycles(fname,r->K27,me2,r);
  strcpy(fname,"me3_t_\0"); strcat(fname,avgfile);
  fprint_t_nCycles(fname,r->K27,me3,r);
  strcpy(fname,"Firing_t_\0"); strcat(fname,avgfile);
  fprint_firing_t_nCycles(fname,r);

  return;
}

void fprintSilacResultsFinalLocus(char *avgfile, record *r) {
  char fname[256]="";

  // light histones
  strcpy(fname,"LIGHT_me0_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me0,r->silac,LIGHT,r);
  strcpy(fname,"LIGHT_me1_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me1,r->silac,LIGHT,r);
  strcpy(fname,"LIGHT_me2_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me2,r->silac,LIGHT,r);
  strcpy(fname,"LIGHT_me3_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me3,r->silac,LIGHT,r);
    
  // heavy histones
  strcpy(fname,"HEAVY_me0_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me0,r->silac,HEAVY,r);
  strcpy(fname,"HEAVY_me1_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me1,r->silac,HEAVY,r);
  strcpy(fname,"HEAVY_me2_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me2,r->silac,HEAVY,r);
  strcpy(fname,"HEAVY_me3_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me3,r->silac,HEAVY,r);

  // unlabelled histones
  strcpy(fname,"UNLABELLED_me0_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me0,r->silac,UNLABELLED,r);
  strcpy(fname,"UNLABELLED_me1_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me1,r->silac,UNLABELLED,r);
  strcpy(fname,"UNLABELLED_me2_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me2,r->silac,UNLABELLED,r);
  strcpy(fname,"UNLABELLED_me3_t_\0"); strcat(fname,avgfile);
  fprint_silac_t_nCycles(fname,r->K27,me3,r->silac,UNLABELLED,r);    

  return;
}

void fprint_transFactorProtein_nCycles(char *fname, record *r) {
  FILE *fptr;
  long unsigned i;
  //fprintf(stderr,"%ld",r->t_outLastSample);
  fptr = fopen(fname,"w");
  fprintf(fptr,"time\tprotein\tRNA\talpha\n");
  for (i=0;i<r->t_outLastSample;i++) {
    fprintf(fptr,"%0.4f\t",r->t_out->el[i]);
    fprintf(fptr,"%ld\t",r->transFactorProtein->el[i]);
    fprintf(fptr,"%ld\t",r->transFactorRNA->el[i]);
    fprintf(fptr,"%0.4f\n",r->alpha->el[i]);
  }
  fclose(fptr);
  return;
}

double tAverageAlpha(record *r) {
  double mean = 0.0, tLast = 0.0;
  long t;

  for (t=0;t<r->K27->cols && t<r->t_outLastSample;t++) { 
    mean += r->alpha->el[t]*(double)(r->t_out->el[t]-r->t_out->el[t-1]);
    tLast = r->t_out->el[t];
  }
  
  return(mean/tLast);
}

double tAverageAlphaSD(record *r, double mean) {
  double var = 0.0, samplePoint = 0.0, sampleFreq;
  long reaction = 0, nSamples;

  // sample the distribution points evenly spaced in time
  nSamples = 10000;
  sampleFreq = r->tMax/(double)nSamples;

  while (samplePoint < r->tMax) {
    var += pow(mean - r->alpha->el[reaction],2.0);
    samplePoint += sampleFreq;
    while (r->t_out->el[reaction] <= samplePoint && reaction < r->t_out->len)
      reaction++;
  }

  return(sqrt(var/(double)nSamples));
}

