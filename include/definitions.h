/* 
   Definitions file for gillespie algorithm simulations of models
   in chromatin-based epigenetics. 
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "scottsmatrices.h"
#include "time.h"
#include "math.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef __APPLE__
#include <mach/mach_time.h>
#endif

#define FALSE 0
#define TRUE 1

#define me0 0
#define me1 1
#define me2 2
#define me3 3

typedef unsigned char logical;

typedef struct {
  long sites;
  long controlSites;
  I_VEC *K27;
} chromatin;

typedef struct {
  logical histone, transcribed;
} flags;

typedef struct {
  // random number parameters
  const gsl_rng_type *gsl_T;
  gsl_rng *gsl_r;

  // PRC2 / transcription parameters
  double noisy_methylate;
  double me0_me1, me1_me2, me2_me3;
  double me2factor, me3factor;
  double firingRateMax, firingRateMin, transcription_demethylate;
  double activation;
  
  // cell cycle parameters
  double firingFactor;
  double cellCycleDuration, G2duration;
  int cellCycles, cellCycleCount;

  // run parameters
  unsigned long loci, reactCount, maxReact;
  unsigned long samples, sampleFreq, sampleCount;
  unsigned long optimSteps;
  logical DNAreplication, resultsLastHourOnly, SILAC, resultsFinalLocus;
} parameters;
 
// specific function pointer typedef
typedef void (*func_ptr_t)( chromatin *, parameters *, flags *, int );

typedef struct {
  D_VEC *propensity;
  I_VEC *doReactionParam;
  I_VEC *methylate_index;
  I_VEC *transcribeDNA_index;
  func_ptr_t *doReaction;
  double t_nextRep, t_nextEndG2;
  flags *update;
} gillespie;

typedef struct {
  I_MAT *K27;
  I_VEC *firing;
  D_VEC *t, *t_out;
  double tMax;
  unsigned long t_outLastSample;
} record;

