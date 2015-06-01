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
  I_VEC *K27;
} chromatin;

typedef struct {
  logical protein, histone, transcribed;
} flags;

typedef struct {
  const gsl_rng_type *gsl_T;
  gsl_rng *gsl_r;
  unsigned long reactCount, maxReact;
  unsigned long optimSteps;
  unsigned long loci;
  
  double noisy_methylate;
  double me0_me1, me1_me2, me2_me3;
  double me2factor, me3factor;
  double firingRateMax, firingRateMin, transcription_demethylate;
  double firingFactor;

  double cellCycleDuration, G2duration;
  int cellCycles, cellCycleCount;
  
  logical testProb, DNAreplication;
  unsigned long samples, sampleFreq, sampleCount;
} parameters;
 
// specific function pointer typedef
typedef void (*func_ptr_t)( chromatin *, parameters *, flags *, int );

typedef struct {
  D_VEC *propensity;
  I_VEC *doReactionParam;
  I_VEC *methylate_index;
  I_VEC *transcribeDNA_index;
  func_ptr_t *doReaction;
  flags *update;
} gillespie;

typedef struct {
  I_MAT *K27;
  I_VEC *firing;
  I_VEC *events;
  D_VEC *t, *t_out;
  double tMax;
  unsigned long t_outLastSample;
} record;

