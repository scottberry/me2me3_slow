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

#define U 0
#define UR 1
#define M 2
#define MR 3

typedef unsigned char logical;

typedef struct {
  unsigned long sites;
  I_VEC *state;
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
  double noisy_Rep_ON, noisy_UR_Rep_OFF, noisy_MR_Rep_OFF, noisy_demethylate;
  double UR_methylate, MR_methylate;
  double firingRateMax, firingRateMin, transcription_RepOFF, transcription_demethylate;

  logical results, testProb;
  unsigned long samples, sampleFreq, sampleCount;
} parameters;
 
// specific function pointer typedef
typedef void (*func_ptr_t)( chromatin *, parameters *, flags *, int );

typedef struct {
  D_VEC *propensity;
  I_VEC *doReactionParam;
  I_VEC *bindRep_index, *unbindRep_index, *methylate_index, *demethylate_index;
  I_VEC *transcribeDNA_index;
  func_ptr_t *doReaction;
  flags *update;
} gillespie;

typedef struct {
  I_MAT *state;
  I_VEC *firing;
  I_VEC *events;
  D_VEC *t, *t_out;
} record;

