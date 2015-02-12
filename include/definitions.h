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
#define M 1
#define MR 2
#define UR 3
#define M_LHP1 4
#define MR_LHP1 5

#define PI 3.1415926535897932

typedef unsigned char logical;

typedef struct {
  unsigned long sites;
  I_VEC *state;
} chromatin;

typedef struct {
  logical protein, histone;
} flags;

// specific function pointer typedef
typedef void (*func_ptr_t)( chromatin *, flags *, int );

typedef struct {
  const gsl_rng_type *gsl_T;
  gsl_rng *gsl_r;
  unsigned long reactCount, maxReact;
  double noisy_RepON, noisy_RepOFF, noisy_LHP1_OFF, noisy_demethylate;
  double stabilised_LHP1_OFF, M_LHP1_ON;
  double U_demethylate, U_LHP1_OFF;
  double UR_methylate, MR_methylate, M_bindRep, LHP1_bindRep;
  logical results, testProb;
  D_VEC *propensity;
  I_VEC *doReactionParam;
  I_VEC *bindRep_index, *unbindRep_index, *methylate_index, *demethylate_index;
  I_VEC *bindLHP1_index, *unbindLHP1_index;
  func_ptr_t *doReaction;
  flags *update;
  
  unsigned long samples, sampleFreq, sampleCount;
} parameters;

typedef struct {
  I_MAT *state;
  I_VEC *events;
  D_VEC *t, *t_out;
} record;

