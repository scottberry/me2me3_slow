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

#define LIGHT 0
#define HEAVY 1
#define UNLABELLED 2

typedef unsigned char logical;

typedef struct {
  long sites;
  long controlSites;
  I_VEC *K27;
  I_VEC *silac;
} chromatin;

typedef struct {
  logical histone, transcribed;
} flags;

typedef struct {
  // random number parameters
  const gsl_rng_type *gsl_T;
  gsl_rng *gsl_r;

  // PRC2 / transcription parameters
  double noisy_me0_me1, noisy_me1_me2, noisy_me2_me3;
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
  logical DNAreplication, resultsLastHourOnly, resultsFinalLocus;

  // silac parameters
  logical silacExperiment;
  int silacLabel;
  long silacLightCycles;
  long silacHeavyCycles;
  double SILAC_0h, SILAC_10h, SILAC_24h, SILAC_48h, SILAC_nextReport;
  int SILAC_report;
  
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

  // test parameters
  FILE *test_fptr;
  logical test;

} gillespie;

typedef struct {
  I_MAT *K27;
  I_MAT *silac;
  I_VEC *firing;
  D_VEC *t, *t_out;
  double tMax;
  unsigned long t_outLastSample;
} record;

/* Function prototypes */

// random.c
double runif(gsl_rng *r);
void rseed(parameters *p);
void setseed(parameters *p, long seed);
void rfree(parameters *p);

// modifications.c
void methylate(chromatin *c, parameters *p, flags *update, int pos);
void demethylate(chromatin *c, parameters *p, flags *update, int pos);
void transcribeDNA(chromatin *c, parameters *p, flags *update, int pos);
void replicateDNA(chromatin *c, parameters *p, flags *update);

// gillespie.c
void initialiseRepressed(chromatin *c);
void initialiseActive(chromatin *c);
void initialiseSilacLight(chromatin *c);
void initialiseRandom(chromatin *c, parameters *p);
double d_vec_sum(D_VEC *d);
void initialiseGillespieFunctions(chromatin *c, gillespie *g);
double frac(I_VEC *vec, int target);
double fracControlRegion_me2me3(chromatin *c);
double enzymaticFactor(chromatin *c, parameters *p, int pos);
double neighboursK27factor(chromatin *c, parameters *p, int pos);
void updatePropensities(chromatin *c, parameters *p, gillespie *g);
double gillespieTimeStep(parameters *p, gillespie *g, double *p_s);
void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r);

// results.c
char *str_replace(char *orig, char *rep, char *with);
void fprint_t_out_nCycles(char *fname, record *r);
void fprint_t_nCycles(char *fname, I_MAT *mat, int target, record *r);
void fprint_silac_t_nCycles(char *fname, I_MAT *mat, int target, I_MAT *silac, int silac_target, record *r);
void fprint_firing_t_nCycles(char *fname, record *r);
double tAverageGap_nCycles(chromatin *c, parameters *p, record *r);
double tAverageGap_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double prob_me2_me3_nCycles(chromatin *c, parameters *p, record *r);
double prob_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double prob_me0_me1_nCycles(chromatin *c, parameters *p, record *r);
double prob_me0_me1_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double tAverage_me2_me3_nCycles(chromatin *c, parameters *p, record *r);
double tAverage_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r);
unsigned long numberHistoneStateFlips(record *r);
double firstPassageTime(record *r, signed char *initial);
int writelog(FILE *fptr, chromatin *c, parameters *p, record *r);
void fprintTripleSILAC(FILE *fptrAbs, FILE *fptrRel, long locus, parameters *p, record *r);
