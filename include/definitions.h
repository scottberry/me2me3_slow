/* 
   Definitions file for Gillespie algorithm simulations of models
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

#define H3_1 1
#define H3_3 3

typedef unsigned char logical;

typedef struct {
  long sites;
  long controlSites;
  I_VEC *K27;
  I_VEC *silac;
  I_VEC *turnover;
  I_VEC *variant;
  logical promoterON;
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
  double noisy_demethylate;
  double me0_me1, me1_me2, me2_me3;
  double me2factor, me3factor;
  double firingRateMax, firingRateMin, transcription_demethylate, transcription_turnover;
  double firingThreshold, firingCap;
  double alpha, beta;

  // cell cycle parameters
  double firingFactor;
  double cellCycleDuration;
  int cellCycles, cellCycleCount;
  int initialCellCycles;
  
  // run parameters
  logical startM, startU, randomSeed;
  long seed;
  unsigned long loci, reactCount, maxReact;
  unsigned long samples, sampleFreq, sampleCount;
  unsigned long optimSteps;
  logical DNAreplication, resultsLastHourOnly, resultsFinalLocus, spatialResults;
  logical checkHistoneTurnover, countFiringEvents;
  char id[128];
  char executable[128];

  // silac parameters
  logical silacExperiment;
  int silacLabel;
  long silacLightCycles;
  long silacHeavyCycles;
  double SILAC_0h, SILAC_10h, SILAC_24h, SILAC_48h, SILAC_nextReport;
  int SILAC_report;
  logical resultsSilacEachLocus;

  // transFactor parameters
  logical stochasticAlpha;
  long stochasticTranslationEfficiency;
  long transFactorRNA;
  long transFactorProtein;
  double k_r, k_p, gamma_r, gamma_p;

  // bursty parameters
  logical burstyFiring;
  double k_onMax, k_onMin, k_off;
  double constFiring;
  
} parameters;
 
// specific function pointer typedef
typedef void (*func_ptr_t)( chromatin *, parameters *, flags *, int );

typedef struct {
  D_VEC *propensity;
  I_VEC *doReactionParam;
  I_VEC *methylate_index;
  I_VEC *demethylate_index;
  I_VEC *transcribeDNA_index;
  func_ptr_t *doReaction;
  double t_nextRep;
  flags *update;

  // transFactor indices
  I_VEC *decreaseProtein_index;
  I_VEC *increaseProtein_index;
  I_VEC *decreaseRNA_index;
  I_VEC *increaseRNA_index;

  // burstyFiring indices
  I_VEC *activatePromoter_index;
  I_VEC *deactivatePromoter_index;
  
  // test parameters
  FILE *test_fptr;
  logical test;

} gillespie;

typedef struct {
  I_MAT *K27;
  I_MAT *silac;
  I_MAT *variant;
  I_VEC *firing;
  I_VEC *promoterON;
  D_MAT *turnover;
  D_VEC *t, *t_out;
  double tMax;
  unsigned long t_outLastSample;

  D_VEC *silacResultsLight_0h, *silacResultsLight_10h, *silacResultsLight_24h, *silacResultsLight_48h;
  D_VEC *silacResultsHeavy_0h, *silacResultsHeavy_10h, *silacResultsHeavy_24h, *silacResultsHeavy_48h;

  I_VEC *transFactorProtein, *transFactorRNA;
  D_VEC *alpha;
  
} record;

typedef struct {
  signed char initial;
  double gap, Mavg, tTot, tTotM, tTotU, tM, tU, lifetime, me3_end;
  double firstPassage, firstPassageM, firstPassageU, fpU, fpM;
  long fh, initM, initU;
  double probM, probU, bistability;
  double totalHistoneTurnover;
  double alphaSD, alphaMean;
  double fracH3_3_M, fracH3_3_U, avgH3_3_M, avgH3_3_U;
  double burstFrequency, burstDuration, quiescentDuration;
  double simulationTime;
  long long bursts, firingEvents, totalFiringEvents;
} quantification;

/* Function prototypes */

// parse.c
void usage(void);
void parseCommandLine(int argc, char *const *argv, chromatin *c, parameters *p);

// random.c
double runif(gsl_rng *r);
void rseed(parameters *p);
void setseed(parameters *p, long seed);
void rfree(parameters *p);

// modifications.c
void transcribeDNA(chromatin *c, parameters *p, flags *update, int pos);
void replicateDNA(chromatin *c, parameters *p, flags *update);
void decreaseRNA(chromatin *c, parameters *p, flags *update, int pos);
void increaseRNA(chromatin *c, parameters *p, flags *update, int pos);
void decreaseProtein(chromatin *c, parameters *p, flags *update, int pos);
void increaseProtein(chromatin *c, parameters *p, flags *update, int pos);
void activatePromoter(chromatin *c, parameters *p, flags *update, int pos);
void deactivatePromoter(chromatin *c, parameters *p, flags *update, int pos);

// nonprocessive.c or processivedemethylation.c or processivemethylation.c
void methylate(chromatin *c, parameters *p, flags *update, int pos);
void demethylate(chromatin *c, parameters *p, flags *update, int pos);

// gillespie.c
void allocateGillespieMemory(chromatin *c, parameters *p, gillespie *g, record *r);
void freeGillespieMemory(chromatin *c, parameters *p, gillespie *g, record *r);
void initialiseRepressed(chromatin *c);
void initialiseActive(chromatin *c);
void initialiseSilacLight(chromatin *c);
void initialiseH3_1(chromatin *c);
void initialiseRandom(chromatin *c, parameters *p);
void initialiseMixed(chromatin *c, parameters *p);
void initialisePromoterOFF(chromatin *c);
double d_vec_sum(D_VEC *d);
void initialiseGillespieFunctions(chromatin *c, gillespie *g);
void initialiseGillespieFunctionsTransFactor(chromatin *c, gillespie *g);
void initialiseGillespieFunctionsBurstyFiring(chromatin *c, gillespie *g);
double frac(I_VEC *vec, int target);
double fracControlRegion_me2me3(chromatin *c);
double enzymaticFactor(chromatin *c, parameters *p, int pos);
double neighboursK27factor(chromatin *c, parameters *p, int pos);
double firingRate(parameters *p, double f_me2_me3);
double k_on(parameters *p, double f_me2_me3);
void updatePropensities(chromatin *c, parameters *p, gillespie *g);
void updatePropensitiesTransFactor(parameters *p, gillespie *g);
double gillespieTimeStep(parameters *p, gillespie *g, double *p_s);
void gillespieStep(chromatin *c, parameters *p, gillespie *g, record *r);

// results.c
void allocateSilacRecordMemory(chromatin *c, parameters *p, record *r);
void freeSilacRecordMemory(record *r);
void incrementSilacReportPoint(parameters *p);
void resetQuantification(quantification *q);
void accumulateQuantification(chromatin *c, parameters *p, record *r, quantification *q);
void averageQuantification(chromatin *c, parameters *p, record *r, quantification *q);
char *parameterDependentBasename(chromatin *c, parameters *p);
void fprintParameterSpaceHeader(FILE *parFile);
void fprintParameterSpaceResults(FILE *parFile, parameters *p, chromatin *c, quantification *q);
void fprintBurstyResultsHeader(FILE *parFile);
void fprintBurstyResults(FILE *parFile, parameters *p, chromatin *c, quantification *q);
char *str_replace(char *orig, char *rep, char *with);
void fprint_t_out_nCycles(char *fname, record *r);
void fprint_t_nCycles(char *fname, I_MAT *mat, int target, record *r);
void fprint_variant_t_nCycles(char *fname, record *r, int variant_target);
void fprint_silac_t_nCycles(char *fname, I_MAT *mat, int target, I_MAT *silac, int silac_target, record *r);
void fprint_silacRelative_t_nCycles(char *fname, I_MAT *mat, int target, I_MAT *silac, int silac_target, record *r);
void fprint_firing_t_nCycles(char *fname, record *r);
double tAverageGap_nCycles(chromatin *c, parameters *p, record *r);
double tAverageGap_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double tAverageVariant_lastHour_nCycles(chromatin *c, parameters *p, record *r, int variant_target);
double tAverageVariant_nCycles(chromatin *c, parameters *p, record *r, int variant_target);
double prob_me2_me3_nCycles(chromatin *c, parameters *p, record *r);
double prob_lowExpression_nCycles(chromatin *c, parameters *p, record *r);
double prob_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double prob_lowExpression_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double prob_me0_me1_nCycles(chromatin *c, parameters *p, record *r);
double prob_highExpression_nCycles(chromatin *c, parameters *p, record *r);
double prob_me0_me1_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double prob_highExpression_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double tAverage_me2_me3_nCycles(chromatin *c, parameters *p, record *r);
double tAverage_me2_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r);
double tAverage_me3_lastHour_nCycles(chromatin *c, parameters *p, record *r);
unsigned long numberHistoneStateFlips(record *r);
double firstPassageTime(record *r, signed char *initial);
double firstPassageTimeExpression(record *r, parameters *p, signed char *initial);
double meanBurstDuration(record *r);
double meanQuiescentDuration(record *r);
double meanBurstFrequency(record *r);
long countBursts(record *r);
long countFiringEvents(record *r);  
int writelog(FILE *fptr, chromatin *c, parameters *p, record *r);
void fprintTripleSILAC_eachLocus(FILE *fptrAbs, FILE *fptrRel, long locus, parameters *p, record *r);
void storeTripleSILAC_me3(long locus, parameters *p, record *r);
void fprintTripleSILAC_average(FILE *fptr, parameters *p, record *r);
void fprintHistoneTurnover(FILE *fptr, parameters *p, record *r);
void fprintResultsFinalLocus(char *avgfile, record *r);
void fprintVariantResultsFinalLocus(char *avgfile, record *r);
void fprintSilacResultsFinalLocus(char *avgfile, record *r);
void fprintSilacResultsRelativeFinalLocus(char *avgfile, record *r);
void fprint_transFactorProtein_nCycles(char *fname, record *r);
void fprint_promoterStatus_nCycles(char *fname, record *r);
void fprint_alphaOnly_nCycles(char *fname, record *r);
double tAverageAlpha(record *r);
double tAverageAlphaSD(record *r);
long countFiringEventsLastCellCycle(parameters *p, record *r);
void resetFiringRecord(record *r);
void fprintSpatialResults(FILE *fptr, record *r);
