# PRC2-based epigenetic silencing directly antagonised by transcription

## Lack of bistability in two-state model (Figure S1)

**Plots:** PlotParameterSpaceTwoState.R

## Bistability in non-processive model (Figures 2A, 2B, S2)

Submit independent simulations with values of histone turnover [0.00025, 0.008]

Figure 2A,2B: histone exchange probability = 0.001

**Submission command:**
`./HistoneTurnover -h{exchange} -r -t1.0`

**Compile-time options flags:**
p.resultsLastHourOnly = TRUE;
p.resultsFinalLocus = FALSE;
p.checkHistoneTurnover = TRUE;
p.stochasticAlpha = FALSE;
p.countFiringEvents = FALSE;

**Other parameters:**
sites: 60
controlSites: 60
loci: 200
maxReact: 200000
samples: 200000
optimSteps: 30
DNAreplication: TRUE
resultsLastHour: TRUE
SILAC: FALSE
cellCycleDuration: 22.00 hours
cellCycles: 50
me2_me3: [0.0032,0.000004]
firingRateMax: 0.0001*2^[1,6]
firingRateMin: 0.0001000000
firingThreshold: 1.0000000000
stochasticAlpha: FALSE
transcription_demethylate: [0.25,0.00032]

Simulation time: 24783.230469 seconds

**Plots:** PlotParameterSpaceTurnover.R

## Histone exchange validation (Figure 3)

**Submission commands:**

`./HistoneTurnover -h0.001 -r -t1.0 -m -im`  
`./HistoneTurnover -h0.001 -r -t1.0 -u -iu`

**Compile-time options flags:**

p.resultsLastHourOnly = TRUE;
p.resultsFinalLocus = TRUE;
p.checkHistoneTurnover = TRUE;
p.stochasticAlpha = FALSE;
p.countFiringEvents = FALSE;

**Other parameters:**
Program name: ./HistoneTurnover
Operating system: Mac OS
sites: 60
controlSites: 60
loci: 1
maxReact: 200000
samples: 200000
optimSteps: 1
DNAreplication: TRUE
resultsLastHour: TRUE
SILAC: FALSE
cellCycleDuration: 22.00 hours
cellCycles: 50
me2_me3: 0.0000300000
firingRateMax: 0.0040000000
firingRateMin: 0.0001000000
firingThreshold: 1.0000000000
stochasticAlpha: FALSE
alpha: 1.0000
beta: 1.0000
transcription_demethylate: 0.0200000000
transcription_turnover: 0.0010000000
noisy_demethylate: 0.0000020000
Simulation time: 0.137020 seconds

**Plots:**
HistoneTurnoverExampleLoci.R (3A/3B)
TurnoverSpecies.R (3C)

## The effect of processivity in methylation or demethylation (Figure S3)

Processivity in methylation and demethylation is incorporated through variant functions methylate() and demethylate() in the files nonprocessive.c, processiveDemethylation.c, processiveMethylation.c. These processive functions are automatically linked by the executables ProcMeth and ProcDemeth, which are compiled with Main.c. Similar parameter searches to the non-processive model can be performed from Main.c using the parameters and compile-time flags listed above.

**Submission commands:**

`./ProcMeth -h{exchange} -r -t1.0 -i Proc_Meth`  
`./ProcDemeth -h{exchange} -r -t1.0 -i Proc_Demeth`

**Plots:**
PlotParameterSpaceProcessive.R

## Fitting SILAC data (Figures 4, S4)

### Silac parameter search

**Compile-time options flags:**
p.silacExperiment = TRUE;
p.silacLightCycles = 5;
p.silacHeavyCycles = 1;
p.resultsLastHourOnly = TRUE;
p.resultsFinalLocus = FALSE;
p.checkHistoneTurnover = FALSE;

**Submission script:** bash script silac.sh peforms simulations initialised in the repressed state to calculate the fit to the Silac data for threshold [0.2,1.0]. The script then also submits simulations initialised in either the active or repressed state randomly (labelled _bal.txt for 'balanced').

**Other parameters:**
loci: 100
maxReact: 30000
samples: 30000
optimSteps: 30
cellCycles: 20

**Plots:**

### Silac example loci

All 1000 loci are individually simulated and processed using an R script. The first 20 are over-plotted and then all 1000 are averaged and plotted together with the experimental data

**Submission script:** individualsSilac.sh
contains the command `./Silac -r -m -s -i $i -t1 -h0.001 > $i.out 2>&1 &`

**Compile-time options flags:**
p.silacExperiment = TRUE;
p.silacLightCycles = 5;
p.silacHeavyCycles = 1;
p.resultsLastHourOnly = TRUE;
p.resultsFinalLocus = TRUE;
p.checkHistoneTurnover = FALSE;

**Other parameters:**
'Fast'
loci: 1
maxReact: 30000
samples: 30000
optimSteps: 1
cellCycles: 20
me2_me3: 0.00003
firingThreshold: 1.0
transcription_demethylate: 0.02
transcription_turnover: 0.001

'Slow'
loci: 1
maxReact: 30000
samples: 30000
optimSteps: 1
cellCycles: 20
me2_me3: 0.000008
firingThreshold: 0.333
transcription_demethylate: 0.004
transcription_turnover: 0.001

**Plots:**
SILAC_Individuals.R

## Noisy inputs (Figures 5, S6, S7)

The main executable Main.c (compiled as me2me3) is used for simulations with different noisy input signals. The noise level is controlled by the command line inputs -n (integer > 1, representing the number of proteins produced by a single RNA molecule).

### Variable noise parameter search

**Submission script:** noisy.sh
contains the command
`./me2me3 -t0.333 -r -h0.001 -n $noisy -i $noisy > noisy$noisy.txt 2>&1 &`

**Compile-time options flags:**
p.resultsLastHourOnly = TRUE;
p.resultsFinalLocus = FALSE;
p.checkHistoneTurnover = FALSE;
p.stochasticAlpha = TRUE;

**Other parameters:**
loci: 200
maxReact: 1000000
samples: 1000000
optimSteps: 20
cellCycles: 20
me2_me3: [0.003,0.000004]
firingThreshold: 0.333
stochasticTranslationEfficiency: [1,2000]
transcription_demethylate: [0.25,0.0003]
transcription_turnover: 0.001

**Plots**
NoisyParameterSpaceCV.R

### Individual examples

Recompile for fast and slow time-scales according to the parameters shown below.

**Submission commands:**
`./me2me3 -t0.333 -r -h0.001 -n1 -i slow_n1_m -m > out.txt 2>&1 &`  
`./me2me3 -t0.333 -r -h0.001 -n1 -i slow_n1_u -u > out.txt 2>&1 &`  
`./me2me3 -t0.333 -r -h0.001 -n1000 -i slow_n1000_m -m > out.txt 2>&1 &`  
`./me2me3 -t0.333 -r -h0.001 -n1000 -i slow_n1000_u -u > out.txt 2>&1 &`

**Compile-time options flags:**
p.resultsLastHourOnly = TRUE;
p.resultsFinalLocus = TRUE;
p.checkHistoneTurnover = FALSE;
p.stochasticAlpha = TRUE;

**Other parameters:**
loci: 1
maxReact: 1000000
samples: 1000000
optimSteps: 1
cellCycles: 20
firingRateMax: 0.004
firingRateMin: 0.0001
firingThreshold: 0.333
transcription_turnover: 0.001

'Fast'
me2_me3: 0.0000447
transcription_demethylate: 0.178

'Slow'
me2_me3: 0.000008
transcription_demethylate: 0.004

**Plots**
NoisyExamples.R

## Dynamic modulation (Figure 6)

The executable Dynamic.c takes alpha and beta as inputs and equilibrated for 5 cell cycles at alpha = beta = 1, before changing alpha and beta to the values specified as command line options.

There are 4 options that need to be simulated:
- StartM_decreaseBeta
- StartU_increaseBeta
- StartM_increaseAlpha
- StartU_decreaseAlpha

Plots are created by simulating 20 loci for each parameter set as individuals. 20 of these are over-plotted. For certain values, average plots over 100 loci need to be performed (alpha or beta = 1/8 or 8).

**Submission script:** individualsDynamic.sh
contains commands such as
`./Dynamic -c 60 -u -t 0.4 -b 1.0 -a $act -i $i -s > out$act$i.out 2>&1 &`

**Compile-time options flags:**
p.resultsLastHourOnly = TRUE;
p.resultsFinalLocus = TRUE;
p.checkHistoneTurnover = FALSE;
p.stochasticAlpha = FALSE;
p.countFiringEvents = TRUE;

**Other parameters:**
maxReact: 200000
samples: 200000
optimSteps: 1
cellCycles: 16
initialCellCycles = 5;
me2_me3: 0.0000080000
firingRateMax: 0.0040000000
firingRateMin: 0.0001000000
firingThreshold: 0.3330000000
transcription_demethylate: 0.0080000000
transcription_turnover: 0.0010000000

**Plots:**
IndividualsSummary_DecreaseBeta.R
IndividualsSummary_IncreaseBeta.R
IndividualsSummary_IncreaseAlpha.R
IndividualsSummary_DecreaseAlpha.R

## Cis memory window (Figure 7) 

