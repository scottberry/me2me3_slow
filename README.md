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

## Dynamic modulation (Figure 6)

## Cis memory window (Figure 7) 