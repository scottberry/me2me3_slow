# PRC2-based epigenetic silencing directly antagonised by transcription

## Lack of bistability in two-state model (Figure S1)
==

## Bistability in non-processive model (Figures 2A, 2B, S2)
==

To run simulations, use the parameter values listed below:
Submit independent simulations with values of histone turnover [0.00025, 0.008]

Figure 2A,2B: histone exchange probability = 0.001

**Submission command:** ./HistoneTurnover -h{exchange} -r -t1.0

Program name: ./HistoneTurnover
Operating system: Unix
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

## The effect of processivity in methylation or demethylation (Figure S3)

## Fitting SILAC data (Figures 4, S4)

## Noisy inputs (Figures 5, S6, S7)

## Dynamic modulation (Figure 6)

## Cis memory window (Figure 7) 