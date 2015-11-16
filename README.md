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
```
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
```

```c
for (p2=0;p2<p.optimSteps;p2++) {
	for (p3=0;p3<p.optimSteps;p3++) {
		P_DEMETHYLATE = pow(10,-0.05*(p2+8));
		P_METHYLATE = pow(10,-0.05*(p3+50));
	}
}
```


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
```
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
transcription_demethylate: 0.02  
transcription_turnover: 0.001  
noisy_demethylate: 0.000002  
Simulation time: 0.137020 seconds  
```
**Plots:**
HistoneTurnoverExampleLoci.R (3A/3B)  
TurnoverSpecies.R (3C)

## The effect of processivity in methylation or demethylation (Figure S3)

Processivity in methylation and demethylation is incorporated through variant functions methylate() and demethylate() in the files nonprocessive.c, processiveDemethylation.c, processiveMethylation.c. These processive functions are automatically linked by the executables ProcMeth and ProcDemeth, which are compiled with Main.c. Similar parameter searches to the non-processive model can be performed from Main.c using the parameters and compile-time flags listed above.

**Submission commands:**

`./ProcMeth -h{exchange} -r -t1.0 -i Proc_Meth`  
`./ProcDemeth -h{exchange} -r -t1.0 -i Proc_Demeth`

Note: need to recompile for ProcMeth as methylation rate range is divided by 10 to account for faster methylation of me0 substrates than me3 substrates.

**Other parameters:**
```
Program name: ./HistoneTurnover  
Operating system: Mac OS  
sites: 60  
controlSites: 60  
loci: 200  
maxReact: 200000  
samples: 200000  
optimSteps: 60  
DNAreplication: TRUE  
resultsLastHour: TRUE  
SILAC: FALSE  
cellCycleDuration: 22.00 hours  
cellCycles: 50  
me2_me3: 10^-0.05*([50,109]);
firingRateMax: 0.0001*pow(2,[1,6]);  
firingRateMin: 0.0001000000  
firingThreshold: 1.0 
stochasticAlpha: FALSE  
alpha: 1.0000  
beta: 1.0000  
transcription_demethylate: 10^-0.05*([8,67]);  
transcription_turnover: 0.001  
noisy_demethylate: 0.000002  
Simulation time: 0.137020 seconds  
```

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

**Submission script:** bash script silac.sh peforms simulations initialised in the repressed state to calculate the fit to the Silac data for threshold [0.2,1.0]. The script then also submits simulations with equal number of loci initialised in either the active or repressed state (labelled _bal.txt for 'balanced').

**Other parameters:**
```
loci: 100  
maxReact: 30000  
samples: 30000  
optimSteps: 60  
cellCycles: 20  
```

```c
for (p2=14;p2<p.optimSteps;p2++) {
	for (p3=24;p3<p.optimSteps;p3++) {
		P_DEMETHYLATE = pow(10,-0.05*(p2+8));
		P_METHYLATE = pow(10,-0.05*(p3+50));
	}
}
```

**Plots:**
`SILAC_FitParameterSpace.R`  
`SILAC_Individuals_Thresh1_00.R`  
`SILAC_Individuals_Thresh0_33.R`

### Silac example loci

All 1000 loci are individually simulated and processed using an R script. The first 20 are over-plotted to show stochastic variation and then all 1000 are averaged and plotted together with the experimental data.

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
```
loci: 1  
maxReact: 30000  
samples: 30000  
optimSteps: 1  
cellCycles: 20  
me2_me3: 0.00003  
firingThreshold: 1.0  
transcription_demethylate: 0.02  
transcription_turnover: 0.001  
```
'Slow'  
```
loci: 1  
maxReact: 30000  
samples: 30000  
optimSteps: 1  
cellCycles: 20  
me2_me3: 0.000008  
firingThreshold: 0.333  
transcription_demethylate: 0.004  
transcription_turnover: 0.001  
```

**Plots:**
SILAC_Individuals.R

## Noisy inputs (Figures 5, S6, S7)

The main executable Main.c (compiled as me2me3) is used for simulations with different noisy input signals. The noise level is controlled by the command line inputs -n (integer > 1, representing the number of proteins produced by a single RNA molecule).

### Variable noise parameter search

#### Figure S6
Parameter search over several values of methylation and demethylation rates to simulate a few noise levels 

**Submission script:** noisy.sh
contains the command
`./me2me3 -t0.333 -r -h0.001 -n $noisy -i $noisy > noisy$noisy.txt 2>&1 &`

**Compile-time options flags:**
p.resultsLastHourOnly = TRUE;  
p.resultsFinalLocus = FALSE;  
p.checkHistoneTurnover = FALSE;  
p.stochasticAlpha = TRUE;  

**Other parameters:**
```
loci: 200  
maxReact: 1000000  
samples: 1000000  
optimSteps: 20  
cellCycles: 25  
me2_me3: [0.003,0.000004]
firingRateMax: 0.004  
firingRateMin: 0.0001  
firingThreshold: 0.333  
stochasticTranslationEfficiency: [1,2000]  
transcription_demethylate: [0.25,0.0003]  
transcription_turnover: 0.001   
```

Started 100 loci each on howard16a and howard16c (estimated 3 hours)

```c
optimSteps: 18
P_DEMETHYLATE = pow(10,-0.15*(p2+4));
P_METHYLATE = pow(10,-0.15*(p3+19));
```

**Plots**
NoisyParameterSpaceCV.R

#### Figure 5A, 5B
Choose several methylation rates logarithmically spaced over range of bistable models (with P_T=0.333) i.e. 4 above 8 x 10^-6 and three below 8 x 10^-6. Then select demethylation rate based on maximum FP. Simulate a large number of loci at high resolution for many different noise strengths.

(methylation, demethylation)

Faster
(0.000128,0.3)  
(0.000064,0.1)  
(0.000032,0.07)  
(0.000016,0.03)  

Fitted
(0.000008,0.004)  

Slower
(0.0000056,0.001)  
(0.000004,0.001)  

Noise levels:
n = 1, SD = 0.04
n = 10, SD = 0.09
n = 100, SD = 0.28
n= 200, SD = 0.4
n = 500, SD = 0.6
n = 1000, SD = 0.9

R script to interpolate these values (CorrelateNoiseWithSD.R) gives rise to the ~ squared relationship between SD and n. Use the following noise values to yield  evenly spaced points on SD-axis.

1    2    9   23   43   71  106  149  200  259  327  404  489  583  687  799 922 1053 1195 1346

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
```
loci: 1  
maxReact: 1000000  
samples: 1000000  
optimSteps: 1  
cellCycles: 20  
firingRateMax: 0.004  
firingRateMin: 0.0001  
firingThreshold: 0.333  
transcription_turnover: 0.001  
```

'Fast'  
```
me2_me3: 0.0000447  
transcription_demethylate: 0.178  
```

'Slow'  
```
me2_me3: 0.000008  
transcription_demethylate: 0.004  
```

**Plots**
NoisyExamples.R

## Dynamic modulation (Figure 6)

The executable Dynamic.c takes alpha and beta as inputs and equilibrated for 5 cell cycles at alpha = beta = 1, before changing alpha and beta to the values specified as command line options.

There are 4 options that need to be simulated:

* StartM_decreaseBeta  
* StartU_increaseBeta  
* StartM_increaseAlpha  
* StartU_decreaseAlpha  

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
```
maxReact: 200000  
samples: 200000  
optimSteps: 1  
cellCycles: 16  
initialCellCycles = 5;  
me2_me3: 0.000008  
firingRateMax: 0.002  
firingRateMin: 0.0001  
firingThreshold: 0.333  
transcription_demethylate: 0.008  
transcription_turnover: 0.001  
```

**Plots:**
IndividualsSummary_DecreaseBeta.R  
IndividualsSummary_IncreaseBeta.R  
IndividualsSummary_IncreaseAlpha.R  
IndividualsSummary_DecreaseAlpha.R  

## Cis memory window (Figure 7) 

### Transcriptional response (Figure 7A)
The executable Dynamic.c takes alpha and beta as inputs and equilibrated for 5 cell cycles at alpha = beta = 1, before changing alpha and beta to the values specified as command line options.

Calculated the average number of transcriptional firing events in the final cell cycle as a function of alpha for both uniform initial states.

**Submission script:** hysteresis.sh
Simulates 81 instances of the simulation, each for a different value of alpha (logarithmically spaced.

```
act=$(echo "scale=6;100*e(-0.05*$i*l(10))" | bc -l)
./Dynamic -r -u -h0.001 -t0.333 -b1.0 -a$act -i startU$i > out$i.out 2>&1 &
```

**Compile-time options flags:**
p.resultsLastHourOnly = TRUE;  
p.resultsFinalLocus = TRUE;  
p.checkHistoneTurnover = FALSE;  
p.stochasticAlpha = FALSE;  
p.countFiringEvents = TRUE;  

**Other parameters:**
```
loci: 1000  
maxReact: 400000  
samples: 400000  
optimSteps: 1  
cellCycles: 25  
me2_me3: 0.000008  
me2factor: 0.1000  
me3factor: 1.0000  
firingRateMax: 0.002  
firingRateMin: 0.0001  
firingThreshold: 0.333  
transcription_demethylate: 0.008  
```

**Plots:**
HysteresisPlots.R

### First passage time (Figure 7B)

Customised Main.c to generate me2me3 with sweep over alpha within the main function  `p.alpha = 100.0*pow(10,-0.05*p2);`. In this case there is no pre-equilibration and alpha remains fixed throughout.

Executed 32 instances in parallel and averaged results for plotting in R script.

**Other parameters:**
```
loci: 20  
maxReact: 10000000  
samples: 10000000  
optimSteps: 81  
cellCycles: 1000  
me2_me3: 0.000008  
firingRateMax: 0.002  
firingRateMin: 0.0001  
firingThreshold: 0.3333  
transcription_demethylate: 0.008  
transcription_turnover: 0.001  
```

**Plots:**
FixedAlphaParamSpace.R


