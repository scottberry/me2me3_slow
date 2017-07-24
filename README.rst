**********
me2me3_slow
**********

Code for performing Gillespie algorithm simulations of chromatin as described in:

Berry, S., Dean, C., and Howard, M. (2017). Slow Chromatin Dynamics Allow Polycomb Target Genes to Filter Fluctuations in Transcription Factor Activity. Cell Systems 4, 445–457.e448.

Disclaimer
----------

The codes were not developed to be particularly user-friendly, so are provided in the interests of reproducibility of research and to stimulate further development in computational simulations of models of chromatin/transcription.

Program structure
==========

Functions with captial letters (including Main.c) are executables (compiled using ```make```) used for various different simulations. These all have common command line options as defined in parse.c:

```
$ ./me2me3
Usage:
 -c <control region>
 -a <alpha>
 -b <beta>
 -t <firing threshold>
 -i <identifier>
 -s (seed based on identifier)
 -r (DNA replication ON)
 -m (start in K27me3 state)
 -u (start in unmodified state)
 -n <translation efficiency (noise)>
 -h <histone turnover rate (per histone per transcription)>
```

In general, other parameter values such as methylation/demethylation rates will need to be set at compile-time. Main.c can also be setup to scan a range of such parameter values.

Submission scripts
=========

Example bash scripts for submitting these compiled functions, and looping over various command-line parameter inputs are provided in the ```scripts``` directory.

Library functions
=========

The code makes use of the ```scottsmatrices.c``` library, which provides get/free inferace functions (e.g. ```i_vec_get()``` / ```i_vec_free()```, ```d_mat_get()``` / ```d_mat_free()```) for memory allocation/deallocation . The allocated struct returned by these interface functions contains also information about their size, and type and functions exist in ```scottsmatrices.c``` to print matrices.

Dependencies
========

```me2me3_slow``` makes use of GNU scientific library for random number generation. It was developed primarily on Mac OS X and generalised to run on a cluster running the CentOS 6 unix distribution.
