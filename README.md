# motifcounter - R package for analysing TFBSs in DNA sequences.

[![Travis-CI Build Status](https://travis-ci.org/wkopp/motifcounter.svg?branch=master)](https://travis-ci.org/wkopp/motifcounter)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/wkopp/motifcounter?branch=master&svg=true)](https://ci.appveyor.com/project/wkopp/motifcounter)
[![Coverage Status](https://img.shields.io/codecov/c/github/wkopp/motifcounter/master.svg)](https://codecov.io/github/wkopp/motifcounter?branch=master)

This software package grew out of the work that I did to obtain my PhD.

If it is of help for your analysis, please cite

```
  @Manual{,
    title = {motifcounter: R package for analysing TFBSs in DNA sequences},
    author = {Wolfgang Kopp},
    year = {2016},
    note = {R package version 0.99.0},
  }
```

## Hallmarks of `motifcounter`

The `motifcounter` package facilitates the analysis of
 transcription factor binding sites (TFBSs) in DNA sequences.
It can be used to scan a set of DNA sequences for known motifs 
(e.g. from TRANSFAC or JASPAR) in order to determine the positions
and enrichment of TFBSs in the sequences.

A major part of the package addresses the improvement and development
of novel algorithms for determining motif hit enrichment,
which is based on two analytic approximations of the 
*distribution of the number of motif hits* that are expected 
in random DNA sequences:

1. A *compound Poisson approximation* 
2. A *combinatorial approximation*

Both methods are threshold-based. 
That is, they require the user to chose a
 desired false positive probability for obtaining motif hits 
in random DNA sequences.

Both approximations yield highly accurate results for stringent
false positive levels, while, the novel *combinatorial approximation*
also yields highly accurate results for rather relaxed false positive levels.

Both approximations take higher-order sequence feature 
 in the background model into account (e.g. as found in CpG islands). 
To this end, `motifcounter` supports the
use of **order-*d* Markov models** for the background models, where
the desired order *d* is prescribed by the user.

Both approximations take the **self-overlapping structure**
of the motifs into account<sup>[1](#footnote)</sup>. 
Examples of which are repeat-like or palindromic
motifs. This facilitates accurate *P-value* estimates
for the enrichment test, 
regardless of the structure of the motif.

Lastly,  both approxmiations can be used to determine the number
of motif hits in both DNA strands, as they also take overlapping
hits on the opposite strand into account.
Optionally, the compound Poisson approximation can also be used
to determine the distribution of the number of motif hits
with respect to scanning only a single strand
which might be useful for studying e.g motif hit enrichment in RNA sequences.
This option will be available also for the combinatorial model in the future.

## Installation
An easy way to install `motifcounter` is by facilitating 
the `devtools` R package.

```R
#install.packages("devtools")
library(devtools)
install_github("wkopp/motifcounter")
```

Alternatively, the package can also be downloaded from this github-rep 
and installed via the `R CMD INSTALL` command.

## Getting started

The `motifcounter` package provides a vignette that walks you through
the analysis of some toy examples, including
1. to determine position- and strand-specific TF motif binding sites,
2. to analyse the profile of motif hit occurrences across a set of 
aligned sequences, and 
3. to test for motif enrichment in a given set of sequences.

The vignette can be found as follows:

```R
library(motifcounter)
vignette("motifcounter")
```

<a name="footnote1">1</>:Self-overlapping motifs induce 
**clumps of motif hits** (that is, mutually
overlapping motif hits) when a DNA sequence is scanned for hits.
As a consequence of **motif clumping**, the distribution of the number of
motif hits, and thus, the enrichment test are affected.
