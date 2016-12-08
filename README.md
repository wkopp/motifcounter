# motifcounter - R package for analysing TFBSs in DNA sequences.

[![Travis-CI Build Status](https://travis-ci.org/wkopp/motifcounter.svg?branch=master)](https://travis-ci.org/wkopp/motifcounter)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/wkopp/motifcounter?branch=master&svg=true)](https://ci.appveyor.com/project/wkopp/motifcounter)
[![Coverage Status](https://img.shields.io/codecov/c/github/wkopp/motifcounter/master.svg)](https://codecov.io/github/wkopp/motifcounter?branch=master)

This software package grew out of the work that I did to obtain my PhD.

## Purpose of the `motifcounter` package

The `motifcounter` package facilitates the analysis of
 transcription factor binding sites (TFBSs) in DNA sequences.
It can be used to scan a set of DNA sequences for known motifs 
(e.g. from TRANSFAC or JASPAR) in order to determine the positions
and enrichment of TFBSs in the sequences.

Towards the end of motif hit enrichment, we developed
two analytic approximations of the *distribution of the number of motif hits*
which are the basis for the enrichment test:

1. A *compound Poisson approximation* 
2. A *combinatorial approximation*

Both methods are threshold-based. 
They require the user to chose a
 desired false positive probability for obtaining motif hits, based
on which the threshold is determined.

While the *compound Poisson approximation* performs
most accurately with stringent false positive probabilities 
(due to "rare hit assumption"),
the *combinatorial approximation* also achieves highly accurate results
with relaxed false positive probabilities.

Both approximations take higher-order sequence feature in the background
model into account. To this end, the user provides a set of sequences
and the desired order *d* which are used to estimate an **order-*d* 
background model**.

Moreover, both approximations take the **self-overlapping structure**
of the motifs into account, which are prevalent in repeat-like or palindromic
motif hits. Such motifs induce **motif clumps** (that is, mutually
overlapping motif hits) when a DNA sequence is scanned for hits.
As a consequence of **motif clumping**, the distribution of the number of
motif hits, and thus, the enrichment test are affected.
To this end, both methods take motif hits on both strands into account. 

## Installation
An easy way to install `motifcounter` is to facilitate the `devtools` R package.

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
aligned sequences, and 3. to test for motif 
enrichment in a given set of sequences.

The vignette can be found as follows:

```R
library(motifcounter)
vignette("motifcounter")
```
