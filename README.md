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
## Usage

```
# Estimate a background model on a set of sequences
bg <- readBackground(sequences, order)

# Normalize a given PFM
new_motif <- normalizeMotif(motif)

# Evaluate the scores along a given sequence
scores <- scoreSequence(sequence, motif, bg)

# Evaluate the motif hits along a given sequence
scores <- motifHits(sequence, motif, bg)

# Evaluate the average score profile
scores <- scoreProfile(sequences, motif, bg)

# Evaluate the average motif hit profile
scores <- motifHitProfile(sequences, motif, bg)

# Compute the motif hit enrichment
scores <- motifEnrichment(sequences, motif, bg)
```

## Hallmarks of `motifcounter`

The `motifcounter` package facilitates the analysis of
 transcription factor binding sites (TFBSs) in DNA sequences.
It can be used to scan a set of DNA sequences for known motifs
(e.g. from TRANSFAC or JASPAR) in order to determine the positions
and enrichment of TFBSs in the sequences.

Therefore, an analysis with `motifcounter` requires as input
1. a position frequency matrix (PFM) which represents the TF affinity towards the DNA
2. a background model, which is estimated from a given DNA sequence and which
serves as a reference for the statistical analysis.
3. a desired false positive level, for identifying putative TFBSs in DNA sequences. For example, a reasonable choice would be to choose a false positive level such that only one in 1000 positions are called TFBSs falsely.
4. a given DNA sequence, which is subject to the TFBS analysis.

The package aims to improve motif hit enrichment analysis. To this end,
the package offers a number of features:
1. `motifcounter` supports the use of **higher-order Markov models**
to account for the sequence composition in unbound DNA segments.
This improves the reliability of the enrichment analysis, because higher-order
sequence features occur commonly in natural DNA sequences (e.g. CpG islands).
2. The package automatically accounts for **self-overlapping** motif
structures<sup><a href="#fn1" id="ref1">1</a></sup>. This aspect is important
for reducing the false positives obtained from the enrichment test, which is
prevalent for repeat-like and palindromic motifs.
`motifcounter` not only determines self-overlapping motif hit occurrences
on a single DNA strand, but (by default)
also with respect to the reverse strand.

### Enrichment model
`motifcounter` implements two analytic approximations of the
*distribution of the number of motif hits*
in random DNA sequences that can optionally be used for the
enrichment test:

1. A *compound Poisson approximation*
2. A *combinatorial approximation*

Both approximations yield highly accurate results for stringent
false positive levels.
Moreover, if you intend to analyse long DNA sequences or
a large set of individual sequences (total sequence length >10kb),
we recommend to use the *compound Poisson approximation*.
On the other hand, we recommend the *combinatorial approximation*
if a relaxed false positive level is prefered to identify TFBSs.




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

The `motifcounter` package contains a tutorial that illustrates:
1. how to determine position- and strand-specific TF motif binding sites,
2. how to analyse the profile of motif hit occurrences across a set of
aligned sequences, and
3. how to test for motif enrichment in a given set of sequences.

The tutorial can be found in the package-vignette:

```R
library(motifcounter)
vignette("motifcounter")
```

## Acknowledgements
Thanks to [matthuska](https://github.com/matthuska) for reviewing and commenting
on the package.
<hr></hr>
<sup id="fn1">1: Self-overlapping motifs induce
**clumps of motif hits** (that is, mutually
overlapping motif hits) when a DNA sequence is scanned for hits.
As a consequence of **motif clumping**, the distribution of the number of
motif hits, and thus, the enrichment test are affected.<a href="#ref1" title="Jump back to footnote 1 in the text.">↩</a></sup>
