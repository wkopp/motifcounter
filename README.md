# motifcounter - R package for determining motif enrichment in DNA sequences.

[![Travis-CI Build Status](https://travis-ci.org/wkopp/motifcounter.svg?branch=master)](https://travis-ci.org/wkopp/motifcounter)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/wkopp/motifcounter?branch=master&svg=true)](https://ci.appveyor.com/project/wkopp/motifcounter)


The package implements functions for  elucidating the statistical enrichment 
of known motif (e.g. from TRANSFAC or JASPAR) in a given DNA sequence of interest (e.g. a gene promoter).

To this end, we provide  two analytic approximations for computing the distribution of the number of motif hits
in a sequence of given length, which are referred to as the *compound Poisson approximation* and the *combinatorial approximation*.
While the *combinatorial approximation* achieves advantageous accuracy if relaxed thresholds for calling motif hits
are used, it comes at the cost of an algorithm whose runtime depends on the length of the DNA sequence that needs to be scanned.
On the other hand, the *compound Poisson approximation* achieves similarly accurate 
results compared to the *combinatorial approximation* when stringent thresholds for calling motif hits are used, however,
using an more efficient algorithm.

Both methods address three important aspects about motif hit enrichment analysis:
- Firstly, it allows for using sophisticated 
 higher-order Markov models as background models, 
 which are more adequate for capturing the content of DNA sequence (e.g. of promoters, enhancers or CpG islands) compared to order-zero background models. Consequently using higher-order background models over order-zero background models
might reduce false positive calls substantially.

- Secondly, the self-overlapping structure of motifs is taken into account
 which is known to affect the distribution of the number of motif hits. This is especially the case for palindromic
motifs (e.g. PPARG) and repetitive motif structures (e.g. SP1).

- Third, the methods take motif hits on both strands of the DNA into account

## Installation
The `motifcounter` package can be installed by cloning the git repository and typing

```R
R CMD INSTALL motifcounter
```
or via the `devtools` package

```R
install.packages("devtools")
library(devtools)
install_github("wkopp/motifcounter")
```

The package contains a vignette that explains the main functionality.
