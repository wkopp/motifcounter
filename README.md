# mdist - R package for determining motif enrichment in DNA sequences.

The package implements functions for statistical enrichment tests to
decide whether a motif (e.g. Transcription factor PWM from TRANSFAC) is
statistically overrepresented in a sequence of interest (e.g. a gene promoter).

The method addresses two important aspects: Firstly, it allows for the usage
of higher-order Markov models as background. This will in particular improve 
enrichment tests i.e. in CpG islands.
Secondly, we take the self-overlapping structure of motifs into account
 which is known to affect the motif hit distribution especially for palindromic
 (e.g. PPARG) or highly repetitive motifs (e.g. SP1). We provide
 two methods for determining the number of motif hits distribution, which
 make use of the compound Poisson model and a Bayesian model, respectively.

## Installation
The `mdist` package can be installed by cloning the git repository and typing

```R
R CMD INSTALL mdist
```
or via the `devtools` package

```R
install.packages("devtools")
library(devtools)
install_github("wkopp/mdist")
```
