# DirichletMultinomial

Dirichlet-multinomial mixture models can be used to describe variability
in microbial metagenomic data. This package is an interface to code
originally made available by Holmes, Harris, and Quince, 2012, PLoS ONE
7(2): 1-15, as discussed further in the man page for this package,
?DirichletMultinomial.

## Installation

Install
[DirichletMultinomial](https://bioconductor.org/packages/DirichletMultinomial)
from [Bioconductor](https://bioconductor.org) with:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DirichletMultinomial")
```

Linux and MacOS with source-level installations require the ‘gsl’ system
dependency. On Debian or Ubuntu

    sudo apt-get install -y libgsl-dev

On Fedora, CentOS or RHEL

    sudo yum install libgsl-devel

On macOS (source installations are not common on macOS, so this step may
is not usually necessary)

    brew install gsl

## Use

See the
[DirichletMultinomial](https://bioconductor.org/packages/DirichletMultinomial)
Bioconductor landing page and
[vignette](https://mtmorgan.github.io/DirichletMultinomial/articles/DirichletMultinomial.html)
for use.
