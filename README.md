
# DirichletMultinomial

Dirichlet-multinomial mixture models can be used to describe
variability in microbial metagenomic data. This package is an
interface to code originally made available by Holmes, Harris, and
Quince, 2012, PLoS ONE 7(2): 1-15, as discussed further in the man
page for this package, ?DirichletMultinomial.

## Installation

Install [DirichletMultinomial][] from [Bioconductor][] with:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DirichletMultinomial")
```

Linux and MacOS with source-level installations require the 'gsl'
system dependency. On Debian or Ubuntu

```
sudo apt-get install -y libgsl-dev
```

On Fedora, CentOS or RHEL

```
sudo yum install libgsl-devel
```

On OS-X

```
brew install gsl
```

## Use

See the [DirichletMultinomial][] landing page and [vignette][] for use.

[DirichletMultinomial]: https://bioconductor.org/packages/DirichletMultinomial
[Bioconductor]: https://bioconductor.org
[vignette]: https://bioconductor.org/packages/release/bioc/vignettes/DirichletMultinomial/inst/doc/DirichletMultinomial.pdf




