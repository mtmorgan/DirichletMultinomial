# Dirichlet-Multinomial generative classifiers.

Fit Dirichlet-Multinomial generative classifiers to groups (rows) within
a sample x taxon count matrix.

## Usage

``` r
dmngroup(count, group, k, ..., simplify = TRUE,
    .lapply = parallel::mclapply)
```

## Arguments

- count:

  [`matrix()`](https://rdrr.io/r/base/matrix.html) of sample x taxon
  counts.

- group:

  [`factor()`](https://rdrr.io/r/base/factor.html) or vector to be
  coerced to a factor, with as many elements as there are rows in
  `count`, indicating the group to which the corresponding sample
  belongs.

- k:

  [`integer()`](https://rdrr.io/r/base/integer.html), the number(s) of
  Dirichlet components to fit.

- ...:

  Additional arguments, passed to
  [`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md).

- simplify:

  Return only the best-fit model for each group?

- .lapply:

  An `lapply`-like function for application of group x k fits.

## Details

This function divided `count` into groups defined by `group`, creates
all combinations of `group` x `k`, and evaluates each using
[`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md).
When `simplify=TRUE`, the best (Laplace) fit is selected for each group.

## Value

An object of class `dmngroup`, a list of fitted models of class
[`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md).
When `simplify=TRUE`, elements are named by the group to which they
correspond.

## References

Holmes I, Harris K, Quince C, 2012 Dirichlet Multinomial Mixtures:
Generative Models for Microbial Metagenomics. PLoS ONE 7(2): e30126.
doi:10.1371/journal.pone.0030126.

## Author

Martin Morgan
[mailto:mtmorgan.xyz@gmail.com](mailto:mtmorgan.xyz@gmail.com)

## See also

[`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md),
[DirichletMultinomial-package](https://mtmorgan.github.io/DirichletMultinomial/reference/DirichletMultinomial-package.md),
[`vignette("DirichletMultinomial")`](https://mtmorgan.github.io/DirichletMultinomial/articles/DirichletMultinomial.md)

## Examples

``` r
## best fit for groups 'Lean' and 'Obese'; full example in vignette.
if (FALSE) bestgrp <- dmngroup(count, pheno, k=1:5, verbose=TRUE, 
                    mc.preschedule=FALSE)
 # \dontrun{}
data(bestgrp)
bestgrp
#> class: DMNGroup 
#> summary:
#>       k samples taxa       NLE   LogDet   Laplace       BIC       AIC
#> Lean  1      61  130  9065.657 162.3513  9027.371  9332.864  9195.657
#> Obese 3     193  130 26769.931 407.4130 26613.414 27801.418 27161.931
bestgrp[["Obese"]]
#> class: DMN 
#> k: 3 
#> samples x taxa: 193 x 130 
#> Laplace: 26613.41 BIC: 27801.42 AIC: 27161.93 
```
