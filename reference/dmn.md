# Fit Dirichlet-Multinomial models to count data.

Fit Dirichlet-Multinomial models to a sample x taxon count matrix.

## Usage

``` r
dmn(count, k, verbose = FALSE, seed = runif(1, 0, .Machine$integer.max))
```

## Arguments

- count:

  [`matrix()`](https://rdrr.io/r/base/matrix.html) of sample x taxon
  counts.

- k:

  `integer(1)`, the number of Dirichlet components to fit.

- verbose:

  `logical(1)` indicating whether progress in fit should be reported.

- seed:

  `numeric(1)` random number seed.

## Details

This implements Dirichlet-multinomial mixture models describe in the
package help page,
[DirichletMultinomial-package](https://mtmorgan.github.io/DirichletMultinomial/reference/DirichletMultinomial-package.md).

## Value

An object of class `dmn`, with elements (elements are usually retrieved
via functions defined in the package, not directly).

- GoodnessOfFit:

  NLE, LogDet, Laplace, AIC, and BIC criteria assessing goodness-of-fit.

- Group:

  `matrix` of dimension samples x `k`, providing the Dirichlet parameter
  vectors.

- Mixture:

  Weight

  :   [`numeric()`](https://rdrr.io/r/base/numeric.html) of length `k`,
      with relative weight of each component.

- Fit:

  Lower

  :   [`matrix()`](https://rdrr.io/r/base/matrix.html) of dimension taxa
      x `k` with 95% lower bounds on Dirichlet component vector
      estimates.

  Estimate

  :   [`matrix()`](https://rdrr.io/r/base/matrix.html) of dimension taxa
      x `k` with Dirichlet component vector estimates.

  Upper

  :   [`matrix()`](https://rdrr.io/r/base/matrix.html) of dimension taxa
      x `k` with 95% upper bounds on Dirichlet component vector
      estimates.

## References

Holmes I, Harris K, Quince C, 2012 Dirichlet Multinomial Mixtures:
Generative Models for Microbial Metagenomics. PLoS ONE 7(2): e30126.
doi:10.1371/journal.pone.0030126.

## Author

Martin Morgan
[mailto:mtmorgan.xyz@gmail.com](mailto:mtmorgan.xyz@gmail.com)

## See also

[DirichletMultinomial-package](https://mtmorgan.github.io/DirichletMultinomial/reference/DirichletMultinomial-package.md),
[`vignette("DirichletMultinomial")`](https://mtmorgan.github.io/DirichletMultinomial/articles/DirichletMultinomial.md)

## Examples

``` r
data(fit)
## k = 1:7; full example in vignette
lplc <- sapply(fit, laplace)
plot(lplc, type="b")

fit[[which.min(lplc)]]
#> class: DMN 
#> k: 4 
#> samples x taxa: 278 x 130 
#> Laplace: 38781.1 BIC: 40425.31 AIC: 39476.69 
```
