# Access model components.

The accessors `mixture` and `mixturewt` return information about the
estimated Dirichlet components of the fitted model. Return values are
described in the Values section, below.

## Usage

``` r
mixture(object, ..., assign=FALSE)
mixturewt(object, ...)
goodnessOfFit(object, ...)
laplace(object, ...)
# S4 method for class 'DMN'
AIC(object, ..., k = 2)
# S4 method for class 'DMN'
BIC(object, ...)

# S4 method for class 'DMN'
fitted(object, ..., scale=FALSE)
# S4 method for class 'DMN'
predict(object, newdata, ..., logevidence=FALSE)
# S4 method for class 'DMNGroup'
fitted(object, ...)
# S4 method for class 'DMNGroup'
predict(object, newdata, ..., assign=FALSE)
# S4 method for class 'DMNGroup'
summary(object, ...)
```

## Arguments

- object:

  An instance of class
  [`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md).

- newdata:

  A `matrix` of new sample x taxon data to be fitted to the model of
  `object`.

- ...:

  Additional arguments, available to methods, when applicable.

- assign:

  `logical(1)` indicating whether the maximum per-sample mixture
  component should be returned (`assign=FALSE`), or the full mixture
  matrix (`assign=TRUE`).

- scale:

  `logical(1)` indicating whether fitted values should be returned
  unscaled (default, `scaled=FALSE`) or scaled by the variability of
  `mixturewt` parameter `theta`.

- logevidence:

  `logical(1)` indicating whether posterior probability (default,
  `logevidence=FALSE`) or log evidence `logical=TRUE` should be
  returned.

- k:

  ignored.

## Value

`mixture` with `assign=FALSE` returns a matrix of sample x Dirichlet
component estimates. With `assign=TRUE` `mixture` returns a named vector
indexing the maximal Dirichlet component of each sample.

`mixturewt` returns a matrix with rows corresponding to mixture
components, and columns `pi` (component weight) and `theta` (component
variability). Small values of `theta` correspond to highly variable
components.

`goodnessOfFit` returns a named numeric vector of measures of goodness
of fit.

`laplace`, `AIC`, and `BIC` return the corresponding measures of
goodness of fit.

## Author

Martin Morgan
[mailto:mtmorgan.xyz@gmail.com](mailto:mtmorgan.xyz@gmail.com)

## Examples

``` r
data(fit)
best <- fit[[4]]
mixturewt(best)
#>          pi    theta
#> 1 0.3108456 52.03706
#> 2 0.1665874 18.72599
#> 3 0.3027727 53.29525
#> 4 0.2197943 30.19582
head(mixture(best), 3)
#>                 [,1]         [,2]         [,3]         [,4]
#> TS1.2   9.999914e-01 2.117284e-11 8.563935e-06 3.306464e-08
#> TS10.2  3.776510e-08 3.268129e-04 9.996731e-01 2.847131e-10
#> TS100.2 7.214444e-09 8.825346e-01 7.953749e-13 1.174654e-01
head(mixture(best, assign=TRUE), 3)
#>   TS1.2  TS10.2 TS100.2 
#>       1       3       2 
goodnessOfFit(best)
#>        NLE     LogDet    Laplace        BIC        AIC 
#> 38953.6920   616.0335 38781.1039 40425.3149 39476.6920 

fl <- system.file(package="DirichletMultinomial", "extdata",
                  "Twins.csv")
count <- t(as.matrix(read.csv(fl, row.names=1)))
data(bestgrp)
bestgrp
#> class: DMNGroup 
#> summary:
#>       k samples taxa       NLE   LogDet   Laplace       BIC       AIC
#> Lean  1      61  130  9065.657 162.3513  9027.371  9332.864  9195.657
#> Obese 3     193  130 26769.931 407.4130 26613.414 27801.418 27161.931
head(predict(bestgrp, count))
#>                 Lean      Obese
#> TS1.2   9.648780e-01 0.03512197
#> TS10.2  1.000058e-03 0.99899994
#> TS100.2 3.522984e-08 0.99999996
#> TS100   3.290371e-05 0.99996710
#> TS101.2 7.349397e-08 0.99999993
#> TS103.2 1.679035e-02 0.98320965
```
