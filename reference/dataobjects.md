# Data objects used for examples and the vignette

These data objects correspond to steps in a typical work flow, as
described in the vignette to this package. `fit` corresponds to `dmn`
fits to different values of `k`. `bestgroup` is the result of the
two-group generative classifier. `xval` summarizes leave-one-out cross
validation of the classifier.

## Usage

``` r
data(fit)
data(bestgrp)
data(xval)
```

## Format

`fit` is a list of seven
[`DMN`](https://mtmorgan.github.io/DirichletMultinomial/reference/DMN-class.md)
objects.

`bestgrp` is a
[`DMNGroup`](https://mtmorgan.github.io/DirichletMultinomial/reference/DMNGroup-class.md)
object.

`xval` is a `data.frame` with columns corresponding to the
cross-validation group membership and the Lean and Obese posterior
probabilities.

## Examples

``` r
data(fit); fit[1:2]
#> [[1]]
#> class: DMN 
#> k: 1 
#> samples x taxa: 278 x 130 
#> Laplace: 39227.36 BIC: 39527.91 AIC: 39292.11 
#> 
#> [[2]]
#> class: DMN 
#> k: 2 
#> samples x taxa: 278 x 130 
#> Laplace: 38872.68 BIC: 39588.93 AIC: 39115.53 
#> 
plot(sapply(fit, laplace), type="b")

data(bestgrp); bestgrp
#> class: DMNGroup 
#> summary:
#>       k samples taxa       NLE   LogDet   Laplace       BIC       AIC
#> Lean  1      61  130  9065.657 162.3513  9027.371  9332.864  9195.657
#> Obese 3     193  130 26769.931 407.4130 26613.414 27801.418 27161.931
data(xval); head(xval, 3)
#>         group         Lean     Obese
#> TS119.2     1 0.0001132595 0.9998867
#> TS5         2 0.0047642534 0.9952357
#> TS72        3 0.6170744622 0.3829255
```
