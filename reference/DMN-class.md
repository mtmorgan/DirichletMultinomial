# Class `"DMN"`

Result from fitting a Dirichlet-Multinomial model.

## Objects from the Class

Objects can be created by calls to
[`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md)..

## Slots

The contents of a slot is usually retrieved via the methods described on
the
[`mixture`](https://mtmorgan.github.io/DirichletMultinomial/reference/fitted.md)
help page.

- goodnessOfFit:

  NLE, LogDet, Laplace, AIC, and BIC criteria assessing goodness-of-fit.

- group:

  `matrix` of dimension samples x `k`, providing the Dirichlet parameter
  vectors.

- mixture:

  Weight

  :   [`numeric()`](https://rdrr.io/r/base/numeric.html) of length `k`,
      with relative weight of each component.

- fit:

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

## Methods

See the
[`mixture`](https://mtmorgan.github.io/DirichletMultinomial/reference/fitted.md)
help page.

## Author

Martin Morgan
[mailto:mtmorgan.xyz@gmail.com](mailto:mtmorgan.xyz@gmail.com)

## See also

[`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md),
[`mixture`](https://mtmorgan.github.io/DirichletMultinomial/reference/fitted.md).

## Examples

``` r
data(fit)
fit[[4]]
#> class: DMN 
#> k: 4 
#> samples x taxa: 278 x 130 
#> Laplace: 38781.1 BIC: 40425.31 AIC: 39476.69 
```
