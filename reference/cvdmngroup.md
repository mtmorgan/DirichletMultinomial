# Cross-validation on Dirichlet-Multinomial classifiers.

Run cross-validation on Dirichlet-Multinomial generative classifiers.

## Usage

``` r
cvdmngroup(ncv, count, k, z, ..., verbose = FALSE,
    .lapply = parallel::mclapply)
```

## Arguments

- ncv:

  `integer(1)` number of cross-validation groups, between 2 and
  `nrow(count)`.

- count:

  `matrix` of sample x taxon counts, subsets of which are used for
  training and cross-validation.

- k:

  named [`integer()`](https://rdrr.io/r/base/integer.html) vector of
  groups and number of Dirichlet components; e.g., `c(Lean=1, Obese=3)`
  performs cross-validation for models with `k=1` Dirichlet components
  for the ‘Lean’ group, `k=3` Dirichlet components for ‘Obese’.

- z:

  True group assignment.

- ...:

  Additional arguments, passed to
  [`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md)
  during each cross-validation.

- verbose:

  `logical(1)` indicating whether progress should be reported

- .lapply:

  A function used to perform the outer cross-vaildation loop, e.g.,
  `lapply` for calculation on a single processor,
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) for
  parallel evaluation.

## Value

A `data.frame` summarizing classifications of test samples in
cross-validation groups. Columns are:

- group:

  The cross-validation group in which the indivdual was used for
  testing.

- additional columns:

  Named after classification groups, giving the posterior probability of
  assignment.

## Author

Martin Morgan
[mailto:mtmorgan.xyz@gmail.com](mailto:mtmorgan.xyz@gmail.com)

## See also

[`dmn`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmn.md),
[DirichletMultinomial-package](https://mtmorgan.github.io/DirichletMultinomial/reference/DirichletMultinomial-package.md),
[`vignette("DirichletMultinomial")`](https://mtmorgan.github.io/DirichletMultinomial/articles/DirichletMultinomial.md)

## Examples
