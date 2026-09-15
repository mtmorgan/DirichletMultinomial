# Helpful utility functions

`csubset` creates a subset of a count matrix, based on identity of
column phenotypes to a specified value.

## Usage

``` r
csubset(val, x, pheno, cidx = TRUE)
```

## Arguments

- val:

  `character(1)` specifying the subset of `phenotype` to select.

- x:

  A matrix of counts, with rows corresponding to samples and columns to
  taxonomic groups.

- pheno:

  A [`character()`](https://rdrr.io/r/base/character.html) vector of
  length equal to the number of rows in `count`, indicating the
  phenotype of the corresponding sample.

- cidx:

  A `logical(1)` indicating whether columns (taxa) with zero counts in
  the count matrix following removal of taxa not satisfying
  `pheno %in% val` should be removed. `cidx=FALSE` removes the 0-count
  columns.

## Value

A `matrix` of counts, with rows satisfying `pheno %in% val` and with
columns equal either to `ncol(x)` (when `cidx=TRUE`) or the number of
columns with non-zero counts after row subsetting (`cidx=FALSE`).

## Author

Martin Morgan
[mailto:mtmorgan.xyz@gmail.com](mailto:mtmorgan.xyz@gmail.com)

## Examples

``` r
## count matrix
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "Twins.csv")
count <- t(as.matrix(read.csv(fl, row.names=1)))

## phenotype
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "TwinStudy.t")
pheno0 <- scan(fl)
lvls <- c("Lean", "Obese", "Overwt")
pheno <- factor(lvls[pheno0 + 1], levels=lvls)
names(pheno) <- rownames(count)

## subset
dim(count)
#> [1] 278 130
sum("Lean" == pheno)
#> [1] 61
dim(csubset("Lean", count, pheno))
#> [1]  61 130
dim(csubset("Lean", count, pheno, cidx=FALSE))
#> [1]  61 106
```
