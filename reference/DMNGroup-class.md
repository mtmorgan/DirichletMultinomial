# Class `"DMNGroup"`

Result from fitting a Dirichlet-Multinomial generative classifier.

## Objects from the Class

Objects can be created by calls to
[`dmngroup`](https://mtmorgan.github.io/DirichletMultinomial/reference/dmngroup.md).

## Slots

All slots in this class are inheritted from `SimpleList`; see ‘Methods’,
below, for information on how to manipulate this object.

## Extends

Class `"SimpleList"`, directly. Class `"List"`, by class "SimpleList",
distance 2. Class `"Vector"`, by class "SimpleList", distance 3. Class
`"Annotated"`, by class "SimpleList", distance 4.

## Methods

See the
[`mixture`](https://mtmorgan.github.io/DirichletMultinomial/reference/fitted.md)
help page for functions that operate on `DMNGroup` and `DMN`.

`DMNGroup` can be manipulated as a list; see `SimpleList ` for a
description of typical list-like functions.

## Author

Martin Morgan
[mailto:mtmorgan.xyz@gmail.com](mailto:mtmorgan.xyz@gmail.com)

## See also

[`mixture`](https://mtmorgan.github.io/DirichletMultinomial/reference/fitted.md),
[`DMN`](https://mtmorgan.github.io/DirichletMultinomial/reference/DMN-class.md),
`SimpleList`.

## Examples

``` r
data(bestgrp)
bestgrp
#> class: DMNGroup 
#> summary:
#>       k samples taxa       NLE   LogDet   Laplace       BIC       AIC
#> Lean  1      61  130  9065.657 162.3513  9027.371  9332.864  9195.657
#> Obese 3     193  130 26769.931 407.4130 26613.414 27801.418 27161.931
bestgrp[[1]]
#> class: DMN 
#> k: 1 
#> samples x taxa: 61 x 130 
#> Laplace: 9027.371 BIC: 9332.864 AIC: 9195.657 
```
