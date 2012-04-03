#ifndef _DIRICHLET_FIT_H_
#define _DIRICHLET_FIT_H_

#include <Rdefines.h>

SEXP dirichlet_fit(SEXP counts, SEXP n_components, SEXP verbose,
                   SEXP seed);

#endif
