#include <R_ext/Rdynload.h>
#include "dirichlet_fit.h"

static const R_CallMethodDef callMethods[] = {
    { ".dirichlet_fit", (DL_FUNC) &dirichlet_fit, 4},
    {NULL, NULL, 0}
};

void R_init_DirichletMultinomial(DllInfo * info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

