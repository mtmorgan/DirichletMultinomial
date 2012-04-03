#ifndef _DIRICHLET_FIT_MAIN_H_
#define _DIRICHLET_FIT_MAIN_H_

#include <R_ext/Boolean.h>

struct data_t {
    Rboolean verbose;
    int N, S, K;
    const int *aanX;
    double* adPi;
    /* result */
    double NLE, LogDet;
    double *group;
    double *mixture_wt;
    double fit_laplace, fit_bic, fit_aic,
        *fit_lower, *fit_mpe, *fit_upper;
};

void dirichlet_fit_main(struct data_t *data, int rseed);

#endif
