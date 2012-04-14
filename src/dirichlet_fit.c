#include "dirichlet_fit_main.h"
#include "dirichlet_fit.h"

SEXP dirichlet_fit(SEXP counts, SEXP n_components, SEXP verbose,
                   SEXP seed)
{
    /* counts: N communities x S taxa */
    struct data_t *data =
        (struct data_t *) R_alloc(1, sizeof(struct data_t));

    /* inputs */
    SEXP dim = Rf_getAttrib(counts, R_DimSymbol),
        dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    data->verbose = LOGICAL(verbose)[0];
    data->N = INTEGER(dim)[0];
    data->S = INTEGER(dim)[1];
    data->K = INTEGER(n_components)[0];
    data->aanX = INTEGER(counts);

    /* results */
    SEXP result, elt, sxp, nms;

    PROTECT(result = Rf_allocVector(VECSXP, 4));
    nms = Rf_allocVector(STRSXP, 4);
    Rf_namesgets(result, nms);
    SET_STRING_ELT(nms, 0, mkChar("GoodnessOfFit"));
    SET_STRING_ELT(nms, 1, mkChar("Group"));
    SET_STRING_ELT(nms, 2, mkChar("Mixture"));
    SET_STRING_ELT(nms, 3, mkChar("Fit"));

    sxp = Rf_allocVector(REALSXP, 5); /* GoodnessOfFit */
    SET_VECTOR_ELT(result, 0, sxp);
    nms = Rf_allocVector(STRSXP, 5);
    Rf_namesgets(sxp, nms);
    SET_STRING_ELT(nms, 0, mkChar("NLE"));
    SET_STRING_ELT(nms, 1, mkChar("LogDet"));
    SET_STRING_ELT(nms, 2, mkChar("Laplace"));
    SET_STRING_ELT(nms, 3, mkChar("BIC"));
    SET_STRING_ELT(nms, 4, mkChar("AIC"));

    sxp = Rf_allocMatrix(REALSXP, data->N, data->K); /* Group */
    SET_VECTOR_ELT(result, 1, sxp);
    nms = Rf_allocVector(VECSXP, 2);
    Rf_dimnamesgets(sxp, nms);
    SET_VECTOR_ELT(nms, 0, VECTOR_ELT(dimnames, 0));
    SET_VECTOR_ELT(nms, 1, R_NilValue);
    data->group = REAL(sxp);

    elt = Rf_allocVector(VECSXP, 1); /* Mixture */
    SET_VECTOR_ELT(result, 2, elt);
    nms = Rf_allocVector(STRSXP, 1);
    Rf_namesgets(elt, nms);
    SET_STRING_ELT(nms, 0, mkChar("Weight"));

    sxp = Rf_allocVector(REALSXP, data->K);
    SET_VECTOR_ELT(elt, 0, sxp);
    data->mixture_wt = REAL(sxp);

    elt = Rf_allocVector(VECSXP, 3); /* Fit */
    SET_VECTOR_ELT(result, 3, elt);
    nms = Rf_allocVector(STRSXP, 3);
    Rf_namesgets(elt, nms);
    SET_STRING_ELT(nms, 0, mkChar("Lower"));
    SET_STRING_ELT(nms, 1, mkChar("Estimate"));
    SET_STRING_ELT(nms, 2, mkChar("Upper"));

    PROTECT(nms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(nms, 0, VECTOR_ELT(dimnames, 1));
    SET_VECTOR_ELT(nms, 1, R_NilValue);
    for (int i = 0; i < 3; ++i) {
        sxp = Rf_allocMatrix(REALSXP, data->S, data->K);
        SET_VECTOR_ELT(elt, i, sxp);
        Rf_dimnamesgets(sxp, nms);
    }
    UNPROTECT(1);
    data->fit_lower = REAL(VECTOR_ELT(elt, 0));
    data->fit_mpe = REAL(VECTOR_ELT(elt, 1));
    data->fit_upper = REAL(VECTOR_ELT(elt, 2));

    dirichlet_fit_main(data, INTEGER(seed)[0]);

    elt = VECTOR_ELT(result , 0);
    REAL(elt)[0] = data->NLE;
    REAL(elt)[1] = data->LogDet;
    REAL(elt)[2] = data->fit_laplace;
    REAL(elt)[3] = data->fit_bic;
    REAL(elt)[4] = data->fit_aic;

    UNPROTECT(1);
    return result;
}
