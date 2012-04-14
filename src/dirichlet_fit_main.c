#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>

#include <Rdefines.h>
#include "dirichlet_fit_main.h"
/* re-map to R transient memory allocation */
#define calloc(_nelm, _elsize) R_alloc(_nelm, _elsize)
#define free(_ptr) (void) _ptr

const double BIG_DBL = 1.0e9;
const double K_MEANS_THRESH = 1.0e-6;
const double SOFT_BETA = 50.0;
const double GAMMA_ITA = 0.1;
const double GAMMA_NU = 0.1;
const int MAX_ITER = 1000;
const size_t MAX_GRAD_ITER = 1000;

static void kmeans(struct data_t *data, gsl_rng *ptGSLRNG,
                   double* adW, double **aadZ, double **aadMu)
{
    const int S = data->S, N = data->N, K = data->K,
        *aanX = data->aanX;
    int i, j, k, iter = 0;

    double **aadY, *adMu;
    double dMaxChange = BIG_DBL;

    if (data->verbose)
        Rprintf("  Soft kmeans\n");

    aadY = (double **) calloc(N, sizeof(double *));
    aadY[0] = (double *) calloc(N * S, sizeof(double));

    adMu = (double *) calloc(S, sizeof(double));

    for (i = 0; i < N; i++) {
        double dTotal = 0.0;
        aadY[i] = aadY[0] + i * S;
        for (j = 0; j < S; j++)
            dTotal += aanX[j * N + i];
        for (j = 0; j < S; j++)
            aadY[i][j] = (aanX[j * N + i]) / dTotal;
    }

    /* initialse */
    for (i = 0; i < N; i++) {
        k = gsl_rng_uniform_int (ptGSLRNG, K);
        for (j = 0; j < K; j++)
            aadZ[j][i] = 0.0;
        aadZ[k][i] = 1.0;
    }

    while (dMaxChange > K_MEANS_THRESH && iter < MAX_ITER) {
        /* update mu */
        dMaxChange = 0.0;
        for (i = 0; i < K; i++){
            double dNormChange = 0.0;
            adW[i] = 0.0;
            for (j = 0; j < N; j++)
                adW[i] += aadZ[i][j];
            for (j = 0; j < S; j++) {
                adMu[j] = 0.0;
                for (k = 0; k < N; k++)
                    adMu[j] += aadZ[i][k]*aadY[k][j];
            }

            for (j = 0; j < S; j++) {
                double dDiff = 0.0;
                adMu[j] /= adW[i];
                dDiff = (adMu[j] - aadMu[i][j]);
                dNormChange += dDiff * dDiff;
                aadMu[i][j] = adMu[j];
            }
            dNormChange = sqrt(dNormChange);
            if (dNormChange > dMaxChange)
                dMaxChange = dNormChange;
        }

        /* calc distances and update Z */
        for (i = 0; i < N; i++) {
            double dNorm = 0.0, adDist[K];
            for (k = 0; k < K; k++) {
                adDist[k] = 0.0;
                for (j = 0; j < S; j++) {
                    double dDiff = (aadMu[k][j] - aadY[i][j]);
                    adDist[k] += dDiff * dDiff;
                }
                adDist[k] = sqrt(adDist[k]);
                dNorm += exp(-SOFT_BETA * adDist[k]);
            }
            for (k = 0; k < K; k++)
                aadZ[k][i] = exp(-SOFT_BETA * adDist[k]) / dNorm;
        }
        iter++;
        if (data->verbose && (iter % 50 == 0))
            Rprintf("    iteration %d change %f\n", iter, dMaxChange);
    }

    free(aadY[0]); free(aadY);
    free(adMu);
}

static double neg_log_evidence_lambda_pi(const gsl_vector *lambda,
                                         void *params)
{
    int i, j;
    const struct data_t *data = (const struct data_t *) params;
    const int S = data->S, N = data->N, *aanX = data->aanX;
    const double *adPi = data->adPi;
    double dLogE = 0.0, dLogEAlpha = 0.0, dSumAlpha = 0.0, dSumLambda = 0.0;
    double adSumAlphaN[N], dWeight = 0.0;

    for (i = 0; i < N; i++) {
        adSumAlphaN[i] = 0.0;
        dWeight += adPi[i];
    }

    for (j = 0; j < S; j++) {
        const double dLambda = gsl_vector_get(lambda, j);
        const double dAlpha = exp(dLambda);
        dLogEAlpha += gsl_sf_lngamma(dAlpha);
        dSumLambda += dLambda;
        dSumAlpha += dAlpha;
        for (i = 0; i < N; i++) {
            const double dN = aanX[j * N + i];
            const double dAlphaN = dAlpha + dN;
            adSumAlphaN[i] += dAlphaN; /*weight by pi*/
            dLogE -= adPi[i] * gsl_sf_lngamma(dAlphaN); /*weight by pi*/
        }
    }
    dLogEAlpha -= gsl_sf_lngamma(dSumAlpha);

    for(i = 0; i < N; i++)
        dLogE += adPi[i] * gsl_sf_lngamma(adSumAlphaN[i]);

    return dLogE + dWeight*dLogEAlpha + GAMMA_NU*dSumAlpha -
        GAMMA_ITA * dSumLambda;
}

static void neg_log_derive_evidence_lambda_pi(const gsl_vector *ptLambda,
                                              void *params, gsl_vector* g)
{
    const struct data_t *data = (const struct data_t *) params;
    const int S = data->S, N = data->N, *aanX = data->aanX;
    const double *adPi = data->adPi;

    int i, j;
    double adDeriv[S], adStore[N], adAlpha[S];
    double dSumStore = 0.0, dStore = 0.0;
    double dWeight = 0;

    for (i = 0; i < N; i++) {
        adStore[i] = 0.0;
        dWeight += adPi[i];
    }

    for (j = 0; j < S; j++) {
        adAlpha[j] = exp(gsl_vector_get(ptLambda, j));
        dStore += adAlpha[j];
        adDeriv[j] = dWeight* gsl_sf_psi(adAlpha[j]);
        for (i = 0; i < N; i++) {
            double dAlphaN = adAlpha[j] + aanX[j * N + i];
            adDeriv[j] -= adPi[i]*gsl_sf_psi (dAlphaN);
            adStore[i] += dAlphaN;
        }
    }

    for (i = 0; i < N; i++)
        dSumStore += adPi[i] * gsl_sf_psi(adStore[i]);
    dStore = dWeight * gsl_sf_psi(dStore);

    for (j = 0; j < S; j++) {
        double value = adAlpha[j] *
            (GAMMA_NU + adDeriv[j] - dStore + dSumStore) - GAMMA_ITA;
        gsl_vector_set(g, j, value);
    }
}

static void neg_log_FDF_lamba_pi(const gsl_vector *x, void *params,
                                 double *f, gsl_vector *g)
{
    *f = neg_log_evidence_lambda_pi(x, params);
    neg_log_derive_evidence_lambda_pi(x, params, g);
}

static void optimise_lambda_k(double *adLambdaK, struct data_t *data,
                              double *adZ)
{
    const int S = data->S;

    int i, status;
    size_t iter = 0;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    gsl_multimin_function_fdf fdf;
    gsl_vector *ptLambda;

    /*initialise vector*/
    ptLambda = gsl_vector_alloc(S);
    for (i = 0; i < S; i++)
        gsl_vector_set(ptLambda, i, adLambdaK[i]);

    /*initialise function to be solved*/
    data->adPi = adZ;
    fdf.n = S;
    fdf.f = neg_log_evidence_lambda_pi;
    fdf.df = neg_log_derive_evidence_lambda_pi;
    fdf.fdf = neg_log_FDF_lamba_pi;
    fdf.params = data;

    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, S);

    gsl_multimin_fdfminimizer_set(s, &fdf, ptLambda, 1.0e-6, 0.1);

    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);
        if (status)
            break;
        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
   } while (status == GSL_CONTINUE && iter < MAX_GRAD_ITER);

    for (i = 0; i < S; i++)
        adLambdaK[i] = gsl_vector_get(s->x, i);

    gsl_vector_free(ptLambda);
    gsl_multimin_fdfminimizer_free(s);
}

static double neg_log_evidence_i(const struct data_t *data,
                                 const int *anX, const double* adLambda)
{
    int j;
    const int S = data->S, N = data->N;
    double dLogE = 0.0, dLogEAlpha = 0.0, dSumAlpha = 0.0,
        dSumAlphaN = 0.0;

    for (j = 0; j < S; j++) {
        const double dAlpha = exp(adLambda[j]);
        const double dAlphaN = anX[j * N] + dAlpha;

        dLogEAlpha += gsl_sf_lngamma(dAlpha);
        dSumAlpha += dAlpha;
        dSumAlphaN += dAlphaN;
        dLogE -= gsl_sf_lngamma(dAlphaN);
    }

    dLogEAlpha -= gsl_sf_lngamma(dSumAlpha);
    dLogE += gsl_sf_lngamma(dSumAlphaN);

    return dLogE + dLogEAlpha;
}

static void calc_z(double **aadZ, const struct data_t *data,
                   const double *adW, double **aadLambda)
{
    int i, k;
    const int N = data->N, K = data->K;
    double adStore[K];

    for (i = 0; i < N; i ++) {
        double dSum = 0.0;
        double dOffset = BIG_DBL;
        for (k = 0; k < K; k++) {
            double dNegLogEviI =
                neg_log_evidence_i(data, data->aanX + i, aadLambda[k]);
            if (dNegLogEviI < dOffset)
                dOffset = dNegLogEviI;
            adStore[k] = dNegLogEviI;
        }
        for (k = 0; k < K; k++) {
            aadZ[k][i] = adW[k] * exp(-(adStore[k] - dOffset));
            dSum += aadZ[k][i];
        }
        for (k = 0; k < K; k++)
            aadZ[k][i] /= dSum;
    }
}

static double neg_log_likelihood(double *adW, double** aadLambda,
                                 const struct data_t *data)
{
    const int S = data->S, N = data->N, K = data->K,
        *aanX = data->aanX;
    int i, j, k;
    double adPi[K], adLogBAlpha[K];
    double dRet = 0.0, dL5 = 0.0, dL6 = 0.0, dL7 = 0.0, dL8 = 0.0;
    double dK = K, dN = N, dS = S;

    for (k = 0; k < K; k++){
        double dSumAlphaK = 0.0;
        adLogBAlpha[k] = 0.0;
        adPi[k] = adW[k]/dN;
        for (j = 0; j < S; j++){
            double dAlpha = exp(aadLambda[k][j]);
            dSumAlphaK += dAlpha;
            adLogBAlpha[k] += gsl_sf_lngamma(dAlpha);
        }
        adLogBAlpha[k] -= gsl_sf_lngamma(dSumAlphaK);
    }

    for (i = 0; i < N; i++) {
        double dProb = 0.0, dFactor = 0.0, dSum = 0.0, adLogStore[K],
            dOffset = -BIG_DBL;

        for (j = 0; j < S; j++) {
            dSum += aanX[j * N + i];
            dFactor += gsl_sf_lngamma(aanX[j * N + i] + 1.0);
        }
        dFactor -= gsl_sf_lngamma(dSum + 1.0);

        for (k = 0; k < K; k++) {
            double dSumAlphaKN = 0.0, dLogBAlphaN = 0.0;
            for (j = 0; j < S; j++) {
                double dAlphaN = exp(aadLambda[k][j]) + aanX[j * N + i];
                dSumAlphaKN += dAlphaN;
                dLogBAlphaN += gsl_sf_lngamma(dAlphaN);
            }
            dLogBAlphaN -= gsl_sf_lngamma(dSumAlphaKN);
            adLogStore[k] = dLogBAlphaN - adLogBAlpha[k] - dFactor;
            if (adLogStore[k] > dOffset)
                dOffset = adLogStore[k];
        }

        for (k = 0; k < K; k++)
            dProb += adPi[k]*exp(-dOffset + adLogStore[k]);
        dRet += log(dProb)+dOffset;
    }
    dL5 = -dS * dK * gsl_sf_lngamma(GAMMA_ITA);
    dL6 = GAMMA_ITA * dK * dS * log(GAMMA_NU);

    for (i = 0; i < K; i++)
        for (j = 0; j < S; j++) {
            dL7 += exp(aadLambda[i][j]);
            dL8 += aadLambda[i][j];
        }
    dL7 *= -GAMMA_NU;
    dL8 *= GAMMA_ITA;

    return -dRet -dL5 - dL6 -dL7 -dL8;
}

static void hessian(gsl_matrix* ptHessian, const double* adLambda,
                    const struct data_t *data)
{
    const int S = data->S, N = data->N, *aanX = data->aanX;
    const double *adPi = data->adPi;

    int i = 0, j = 0;
    double adAlpha[S], adAJK[S], adCJK[S], adAJK0[S], adCJK0[S];
    double dCK0 = 0.0, dAK0;
    double dCSum, dAlphaSum = 0.0, dW = 0.0, dCK = 0.0, dAK;

    for (j = 0; j < S; j++) {
        adAlpha[j] = exp(adLambda[j]);
        dAlphaSum += adAlpha[j];
        adAJK0[j] = adAJK[j] = adCJK0[j] = adCJK[j] = 0.0;
        for (i = 0; i < N; i++) {
            adCJK0[j] += adPi[i] * gsl_sf_psi(adAlpha[j] + aanX[j * N + i]);
            adAJK0[j] += adPi[i] * gsl_sf_psi(adAlpha[j]);
            adCJK[j] += adPi[i] * gsl_sf_psi_1(adAlpha[j] + aanX[j * N + i]);
            adAJK[j] += adPi[i] * gsl_sf_psi_1(adAlpha[j]);
        }
    }

    for (i = 0; i < N; i++) {
        dW += adPi[i];
        dCSum = 0.0;
        for (j = 0; j < S; j++)
            dCSum += adAlpha[j] + aanX[j * N + i];
        dCK  += adPi[i]*gsl_sf_psi_1(dCSum);
        dCK0 += adPi[i]*gsl_sf_psi(dCSum);
    }

    dAK = dW * gsl_sf_psi_1(dAlphaSum);
    dAK0 = dW * gsl_sf_psi(dAlphaSum);
    for (i = 0; i < S; i++)
        for (j = 0; j < S; j++) {
            double dVal = 0.0;
            if (i == j) {
                double dG1 = -adAlpha[i] *
                    (dAK0 - dCK0 + adCJK0[i] - adAJK0[i]);
                double dG2 = -adAlpha[i] *
                    adAlpha[i]*(dAK - dCK + adCJK[i] - adAJK[i]);
                double dG3 = adAlpha[i]*GAMMA_NU;
                dVal = dG1 + dG2 + dG3;
            } else
                dVal = -adAlpha[i] * adAlpha[j] * (dAK - dCK);
            gsl_matrix_set(ptHessian, i, j, dVal);
        }
}

static void group_output(struct data_t *data, double** aadZ)
{
    const int N = data->N, K = data->K;
    int i, k;
    for(k = 0; k < K; k++)
        for (i = 0; i < N; i++)
            data->group[k * N + i] = aadZ[k][i];
}

static void mixture_output(struct data_t *data, double *adW,
                           double** aadLambda, double **aadErr)
{
    const int N = data->N, S = data->S;
    int i, k;
    for (k = 0; k < data->K; k++)
        data->mixture_wt[k] = adW[k] / N;
    for (i = 0; i < data->S; i++) {
        for (k = 0; k < data->K; k++) {
            double dErr = aadErr[k][i], dL = 0.0, dU = 0.0;
            int bIll = FALSE;
            if (dErr >= 0.0) {
                dErr = sqrt(dErr);
                if (dErr < 100.0) {
                    dL =  exp(aadLambda[k][i] - 2.0*dErr);
                    dU =  exp(aadLambda[k][i] + 2.0*dErr);
                } else bIll = TRUE;
            } else bIll = TRUE;

            if (bIll)
                dL = dU = R_NaN;
            data->fit_lower[k * S + i] = dL;
            data->fit_mpe[k * S + i] = exp(aadLambda[k][i]);
            data->fit_upper[k * S + i] = dU;
        }
    }
}

void dirichlet_fit_main(struct data_t *data, int rseed)
{
    const int N = data->N, S = data->S, K = data->K;
    int i, j, k;

    gsl_rng *ptGSLRNG;
    gsl_rng_env_setup();
    gsl_set_error_handler_off();
    ptGSLRNG = gsl_rng_alloc(gsl_rng_default);
    gsl_set_error_handler_off();
    gsl_rng_set(ptGSLRNG, rseed);

    /* allocate matrices */
    double **aadZ, **aadLambda, **aadErr, *adW;
    adW = (double *) calloc(K, sizeof(double));

    aadZ = (double **) calloc(K, sizeof(double *));
    aadLambda = (double **) calloc(K, sizeof(double *));
    aadErr = (double **) calloc(K, sizeof(double*));

    aadZ[0] = (double *) calloc(K * N, sizeof(double));
    aadLambda[0] = (double *) calloc(K * S, sizeof(double));
    aadErr[0] = (double *) calloc(K * S, sizeof(double));

    for (k = 1; k < K; k++) {
        aadZ[k] = aadZ[0] + k * N;
        aadLambda[k] = aadLambda[0] + k * S;
        aadErr[k] = aadErr[0] + k * S;
    }

    /* soft k means initialiser */
    kmeans(data, ptGSLRNG, adW, aadZ, aadLambda);
    for (k = 0; k < K; k++) {
        adW[k] = 0.0;
        for (i = 0; i < N; i++)
            adW[k] += aadZ[k][i];
    }

    if (data->verbose)
        Rprintf("  Expectation Maximization setup\n");
    for (k = 0; k < K; k++) {
        for (j = 0; j < S; j++) {
            const double x = aadLambda[k][j];
            aadLambda[k][j] = (x > 0.0) ? log(x) : -10;
        }
        optimise_lambda_k(aadLambda[k], data, aadZ[k]);
    }

    /* simple EM algorithm */
    int iter = 0;
    double dNLL = 0.0, dNew, dChange = BIG_DBL;

    if (data->verbose)
        Rprintf("  Expectation Maximization\n");
    while (dChange > 1.0e-6 && iter < 100) {
        calc_z(aadZ, data, adW, aadLambda); /* latent var expectation */
        for (k = 0; k < K; k++) /* mixture components, given pi */
            optimise_lambda_k(aadLambda[k], data, aadZ[k]);
        for (k = 0; k < K; k++) { /* current likelihood & weights */
            adW[k] = 0.0;
            for(i = 0; i < N; i++)
                adW[k] += aadZ[k][i];
        }

        dNew = neg_log_likelihood(adW, aadLambda, data);
        dChange = fabs(dNLL - dNew);
        dNLL = dNew;
        iter++;
        R_CheckUserInterrupt();
        if (data->verbose && (iter % 10) == 0)
            Rprintf("    iteration %d change %f\n", iter, dChange);
    }

    /* hessian */
    if (data->verbose)
        Rprintf("  Hessian\n");
    gsl_matrix *ptHessian = gsl_matrix_alloc(S, S),
        *ptInverseHessian = gsl_matrix_alloc(S, S);
    gsl_permutation *p = gsl_permutation_alloc(S);
    double dLogDet = 0., dTemp;
    int signum, status;

    for (k = 0; k < K; k++) {
        data->adPi = aadZ[k];
        if (k > 0)
            dLogDet += 2.0 * log(N) - log(adW[k]);
        hessian(ptHessian, aadLambda[k], data);

        status = gsl_linalg_LU_decomp(ptHessian, p, &signum);
        gsl_linalg_LU_invert(ptHessian, p, ptInverseHessian);
        for (j = 0; j < S; j++) {
            aadErr[k][j] = gsl_matrix_get(ptInverseHessian, j, j);
            dTemp = gsl_matrix_get(ptHessian, j, j);
            dLogDet += log(fabs(dTemp));
        }
    }

    gsl_matrix_free(ptHessian);
    gsl_matrix_free(ptInverseHessian);
    gsl_permutation_free(p);

    /* results */
    double dP = K * S + K - 1;
    data->NLE = dNLL; data->LogDet = dLogDet;
    data->fit_laplace = dNLL + 0.5 * dLogDet - 0.5 * dP * log(2. * M_PI);
    data->fit_bic = dNLL + 0.5 * log(N) * dP;
    data->fit_aic = dNLL + dP;

    group_output(data, aadZ);
    mixture_output(data, adW, aadLambda, aadErr);

    free(aadErr[0]); free(aadErr);
    free(aadLambda[0]); free(aadLambda);
    free(aadZ[0]); free(aadZ);
    free(adW);
}
