#include <R.h>
#include <Rinternals.h>
#include <complex.h>


/*
  Discrete fourier transform, moving heavy stuff of fourier.transform into C code

  Inputs:
    - z: an (C-) Array of (implicit) dimension dim_z1 times dim_z2 times n_lags
    - freq: a vector of real-valued frequencies (values on [-pi, pi])
    - lags: an integer vector of lags

  Returns:
    A numeric vector of length dim_z1 * dim_z2 * n_freq



*/
SEXP fourier_transform(SEXP z, SEXP dim_z1, SEXP dim_z2,
                       SEXP freq, SEXP n_freq,
                       SEXP lags, SEXP n_lags)
{

    int res_array_ln = *INTEGER(dim_z1) * *INTEGER(dim_z2) * *INTEGER(n_freq);
    int sub_dim = *INTEGER(dim_z1) * *INTEGER(dim_z2);

    // initialize res to zero
    SEXP res = PROTECT(allocVector(CPLXSXP, res_array_ln));
    Rcomplex cpl;
    cpl.r = 0.0; cpl.i = 0.0;
    for (int i = 0; i < res_array_ln; i++) {
        COMPLEX(res)[i] = cpl;
    }

    double complex w;

    for (int i = 0; i< *INTEGER(n_freq); i++) {
        for (int j = 0; j < *INTEGER(n_lags); j++) {
            w = cexp(0.0 - INTEGER(lags)[j] * REAL(freq)[i] * I);

            for (int k = 0; k < sub_dim; k++) {
                (COMPLEX(res)[i * sub_dim + k]).r += creal(w) * REAL(z)[j * sub_dim + k];
                (COMPLEX(res)[i * sub_dim + k]).i += cimag(w) * REAL(z)[j * sub_dim + k];
            }
        }
    }


    /* double complex accum = 0 + 0 * I; */
    /* int res_idx = 0; */
    /* while (res_idx < res_array_ln) { */
    /*     temp_z.r = creal(accum); temp_z.i = cimag(accum); */
    /*     COMPLEX(res)[res_idx] = temp_z; */
    /*     res_idx +=1; */
    /* } */

    UNPROTECT(1);
    return res;
}


static const R_CallMethodDef CallEntries[] = {
        {"fourier_transform", (DL_FUNC) &fourier_transform, 7},
        {NULL, NULL, 0}
};



void R_init_freqdom(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
