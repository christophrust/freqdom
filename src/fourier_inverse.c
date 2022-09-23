#include <R.h>
#include <Rinternals.h>
#include <complex.h>






SEXP fourier_inverse(SEXP f, SEXP dim_f1, SEXP dim_f2,
                     SEXP lags, SEXP n_lags,
                     SEXP freqs, SEXP n_freqs) {

    int res_array_len = *INTEGER(dim_f1) * *INTEGER(dim_f2) * *INTEGER(n_lags);
    int sub_dim = *INTEGER(dim_f1) * *INTEGER(dim_f2);

    // initialize res object
    SEXP res = PROTECT(allocVector(REALSXP, res_array_len));
    double complex *tmp_cmplx_array;
    tmp_cmplx_array = (double complex *) Calloc(res_array_len, double complex);


    // temporary objects
    double complex z;

    for (int i = 0; i < *INTEGER(n_lags); i++){
        for (int j = 0; j < *INTEGER(n_freqs); j++) {
            z = cexp(REAL(freqs)[j] * (double)INTEGER(lags)[i] * I);
            for (int k = 0; k < sub_dim; k++) {
                tmp_cmplx_array[i * sub_dim + k] +=
                    z * (COMPLEX(f)[j * sub_dim + k].r + COMPLEX(f)[j * sub_dim + k].i * I);
            }
        }
    }

    //check that all values of tmp_cmplx_array are real and copy into res
    double accum = 0.0;
    for (int i = 0; i < res_array_len; i++){
        accum += fabs(cimag(tmp_cmplx_array[i]));
        REAL(res)[i] = creal(tmp_cmplx_array[i]);
    }

    if (accum > 1e-9) {
        warning("The imaginary part of the coefficients was not zero, probably due to an assymmetric spectrum!");
    }

    UNPROTECT(1);
    return res;
}
