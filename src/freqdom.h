#ifndef FREQDOM_H_

#include <R.h>
#include <Rinternals.h>
#include <complex.h>

SEXP fourier_inverse(SEXP f, SEXP dim_f1, SEXP dim_f2,
                     SEXP lags, SEXP n_lags,
                     SEXP freqs, SEXP n_freqs);

SEXP fourier_transform(SEXP z, SEXP dim_z1, SEXP dim_z2,
                       SEXP freq, SEXP n_freq,
                       SEXP lags, SEXP n_lags);
#define FREQDOM_H_




#endif // FREQDOM_H_
