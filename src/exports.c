#include "freqdom.h"


static const R_CallMethodDef CallEntries[] = {
    {"fourier_transform", (DL_FUNC) &fourier_transform, 7},
    {"fourier_inverse", (DL_FUNC) &fourier_transform, 7},
    {NULL, NULL, 0}
};



void R_init_freqdom(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
