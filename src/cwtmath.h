#ifndef CWTMATH_H_
#define CWTMATH_H_

#include "hsfft.h"
#include "wtmath.h"

/// \brief lb -lower bound, ub - upper bound, w - time or frequency grid (Size N)
void nsfft_exec(fft_object obj, fft_data* inp, fft_data* oup, double lb, double ub, double* w);

auto cwt_gamma(double x) -> double;

auto nint(double N) -> int;

#endif /* WAVELIB_H_ */
