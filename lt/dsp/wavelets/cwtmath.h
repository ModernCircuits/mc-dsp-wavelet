#ifndef CWTMATH_H_
#define CWTMATH_H_

#include "FFT.hpp"
#include "wtmath.h"

/// \brief lb -lower bound, ub - upper bound, w - time or frequency grid (Size N)
auto nsfftExec(FftSet* obj, FftData* inp, FftData* oup, double lb, double ub, double* w) -> void;

auto cwtGamma(double x) -> double;

#endif /* WAVELIB_H_ */
