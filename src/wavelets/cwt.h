#ifndef CWT_H_
#define CWT_H_

#include "wavefunc.h"

auto cwavelet(double const* y, int n, double dt, int mother, double param, double s0, double dj, int jtot, int npad,
    double* wave, double const* scale, double* period, double* coi) -> void;

auto psi0(int mother, double param, double* val, int* real) -> void;

auto cdelta(int mother, double param, double psi0) -> double;

auto icwavelet(double const* wave, int n, double* scale, int jtot, double dt, double dj, double cdelta, double psi0, double* oup) -> void;

#endif /* WAVELIB_H_ */
