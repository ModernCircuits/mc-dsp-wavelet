#ifndef CWT_H_
#define CWT_H_

#include "wavefunc.h"

void cwavelet(double const* y, int N, double dt, int mother, double param, double s0, double dj, int jtot, int npad,
    double* wave, double const* scale, double* period, double* coi);

void psi0(int mother, double param, double* val, int* real);

auto cdelta(int mother, double param, double psi0) -> double;

void icwavelet(double const* wave, int N, double* scale, int jtot, double dt, double dj, double cdelta, double psi0, double* oup);

#endif /* WAVELIB_H_ */
