#ifndef CWT_H_
#define CWT_H_

#include "wavefunc.h"

#ifdef __cplusplus
extern "C" {
#endif

void cwavelet(const double* y, int N, double dt, int mother, double param, double s0, double dj, int jtot, int npad,
    double* wave, const double* scale, double* period, double* coi);

void psi0(int mother, double param, double* val, int* real);

auto factorial(int N) -> double;

auto cdelta(int mother, double param, double psi0) -> double;

void icwavelet(const double* wave, int N, double* scale, int jtot, double dt, double dj, double cdelta, double psi0, double* oup);

#ifdef __cplusplus
}
#endif

#endif /* WAVELIB_H_ */
