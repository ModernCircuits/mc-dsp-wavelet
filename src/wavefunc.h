#ifndef WAVEFUNC_H_
#define WAVEFUNC_H_

#include "cwtmath.h"

void meyer(int n, double lb, double ub, double* phi, double* psi, double* tgrid);
void gauss(int n, int p, double lb, double ub, double* psi, double* t);
void mexhat(int n, double lb, double ub, double* psi, double* t);
void morlet(int n, double lb, double ub, double* psi, double* t);

#endif /* WAVEFUNC_H_ */
