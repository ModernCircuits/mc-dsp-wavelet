#pragma once

auto meyer(int n, double lb, double ub, double* phi, double* psi, double* tgrid) -> void;
auto gauss(int n, int p, double lb, double ub, double* psi, double* t) -> void;
auto mexhat(int n, double lb, double ub, double* psi, double* t) -> void;
auto morlet(int n, double lb, double ub, double* psi, double* t) -> void;

auto cwtGamma(double x) -> double;