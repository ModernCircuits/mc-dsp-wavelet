#pragma once

#include "lt/cstddef.hpp"

namespace lt
{
namespace dsp
{

auto dwtPerStride(float const* inp, int n, float const* lpd, float const* hpd, int lpdLen, float* cA, int lenCA,
                  float* cD, int istride, int ostride) -> void;

auto dwtSymStride(float const* inp, int n, float const* lpd, float const* hpd, int lpdLen, float* cA, int lenCA,
                  float* cD, int istride, int ostride) -> void;

auto modwtPerStride(int m, float const* inp, int n, float const* filt, int lpdLen, float* cA, int lenCA, float* cD,
                    int istride, int ostride) -> void;

auto swtPerStride(int m, float const* inp, int n, float const* lpd, float const* hpd, int lpdLen, float* cA, int lenCA,
                  float* cD, int istride, int ostride) -> void;

auto idwtPerStride(float const* cA, int lenCA, float const* cD, float const* lpr, float const* hpr, int lprLen,
                   float* x, int istride, int ostride) -> void;

auto idwtSymStride(float const* cA, int lenCA, float const* cD, float const* lpr, float const* hpr, int lprLen,
                   float* x, int istride, int ostride) -> void;

auto testSWTlength(int n, int j) -> int;

auto maxIterations(std::size_t sigLen, std::size_t filtLen) -> std::size_t;

}  // namespace dsp
}  // namespace lt