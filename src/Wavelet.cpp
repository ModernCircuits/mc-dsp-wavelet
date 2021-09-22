#include "Wavelet.hpp"

#include "wavefilt.h"

Wavelet::Wavelet(char const* name)
    : name_ { name }
    , size_ { static_cast<std::size_t>(::filtlength(name)) }
    , params_ { std::make_unique<double[]>(4 * size_) }
    , lpd_ { &params_[0], size_ }
    , hpd_ { &params_[size_], size_ }
    , lpr_ { &params_[2 * size_], size_ }
    , hpr_ { &params_[3 * size_], size_ }
{
    auto* p = params_.get();
    if (name != nullptr) {
        filtcoef(name, p, p + size_, p + 2 * size_, p + 3 * size_);
    }
}