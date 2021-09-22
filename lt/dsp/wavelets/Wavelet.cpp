#include "Wavelet.hpp"

#include "filters/utility.hpp"

Wavelet::Wavelet(char const* name)
    : name_ { name }
    , size_ { static_cast<std::size_t>(::waveletFilterLength(name)) }
    , params_ { std::make_unique<double[]>(4 * size_) }
    , lpd_ { &params_[0], size_ }
    , hpd_ { &params_[size_], size_ }
    , lpr_ { &params_[2 * size_], size_ }
    , hpr_ { &params_[3 * size_], size_ }
{
    auto* p = params_.get();
    if (name != nullptr) {
        waveletFilterCoefficients(name, p, p + size_, p + 2 * size_, p + 3 * size_);
    }
}

auto summary(Wavelet const& obj) -> void
{
    auto const n = obj.size();
    printf("\n");
    printf("Wavelet Name: %s \n", obj.name().c_str());
    printf("\n");
    printf("Wavelet Filters \n");
    printf("lpd: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.lpd()[i]);
    }
    printf("%g", obj.lpd()[n - 1]);
    printf("] \n");
    printf("hpd: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.hpd()[i]);
    }
    printf("%g", obj.hpd()[n - 1]);
    printf("] \n");
    printf("lpr: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.lpr()[i]);
    }
    printf("%g", obj.lpr()[n - 1]);
    printf("] \n");
    printf("hpr: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.hpr()[i]);
    }
    printf("%g", obj.hpr()[n - 1]);
    printf("] \n");
}