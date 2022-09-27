# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.15)
project(mc-wavelet)

find_package(mc-dsp REQUIRED)
find_package(FFTW3f REQUIRED)
find_package(pffft REQUIRED)

add_library(mc-wavelet)
add_library(mc::wavelet ALIAS mc-wavelet)
target_include_directories(mc-wavelet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mc-wavelet
    PUBLIC
        mc-dsp::mc-dsp

        FFTW3::fftw3f
        pffft::pffft

        mc::CompilerOptions
        mc::CompilerWarnings
        mc::CodeCoverage
)

target_sources(mc-wavelet
    PRIVATE
        "mc/audio/audio_file.hpp"

        "mc/dsp/algorithm/abs_max.hpp"
        "mc/dsp/algorithm/abs_max.cpp"
        "mc/dsp/algorithm/corrcoef.hpp"
        "mc/dsp/algorithm/corrcoef.cpp"
        "mc/dsp/algorithm/down_sample.hpp"
        "mc/dsp/algorithm/mean.hpp"
        "mc/dsp/algorithm/median.hpp"
        "mc/dsp/algorithm/mode.hpp"
        "mc/dsp/algorithm/relative_error.hpp"
        "mc/dsp/algorithm/relative_error.cpp"
        "mc/dsp/algorithm/rms_error.hpp"
        "mc/dsp/algorithm/rms_error.cpp"
        "mc/dsp/algorithm/signal_extension.hpp"
        "mc/dsp/algorithm/signal_extension.cpp"
        "mc/dsp/algorithm/spectral_convolution.hpp"
        "mc/dsp/algorithm/spectral_convolution.cpp"
        "mc/dsp/algorithm/spectral_correlation.hpp"
        "mc/dsp/algorithm/spectral_correlation.cpp"
        "mc/dsp/algorithm/up_sample_even.hpp"
        "mc/dsp/algorithm/up_sample.hpp"

        "mc/dsp/convolution.hpp"
        "mc/dsp/convolution/convolution_method.hpp"
        "mc/dsp/convolution/fft_convolver.cpp"
        "mc/dsp/convolution/fft_convolver.hpp"
        "mc/dsp/convolution/overlap_save_convolver.cpp"
        "mc/dsp/convolution/overlap_save_convolver.hpp"

        "mc/dsp/fft.hpp"
        "mc/dsp/fft/fft.hpp"
        "mc/dsp/fft/fft.cpp"
        "mc/dsp/fft/rfft.hpp"
        "mc/dsp/fft/rfft.cpp"

        "mc/dsp/fft/backend/pffft.hpp"
        "mc/dsp/fft/backend/pffft.cpp"

        "mc/dsp/wavelet.hpp"
        "mc/dsp/wavelet/wavelet.cpp"
        "mc/dsp/wavelet/wavelet.hpp"

        "mc/dsp/wavelet/filters/coif.hpp"
        "mc/dsp/wavelet/filters/daubechies.hpp"
        "mc/dsp/wavelet/filters/sym.hpp"

        "mc/dsp/wavelet/transform/common.cpp"
        "mc/dsp/wavelet/transform/common.hpp"
        "mc/dsp/wavelet/transform/wavelet_packet_transform.cpp"
        "mc/dsp/wavelet/transform/wavelet_packet_transform.hpp"
        "mc/dsp/wavelet/transform/wavelet_transform_2d.cpp"
        "mc/dsp/wavelet/transform/wavelet_transform_2d.hpp"
        "mc/dsp/wavelet/transform/wavelet_transform.cpp"
        "mc/dsp/wavelet/transform/wavelet_transform.hpp"
        "mc/dsp/wavelet/transform/wavelet_tree.cpp"
        "mc/dsp/wavelet/transform/wavelet_tree.hpp"


        "mc/dsp//wavelet/widget/denoise.hpp"
        "mc/dsp//wavelet/widget/denoise.cpp"
        "mc/dsp//wavelet/widget/tempo_detect.cpp"
        "mc/dsp//wavelet/widget/tempo_detect.hpp"

        "mc/testing/test.cpp"
        "mc/testing/test.hpp"
)

find_package(Catch2 REQUIRED)
include(Catch)

add_executable(mc-fft_tests)
catch_discover_tests(mc-fft_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
target_link_libraries(mc-fft_tests PRIVATE mc-wavelet Catch2::Catch2)
target_sources(mc-fft_tests
    PRIVATE
        "mc/dsp/convolution/convolute.test.cpp"
        "mc/dsp/convolution/overlap_save_convolver.test.cpp"

        "mc/dsp/fft/rfft.test.cpp"
)


add_executable(mc-wavelet_tests)
catch_discover_tests(mc-wavelet_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
target_link_libraries(mc-wavelet_tests PRIVATE mc-wavelet Catch2::Catch2)
target_sources(mc-wavelet_tests
    PRIVATE
        "mc/dsp/wavelet/wavelet.test.cpp"
        "mc/dsp/wavelet/transform/wavelet_packet_transform.test.cpp"
        "mc/dsp/wavelet/transform/wavelet_transform.test.cpp"
        "mc/dsp/wavelet/transform/wavelet_transform_2d.test.cpp"

        "mc/dsp/wavelet/widget/tempo_detect.test.cpp"
)
