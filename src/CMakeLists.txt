# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.15)
project(mc-wavelet)

find_package(mc-core REQUIRED)
find_package(FFTW3f REQUIRED)
find_package(pffft REQUIRED)

add_library(mc-wavelet)
add_library(mc::wavelet ALIAS mc-wavelet)
target_include_directories(mc-wavelet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mc-wavelet
    PUBLIC
        mc-core::mc-core

        FFTW3::fftw3f
        pffft::pffft

        mc::CompilerOptions
        mc::CompilerWarnings
        mc::CodeCoverage
)

target_sources(mc-wavelet
    PRIVATE
        "mc/audio/audio_file.hpp"

        "mc/algorithm/abs_max.hpp"
        "mc/algorithm/abs_max.cpp"
        "mc/algorithm/corrcoef.hpp"
        "mc/algorithm/corrcoef.cpp"
        "mc/algorithm/down_sample.hpp"
        "mc/algorithm/mean.hpp"
        "mc/algorithm/median.hpp"
        "mc/algorithm/mode.hpp"
        "mc/algorithm/relative_error.hpp"
        "mc/algorithm/relative_error.cpp"
        "mc/algorithm/rms_error.hpp"
        "mc/algorithm/rms_error.cpp"
        "mc/algorithm/signal_extension.hpp"
        "mc/algorithm/signal_extension.cpp"
        "mc/algorithm/spectral_convolution.hpp"
        "mc/algorithm/spectral_convolution.cpp"
        "mc/algorithm/spectral_correlation.hpp"
        "mc/algorithm/spectral_correlation.cpp"
        "mc/algorithm/up_sample_even.hpp"
        "mc/algorithm/up_sample.hpp"

        "mc/fft/convolution.hpp"
        "mc/fft/convolution/convolution_method.hpp"
        "mc/fft/convolution/fft_convolver.cpp"
        "mc/fft/convolution/fft_convolver.hpp"
        "mc/fft/convolution/overlap_save_convolver.cpp"
        "mc/fft/convolution/overlap_save_convolver.hpp"

        "mc/fft.hpp"
        "mc/fft/fft.hpp"
        "mc/fft/fft.cpp"
        "mc/fft/rfft.hpp"
        "mc/fft/rfft.cpp"

        "mc/fft/backend/pffft.hpp"
        "mc/fft/backend/pffft.cpp"

        "mc/wavelet.hpp"
        "mc/wavelet/wavelet.cpp"
        "mc/wavelet/wavelet.hpp"

        "mc/wavelet/filters/coif.hpp"
        "mc/wavelet/filters/daubechies.hpp"
        "mc/wavelet/filters/sym.hpp"

        "mc/wavelet/transform/common.cpp"
        "mc/wavelet/transform/common.hpp"
        "mc/wavelet/transform/wavelet_packet_transform.cpp"
        "mc/wavelet/transform/wavelet_packet_transform.hpp"
        "mc/wavelet/transform/wavelet_transform_2d.cpp"
        "mc/wavelet/transform/wavelet_transform_2d.hpp"
        "mc/wavelet/transform/wavelet_transform.cpp"
        "mc/wavelet/transform/wavelet_transform.hpp"
        "mc/wavelet/transform/wavelet_tree.cpp"
        "mc/wavelet/transform/wavelet_tree.hpp"


        "mc/wavelet/widget/denoise.hpp"
        "mc/wavelet/widget/denoise.cpp"
        "mc/wavelet/widget/tempo_detect.cpp"
        "mc/wavelet/widget/tempo_detect.hpp"

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
        "mc/fft/convolution/convolute.test.cpp"
        "mc/fft/convolution/overlap_save_convolver.test.cpp"

        "mc/fft/rfft.test.cpp"
)


add_executable(mc-wavelet_tests)
catch_discover_tests(mc-wavelet_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
target_link_libraries(mc-wavelet_tests PRIVATE mc-wavelet Catch2::Catch2)
target_sources(mc-wavelet_tests
    PRIVATE
        "mc/wavelet/wavelet.test.cpp"
        "mc/wavelet/transform/wavelet_packet_transform.test.cpp"
        "mc/wavelet/transform/wavelet_transform.test.cpp"
        "mc/wavelet/transform/wavelet_transform_2d.test.cpp"

        "mc/wavelet/widget/tempo_detect.test.cpp"
)
