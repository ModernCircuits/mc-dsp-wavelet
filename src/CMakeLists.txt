cmake_minimum_required(VERSION 3.15)
project(mc-dsp-wavelet)

find_package(mc-dsp REQUIRED)
find_package(FFTW3f REQUIRED)
find_package(kissfft REQUIRED)
find_package(pffft REQUIRED)

add_library(mc-dsp-wavelet)
add_library(mc::dsp-wavelet ALIAS mc-dsp-wavelet)
target_include_directories(mc-dsp-wavelet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mc-dsp-wavelet
    PUBLIC
        mc-dsp::mc-dsp

        FFTW3::fftw3f
        kissfft::kissfft
        pffft::pffft

        mc::CompilerOptions
        mc::CompilerWarnings
        mc::CodeCoverage
)

target_sources(mc-dsp-wavelet
    PRIVATE
        "mc/core/AudioFile.h"

        "mc/dsp/convolution.hpp"
        "mc/dsp/convolution/FFTConvolver.cpp"
        "mc/dsp/convolution/FFTConvolver.hpp"
        "mc/dsp/convolution/ConvolutionMethod.hpp"
        "mc/dsp/convolution/OverlapSaveConvolver.cpp"
        "mc/dsp/convolution/OverlapSaveConvolver.hpp"

        "mc/dsp/fft.hpp"
        "mc/dsp/fft/FFT.hpp"
        "mc/dsp/fft/FFT.cpp"
        "mc/dsp/fft/STFT.hpp"
        "mc/dsp/fft/STFT.cpp"
        "mc/dsp/fft/Window.hpp"
        "mc/dsp/fft/WindowFunction.hpp"
        "mc/dsp/fft/WindowFunction.cpp"

        "mc/dsp/wavelets.hpp"
        "mc/dsp/wavelets/Wavelet.cpp"
        "mc/dsp/wavelets/Wavelet.hpp"

        "mc/dsp/wavelets/filters/coif.hpp"
        "mc/dsp/wavelets/filters/daubechies.hpp"
        "mc/dsp/wavelets/filters/sym.hpp"

        "mc/dsp/wavelets/transform/common.cpp"
        "mc/dsp/wavelets/transform/common.hpp"
        "mc/dsp/wavelets/transform/wavelet_packet_transform.cpp"
        "mc/dsp/wavelets/transform/wavelet_packet_transform.hpp"
        "mc/dsp/wavelets/transform/wavelet_transform_2d.cpp"
        "mc/dsp/wavelets/transform/wavelet_transform_2d.hpp"
        "mc/dsp/wavelets/transform/wavelet_transform.cpp"
        "mc/dsp/wavelets/transform/wavelet_transform.hpp"
        "mc/dsp/wavelets/transform/wavelet_tree.cpp"
        "mc/dsp/wavelets/transform/wavelet_tree.hpp"


        "mc/dsp//wavelets/widget/denoise.hpp"
        "mc/dsp//wavelets/widget/denoise.cpp"
        "mc/dsp//wavelets/widget/tempo_detect.cpp"
        "mc/dsp//wavelets/widget/tempo_detect.hpp"

        "mc/testing/test.cpp"
        "mc/testing/test.hpp"
)

find_package(Catch2 REQUIRED)
include(Catch)

add_executable(mc-dsp-wavelet_tests)
catch_discover_tests(mc-dsp-wavelet_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
target_link_libraries(mc-dsp-wavelet_tests PRIVATE mc-dsp-wavelet Catch2::Catch2)
target_sources(mc-dsp-wavelet_tests
    PRIVATE
        "mc/dsp/convolution/convolute.test.cpp"
        "mc/dsp/convolution/OverlapSaveConvolver.test.cpp"

        "mc/dsp/fft/FFT.test.cpp"
        "mc/dsp/fft/Window.test.cpp"

        "mc/dsp/wavelets/wavelet.test.cpp"
        "mc/dsp/wavelets/transform/wavelet_packet_transform.test.cpp"
        "mc/dsp/wavelets/transform/wavelet_transform.test.cpp"
        "mc/dsp/wavelets/transform/wavelet_transform_2d.test.cpp"

        "mc/dsp/wavelets/widget/tempo_detect.test.cpp"
)
