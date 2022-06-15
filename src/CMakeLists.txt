
add_library(mc_dsp_wavelet)
target_include_directories(mc_dsp_wavelet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mc_dsp_wavelet
    PUBLIC
    mc::Span
    FFTW3::fftw3f
    kissfft::kissfft
    pffft::pffft
    fmt::fmt
    mc::CompilerOptions
    mc::CompilerWarnings
    mc::CodeCoverage
)

target_sources(mc_dsp_wavelet
    PRIVATE
    "mc/algorithm.hpp"
    "mc/array.hpp"
    "mc/bit.hpp"
    "mc/cassert.hpp"
    "mc/cfloat.hpp"
    "mc/cmath.hpp"
    "mc/complex.hpp"
    "mc/cstddef.hpp"
    "mc/cstdlib.hpp"
    "mc/cstring.hpp"
    "mc/filesystem.hpp"
    "mc/format.hpp"
    "mc/fstream.hpp"
    "mc/functional.hpp"
    "mc/iterator.hpp"
    "mc/limits.hpp"
    "mc/memory.hpp"
    "mc/numbers.hpp"
    "mc/numeric.hpp"
    "mc/optional.hpp"
    "mc/random.hpp"
    "mc/span.hpp"
    "mc/sstream.hpp"
    "mc/string.hpp"
    "mc/string_view.hpp"
    "mc/type_traits.hpp"
    "mc/unordered_map.hpp"
    "mc/utility.hpp"
    "mc/vector.hpp"

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
    "mc/dsp/wavelets/common.cpp"
    "mc/dsp/wavelets/common.hpp"
    "mc/dsp/wavelets/ContinuousWaveletTransform.cpp"
    "mc/dsp/wavelets/ContinuousWaveletTransform.hpp"
    "mc/dsp/wavelets/Denoise.hpp"
    "mc/dsp/wavelets/Denoise.cpp"
    "mc/dsp/wavelets/SignalExtension.hpp"
    "mc/dsp/wavelets/Wavelet.cpp"
    "mc/dsp/wavelets/Wavelet.hpp"
    "mc/dsp/wavelets/WaveletPacketTransform.cpp"
    "mc/dsp/wavelets/WaveletPacketTransform.hpp"
    "mc/dsp/wavelets/WaveletTransform2D.cpp"
    "mc/dsp/wavelets/WaveletTransform2D.hpp"
    "mc/dsp/wavelets/WaveletTransform.cpp"
    "mc/dsp/wavelets/WaveletTransform.hpp"
    "mc/dsp/wavelets/WaveletTree.cpp"
    "mc/dsp/wavelets/WaveletTree.hpp"

    "mc/dsp/wavelets/filters/coif.hpp"
    "mc/dsp/wavelets/filters/daubechies.hpp"
    "mc/dsp/wavelets/filters/h.hpp"
    "mc/dsp/wavelets/filters/meyer.hpp"
    "mc/dsp/wavelets/filters/sym.hpp"

    "mc/dsp/widget/TempoDetect.cpp"
    "mc/dsp/widget/TempoDetect.hpp"
)

add_executable(mc_dsp_wavelet_tests "mc/bit.test.cpp" "mc/dsp/widget/TempoDetect.test.cpp")
catch_discover_tests(mc_dsp_wavelet_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
target_link_libraries(mc_dsp_wavelet_tests PRIVATE mc_dsp_wavelet Catch2::Catch2)