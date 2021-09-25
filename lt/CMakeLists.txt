add_library(${CMAKE_PROJECT_NAME})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(${CMAKE_PROJECT_NAME}
    PUBLIC
        lt::Span
        Boost::headers
        FFTW3::fftw3
        coz::coz
        lt::CompilerOptions
        lt::CompilerWarnings
        lt::CodeCoverage
)
target_sources(${CMAKE_PROJECT_NAME}
    PRIVATE
        core/AudioFile.h

        dsp/convolution.hpp
        dsp/convolution/FFTConvolver.cpp
        dsp/convolution/FFTConvolver.hpp
        dsp/convolution/ConvolutionMethod.hpp
        dsp/convolution/OverlapSaveConvolver.cpp
        dsp/convolution/OverlapSaveConvolver.hpp

        dsp/fft.hpp
        dsp/fft/FFT.hpp
        dsp/fft/FFT.cpp

        dsp/wavelets.hpp
        dsp/wavelets/common.cpp
        dsp/wavelets/common.hpp
        dsp/wavelets/ComplexWaveletTransform.cpp
        dsp/wavelets/ComplexWaveletTransform.hpp
        dsp/wavelets/Denoise.hpp
        dsp/wavelets/Denoise.cpp
        dsp/wavelets/SignalExtension.hpp
        dsp/wavelets/Wavelet.cpp
        dsp/wavelets/Wavelet.hpp
        dsp/wavelets/WaveletPacketTransform.cpp
        dsp/wavelets/WaveletPacketTransform.hpp
        dsp/wavelets/WaveletTransform2D.cpp
        dsp/wavelets/WaveletTransform2D.hpp
        dsp/wavelets/WaveletTransform.cpp
        dsp/wavelets/WaveletTransform.hpp
        dsp/wavelets/WaveletTree.cpp
        dsp/wavelets/WaveletTree.hpp

        dsp/wavelets/filters/coif.hpp
        dsp/wavelets/filters/daubechies.hpp
        dsp/wavelets/filters/h.hpp
        dsp/wavelets/filters/meyer.hpp
        dsp/wavelets/filters/sym.hpp
)