add_library(${CMAKE_PROJECT_NAME})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(${CMAKE_PROJECT_NAME}
    PUBLIC
        lt::Span
        FFTW3::fftw3
        coz::coz
        lt::CompilerOptions
        lt::CompilerWarnings
        lt::CodeCoverage
)
target_sources(${CMAKE_PROJECT_NAME}
    PRIVATE
        core/AudioFile.h

        dsp/convolution/Convolution.cpp
        dsp/convolution/Convolution.hpp
        dsp/convolution/ConvolutionMethod.hpp
        dsp/convolution/OverlapSave.cpp
        dsp/convolution/OverlapSave.hpp

        dsp/fft/FFT.hpp
        dsp/fft/FFT.cpp

        dsp/wavelets.hpp
        dsp/wavelets/Denoise.hpp
        dsp/wavelets/Denoise.cpp
        dsp/wavelets/SignalExtension.hpp
        dsp/wavelets/Wavelet.cpp
        dsp/wavelets/Wavelet.hpp

        dsp/wavelets/filters/coif.hpp
        dsp/wavelets/filters/daubechies.hpp
        dsp/wavelets/filters/h.hpp
        dsp/wavelets/filters/meyer.hpp
        dsp/wavelets/filters/sym.hpp

        dsp/wavelets/transforms/common.cpp
        dsp/wavelets/transforms/common.hpp
        dsp/wavelets/transforms/ComplexWaveletTransform.cpp
        dsp/wavelets/transforms/ComplexWaveletTransform.hpp
        dsp/wavelets/transforms/WaveletPacketTransform.cpp
        dsp/wavelets/transforms/WaveletPacketTransform.hpp
        dsp/wavelets/transforms/WaveletTransform2D.cpp
        dsp/wavelets/transforms/WaveletTransform2D.hpp
        dsp/wavelets/transforms/WaveletTransform.cpp
        dsp/wavelets/transforms/WaveletTransform.hpp
        dsp/wavelets/transforms/WaveletTree.cpp
        dsp/wavelets/transforms/WaveletTree.hpp
)