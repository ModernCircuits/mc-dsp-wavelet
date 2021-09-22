add_library(${CMAKE_PROJECT_NAME})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
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
        AudioFile.h
        OverlapSaveConvolver.cpp
        OverlapSaveConvolver.hpp

        wavelets.hpp
        wavelets/Convolution.cpp
        wavelets/Convolution.hpp
        wavelets/ConvolutionMethod.hpp
        wavelets/cwtmath.cpp
        wavelets/cwtmath.h
        wavelets/denoise.cpp
        wavelets/FFT.hpp
        wavelets/FFT.cpp
        wavelets/real.cpp
        wavelets/real.h
        wavelets/SignalExtension.hpp
        wavelets/waux.cpp
        wavelets/waux.h
        wavelets/wauxlib.h
        wavelets/wavefilt.cpp
        wavelets/wavefilt.h
        wavelets/wavefunc.cpp
        wavelets/wavefunc.h
        wavelets/Wavelet.cpp
        wavelets/Wavelet.hpp
        wavelets/transforms/ComplexWaveletTransform.cpp
        wavelets/transforms/ComplexWaveletTransform.hpp
        wavelets/transforms/WaveletPacketTransform.cpp
        wavelets/transforms/WaveletPacketTransform.hpp
        wavelets/transforms/WaveletTransform2D.cpp
        wavelets/transforms/WaveletTransform2D.hpp
        wavelets/transforms/WaveletTransform.cpp
        wavelets/transforms/WaveletTransform.hpp
        wavelets/transforms/WaveletTree.cpp
        wavelets/transforms/WaveletTree.hpp
        wavelets/wtmath.cpp
        wavelets/wtmath.h
)

add_library(wauxlib)
target_sources(wauxlib PRIVATE wavelets/denoise.cpp wavelets/waux.h wavelets/waux.cpp)
target_link_libraries(wauxlib PRIVATE wavelib)
target_include_directories(wauxlib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})