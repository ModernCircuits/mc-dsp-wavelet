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
        conv.h
        conv.cpp
        cwt.h
        cwt.cpp
        cwtmath.h
        cwtmath.cpp
        hsfft.h
        hsfft.cpp
        overlap_save_convolver.hpp
        overlap_save_convolver.cpp
        real.h
        real.cpp
        wavefilt.h
        wavefilt.cpp
        wavefunc.h
        wavefunc.cpp
        wavelib.cpp
        wtmath.h
        wtmath.cpp
)

add_library(wauxlib)
target_sources(wauxlib PRIVATE denoise.cpp waux.h waux.cpp)
target_link_libraries(wauxlib PRIVATE wavelib)
target_include_directories(wauxlib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})