# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.15)
project(mc-wavelet-dev VERSION 0.1.0)

find_program(CCACHE ccache)
if(CCACHE)
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE})
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if (NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.14")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CTest)
enable_testing()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(mcCompilerOptions)
include(mcCompilerWarnings)
include(mcCodeCoverage)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})
include(mcConan)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)
add_subdirectory(benchmark EXCLUDE_FROM_ALL)

find_package(Catch2 REQUIRED)
include(Catch)

add_executable(mc-fft_tests)
catch_discover_tests(mc-fft_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
target_link_libraries(mc-fft_tests PRIVATE mc::fft mc::testing Catch2::Catch2)
target_sources(mc-fft_tests
    PRIVATE
        "src/mc/fft/convolution/convolute.test.cpp"
        "src/mc/fft/convolution/overlap_save_convolver.test.cpp"

        "src/mc/fft/transform/rfft.test.cpp"
)


add_executable(mc-wavelet_tests)
catch_discover_tests(mc-wavelet_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
target_link_libraries(mc-wavelet_tests PRIVATE mc::wavelet mc::testing Catch2::Catch2)
target_sources(mc-wavelet_tests
    PRIVATE
        "src/mc/wavelet/wavelet.test.cpp"
        "src/mc/wavelet/transform/wavelet_packet_transform.test.cpp"
        "src/mc/wavelet/transform/wavelet_transform.test.cpp"
        "src/mc/wavelet/transform/wavelet_transform_2d.test.cpp"

        "src/mc/wavelet/widget/tempo_detect.test.cpp"
)
