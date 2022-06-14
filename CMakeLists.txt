cmake_minimum_required(VERSION 3.18)
project(mc_dsp_wavelet VERSION 0.1.0)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    find_program(CCACHE ccache)

    if(CCACHE)
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
    endif()

    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.14")

    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    include(CTest)
    enable_testing()

    include(CheckCXXCompilerFlag)

    if(MSVC)
    # check_cxx_compiler_flag("/arch:AVX" COMPILER_AVX_SUPPORTED)
    # if (COMPILER_AVX_SUPPORTED)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
    # endif ()
    else()
        # check_cxx_compiler_flag("-msse2" COMPILER_SSE2_SUPPORTED)
        # if (COMPILER_SSE2_SUPPORTED)
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
        # endif ()
    endif()

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    include(mcCompilerOptions)
    include(mcCompilerWarnings)
    include(mcCodeCoverage)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})
include(mcConan)

find_package(FFTW3f REQUIRED)
find_package(fmt REQUIRED)
find_package(kissfft REQUIRED)
find_package(pffft REQUIRED)

find_package(Catch2 REQUIRED)
include(Catch)

add_subdirectory(3rd_party/span)

add_subdirectory(mc)
add_subdirectory(test)
