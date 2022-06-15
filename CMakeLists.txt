cmake_minimum_required(VERSION 3.15)
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

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    include(mcCompilerOptions)
    include(mcCompilerWarnings)
    include(mcCodeCoverage)

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})
    include(mcConan)

    find_package(Catch2 REQUIRED)
    include(Catch)
endif()

find_package(FFTW3f REQUIRED)
find_package(fmt REQUIRED)
find_package(kissfft REQUIRED)
find_package(pffft REQUIRED)

add_subdirectory(3rd_party/span)
add_subdirectory(src)

# add_subdirectory(test)
