cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)
project(wavelib VERSION 0.1.0)

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
    include(ltCompilerOptions)
    include(ltCompilerWarnings)
    include(ltCodeCoverage)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})

if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
        "${CMAKE_BINARY_DIR}/conan.cmake"
        TLS_VERIFY ON)
endif()

include(${CMAKE_BINARY_DIR}/conan.cmake)
conan_add_remote(NAME conancenter URL https://center.conan.io)
conan_cmake_configure(
    REQUIRES
    fftw/3.3.9
    fmt/8.0.1
    kissfft/131.1.0
    pffft/cci.20210511
    GENERATORS
    cmake_find_package
    OPTIONS
    fftw:precision=single
    fftw:simd=False
    kissfft:datatype=float
    pffft:disable_simd=True
)

conan_cmake_autodetect(settings)
conan_cmake_install(PATH_OR_REFERENCE . BUILD outdated REMOTE conancenter SETTINGS ${settings})

find_package(FFTW3f REQUIRED)
find_package(fmt REQUIRED)
find_package(kissfft REQUIRED)
find_package(pffft REQUIRED)

add_subdirectory(3rd_party/span)

add_subdirectory(lt)
add_subdirectory(test)
