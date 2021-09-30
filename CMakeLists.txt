cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)
project(wavelib VERSION 0.1.0)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    find_program(CCACHE ccache)
    if (CCACHE)
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
    endif ()

    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)

    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.12")

    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    include(CTest)
    enable_testing()

    include(CheckCXXCompilerFlag)

    if (MSVC)
        # check_cxx_compiler_flag("/arch:AVX" COMPILER_AVX_SUPPORTED)
        # if (COMPILER_AVX_SUPPORTED)
        #     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
        # endif ()
    else ()
        # check_cxx_compiler_flag("-msse2" COMPILER_SSE2_SUPPORTED)
        # if (COMPILER_SSE2_SUPPORTED)
        #     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
        # endif ()
    endif ()

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    include(ltCompilerOptions)
    include(ltCompilerWarnings)
    include(ltCodeCoverage)
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})

if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/v0.16.1/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan.cmake"
                EXPECTED_HASH SHA256=396e16d0f5eabdc6a14afddbcfff62a54a7ee75c6da23f32f7a31bc85db23484
                TLS_VERIFY ON)
endif()

include(${CMAKE_BINARY_DIR}/conan.cmake)
conan_add_remote(NAME conancenter URL https://center.conan.io)
conan_cmake_configure(
    REQUIRES
        boost/1.77.0
        fftw/3.3.9
        fmt/8.0.1
        kissfft/131.1.0
    GENERATORS
        cmake_find_package
    OPTIONS
        boost:without_context=True
        boost:without_contract=True
        boost:without_coroutine=True
        boost:without_fiber=True
        boost:without_filesystem=False
        boost:without_graph=True
        boost:without_graph_parallel=True
        boost:without_iostreams=True
        boost:without_json=False
        boost:without_log=True
        boost:without_regex=True
        boost:without_serialization=True
        boost:without_test=True
        boost:without_type_erasure=True
        boost:without_wave=True
        boost:system_no_deprecated=True
        boost:asio_no_deprecated=True
        boost:filesystem_no_deprecated=True
        fftw:precision=single
        fftw:simd=False
        kissfft:datatype=float
)

conan_cmake_autodetect(settings)

conan_cmake_install(
    PATH_OR_REFERENCE
        .
    BUILD
        outdated
    REMOTE
        conancenter
    SETTINGS
        ${settings}
    SETTINGS
        compiler.cppstd=14
)

find_package(Boost 1.77 REQUIRED COMPONENTS container filesystem json)
find_package(FFTW3f REQUIRED)
find_package(fmt REQUIRED)
find_package(kissfft REQUIRED)

add_subdirectory(3rd_party/span)

add_subdirectory(lt)
add_subdirectory(test)
