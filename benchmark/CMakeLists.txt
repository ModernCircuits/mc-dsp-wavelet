# SPDX-License-Identifier: BSL-1.0

project(mc-wavelet-benchmark)

find_package(benchmark REQUIRED)

add_executable(benchmark-fft src/fft.cpp)
target_link_libraries(benchmark-fft benchmark::benchmark mc::wavelet)
