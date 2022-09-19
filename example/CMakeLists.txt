# SPDX-License-Identifier: BSL-1.0

function(mc_dsp_wavelet_add_example _target)
    add_executable(${_target} "${_target}.cpp")
    add_test(NAME ${_target} COMMAND ${_target} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    target_link_libraries(${_target} PRIVATE mc::dsp-wavelet mc::CodeCoverage)
    target_compile_definitions(${_target} PRIVATE _CRT_SECURE_NO_WARNINGS)
endfunction()

mc_dsp_wavelet_add_example(wavelet)
