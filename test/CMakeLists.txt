add_library(mc_dsp_wavelet_testdata ${CMAKE_SOURCE_DIR}/mc/testing/test.cpp ${CMAKE_SOURCE_DIR}/mc/testing/test.hpp)
target_include_directories(mc_dsp_wavelet_testdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mc_dsp_wavelet_testdata PUBLIC mc_dsp_wavelet)
add_custom_command(
    TARGET
    mc_dsp_wavelet_testdata
    COMMAND
    ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/testData
    $<TARGET_FILE_DIR:mc_dsp_wavelet_testdata>/testData
    POST_BUILD
)

function(mc_dsp_wavelet_add_test _target)
    add_executable(${_target} "${_target}.cpp")
    add_test(NAME ${_target} COMMAND ${_target})
    target_link_libraries(${_target} PRIVATE mc_dsp_wavelet mc_dsp_wavelet_testdata mc::CodeCoverage)
endfunction()

add_executable(bpm_detector "bpm_detector.cpp")
target_link_libraries(bpm_detector PRIVATE mc_dsp_wavelet mc_dsp_wavelet_testdata mc::CodeCoverage)

mc_dsp_wavelet_add_test(FFTTest)
mc_dsp_wavelet_add_test(DirectConvolverTest)
mc_dsp_wavelet_add_test(OverlapSaveConvolverTest)
mc_dsp_wavelet_add_test(WindowTest)
mc_dsp_wavelet_add_test(StdHeadersTest)

mc_dsp_wavelet_add_test(cwttest)
mc_dsp_wavelet_add_test(dwttest)
mc_dsp_wavelet_add_test(swttest)
mc_dsp_wavelet_add_test(modwttest)
mc_dsp_wavelet_add_test(dwpttest)
mc_dsp_wavelet_add_test(wtreetest)
mc_dsp_wavelet_add_test(dwt2test)
mc_dsp_wavelet_add_test(swt2test)
mc_dsp_wavelet_add_test(modwt2test)
mc_dsp_wavelet_add_test(tst_dwt)
mc_dsp_wavelet_add_test(dwt_reconstruction)
mc_dsp_wavelet_add_test(swt_reconstruction)
mc_dsp_wavelet_add_test(denoisetest)
mc_dsp_wavelet_add_test(modwtdenoisetest)
