add_library(wavelib_testdata helper.hpp helper.cpp)
add_custom_command(
    TARGET
        wavelib_testdata
    COMMAND
        ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/testData
        $<TARGET_FILE_DIR:wavelib_testdata>/testData
    POST_BUILD
)

function(wavelib_add_test _target)
    add_executable(${_target} "${_target}.cpp")
    add_test(NAME ${_target} COMMAND ${_target})
    target_link_libraries(${_target} PRIVATE wavelib wavelib_testdata lt::CodeCoverage)
endfunction()

wavelib_add_test(cwttest)
wavelib_add_test(dwttest)
wavelib_add_test(swttest)
wavelib_add_test(modwttest)
wavelib_add_test(dwpttest)
wavelib_add_test(wtreetest)
wavelib_add_test(dwt2test)
wavelib_add_test(swt2test)
wavelib_add_test(modwt2test)
wavelib_add_test(tst_dwt)
wavelib_add_test(dwt_reconstruction)
wavelib_add_test(swt_reconstruction)

wavelib_add_test(denoisetest)
target_link_libraries(denoisetest PRIVATE wauxlib)

wavelib_add_test(modwtdenoisetest)
target_link_libraries(modwtdenoisetest PRIVATE wauxlib)

