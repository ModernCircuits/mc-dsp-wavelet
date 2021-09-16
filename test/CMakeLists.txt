add_executable(cwttest cwttest.cpp)
target_link_libraries(cwttest wavelib)

add_executable(dwttest dwttest.cpp)
target_link_libraries(dwttest wavelib)

add_executable(swttest swttest.cpp)
target_link_libraries(swttest wavelib)

add_executable(modwttest modwttest.cpp)
target_link_libraries(modwttest wavelib)

add_executable(dwpttest dwpttest.cpp)
target_link_libraries(dwpttest wavelib)

add_executable(wtreetest wtreetest.cpp)
target_link_libraries(wtreetest wavelib)

add_executable(denoisetest denoisetest.cpp)
target_link_libraries(denoisetest wauxlib wavelib)

add_executable(modwtdenoisetest modwtdenoisetest.cpp)
target_link_libraries(modwtdenoisetest wauxlib wavelib)

add_executable(dwt2test dwt2test.cpp)
target_link_libraries(dwt2test wavelib)

add_executable(swt2test swt2test.cpp)
target_link_libraries(swt2test wavelib)

add_executable(modwt2test modwt2test.cpp)
target_link_libraries(modwt2test wavelib)

set_target_properties(cwttest dwttest swttest modwttest dwpttest wtreetest denoisetest modwtdenoisetest dwt2test swt2test modwt2test
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
        )
