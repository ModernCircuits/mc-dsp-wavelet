cmake_minimum_required(VERSION 3.15)
project(wavelib C)

option(BUILD_UT "Enable Unit test" ON)

add_library(${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(${PROJECT_NAME} PUBLIC m)
target_sources(${PROJECT_NAME}
    PRIVATE
        src/conv.h
        src/conv.c
        src/cwt.h
        src/cwt.c
        src/cwtmath.h
        src/cwtmath.c
        src/hsfft.h
        src/hsfft.c
        src/real.h
        src/real.c
        src/wavefilt.h
        src/wavefilt.c
        src/wavefunc.h
        src/wavefunc.c
        src/wavelib.c
        src/wtmath.h
        src/wtmath.c
)

add_subdirectory(auxiliary)
add_subdirectory(test)

if(BUILD_UT)
    include(CTest)
    enable_testing()
    enable_language(CXX)
    add_subdirectory(unitTests)
endif()


# add_executable(${PROJECT_NAME}-tests)
# add_test(NAME "PluginTests" COMMAND ${PROJECT_NAME}-tests)
# target_sources(${PROJECT_NAME}-tests PRIVATE src/PluginProcessorTest.cpp)
# target_link_libraries(${PROJECT_NAME}-tests PRIVATE ${PROJECT_NAME})
# add_custom_command(
#     TARGET
#         ${PROJECT_NAME}-tests
#     COMMAND
#         ${CMAKE_COMMAND} -E copy_directory
#         ${CMAKE_CURRENT_SOURCE_DIR}/testData
#         $<TARGET_FILE_DIR:${PROJECT_NAME}-tests>/testData
#     POST_BUILD
# )