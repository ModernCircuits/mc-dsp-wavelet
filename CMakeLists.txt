cmake_minimum_required(VERSION 3.15)
cmake_policy(VERSION 3.15)
project(wavelib VERSION 0.1.0)


# Only do these if this is the main project,
# and not if it is included through add_subdirectory.
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Caches build artifacts for faster builds
    find_program(CCACHE ccache)
    if (CCACHE)
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
    endif ()

    # Set C11 with no extensions enabled. Fail if the
    # standard is not supported.
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)

    # Set C++17 with no extensions enabled. Fail if the
    # standard is not supported.
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # misc
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif ()

option(BUILD_UT "Enable Unit test" ON)

add_library(${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(${PROJECT_NAME} PUBLIC m)
target_sources(${PROJECT_NAME}
    PRIVATE
        src/conv.h
        src/conv.cpp
        src/cwt.h
        src/cwt.cpp
        src/cwtmath.h
        src/cwtmath.cpp
        src/hsfft.h
        src/hsfft.cpp
        src/real.h
        src/real.cpp
        src/wavefilt.h
        src/wavefilt.cpp
        src/wavefunc.h
        src/wavefunc.cpp
        src/wavelib.cpp
        src/wtmath.h
        src/wtmath.cpp
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
# target_sources(${PROJECT_NAME}-tests PRIVATE src/PluginProcessorTest.cpppp)
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