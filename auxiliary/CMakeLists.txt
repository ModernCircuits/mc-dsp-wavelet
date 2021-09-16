project(wauxlib)

add_library(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE denoise.cpp waux.h waux.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE wavelib)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
