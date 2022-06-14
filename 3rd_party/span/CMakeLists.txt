project(span)

add_library(${PROJECT_NAME} INTERFACE)
add_library(mc::Span ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR})
target_compile_definitions(${PROJECT_NAME} INTERFACE TCB_SPAN_NAMESPACE_NAME=mc)
