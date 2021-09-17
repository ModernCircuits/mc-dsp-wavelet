# Code Coverage Configuration
add_library(lt_coverage INTERFACE)
add_library(lt::CodeCoverage ALIAS lt_coverage)

if(LT_BUILD_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(lt_coverage INTERFACE -O0 -g --coverage)
  target_link_libraries(lt_coverage INTERFACE --coverage)
  target_link_options(lt_coverage INTERFACE --coverage)
endif()