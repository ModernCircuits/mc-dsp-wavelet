add_library(${CMAKE_PROJECT_NAME})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC lt::CodeCoverage)
target_sources(${CMAKE_PROJECT_NAME}
    PRIVATE
        conv.h
        conv.cpp
        cwt.h
        cwt.cpp
        cwtmath.h
        cwtmath.cpp
        hsfft.h
        hsfft.cpp
        real.h
        real.cpp
        wavefilt.h
        wavefilt.cpp
        wavefunc.h
        wavefunc.cpp
        wavelib.cpp
        wtmath.h
        wtmath.cpp
)
