add_library(lt_compiler_warnings INTERFACE)
add_library(lt::CompilerWarnings ALIAS lt_compiler_warnings)

if (MSVC)
    if (LT_BUILD_WERROR)
        target_compile_options(lt_compiler_warnings INTERFACE /WX)
    endif (LT_BUILD_WERROR)
    target_compile_options(lt_compiler_warnings INTERFACE /W3)
else ()
    if (LT_BUILD_WERROR)
        target_compile_options(lt_compiler_warnings INTERFACE -Werror)
    endif (LT_BUILD_WERROR)
    target_compile_options(lt_compiler_warnings
            INTERFACE
                    -Wall
                    -Wextra
                    -Wpedantic
                    -Wcast-align
                    -Wshadow
                    -Wunused-parameter
                    -Wnarrowing
                    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
                        #-Weverything
                    >
                    $<$<CXX_COMPILER_ID:AppleClang>:
                        -Wno-poison-system-directories
                    >
                    $<$<CXX_COMPILER_ID:GNU>:
                        #-Wlogical-op
                    >
            )

endif (MSVC)
