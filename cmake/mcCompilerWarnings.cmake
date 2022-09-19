# SPDX-License-Identifier: BSL-1.0

add_library(mc_compiler_warnings INTERFACE)
add_library(mc::CompilerWarnings ALIAS mc_compiler_warnings)

if(MSVC)
    if(MC_BUILD_WERROR)
        target_compile_options(mc_compiler_warnings INTERFACE /WX)
    endif(MC_BUILD_WERROR)

    target_compile_options(mc_compiler_warnings INTERFACE /W3)
else()
    if(MC_BUILD_WERROR)
        target_compile_options(mc_compiler_warnings INTERFACE -Werror)
    endif(MC_BUILD_WERROR)

    target_compile_options(mc_compiler_warnings
        INTERFACE
            -Wall
            -Wextra
            -Wpedantic
            -Wcast-align
            -Wshadow
            -Wunused-parameter
            -Wnarrowing
            -Wshadow
            -Wconversion

        # -Wno-sign-compare
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
            -Wshadow-all

            # -Weverything
            # -Wno-c++98-compat-pedantic
            # -Wno-documentation-unknown-command
            # -Wno-newline-eof
            # -Wno-float-equal
            # -Wno-global-constructors
            # -Wno-padded
            # -Wno-missing-noreturn
            # -Wno-disabled-macro-expansion
            # -Wno-ctad-maybe-unsupported
            # -Wno-unused-member-function
            # -Wno-old-style-cast
            # -Wno-implicit-int-float-conversion
        >
        $<$<CXX_COMPILER_ID:AppleClang>:
        >
        $<$<CXX_COMPILER_ID:GNU>:
        >
    )
endif(MSVC)
