cmake_minimum_required(VERSION 3.14)

project(
    libgtoint
    VERSION 1.0.0.0
    DESCRIPTION "An analytical GTO integral library for C and Fortran"
    HOMEPAGE_URL https://github.com/arithy/libgtoint
    LANGUAGES C
)

set(Fortran ON CACHE BOOL "Support Fortran interface")

set(CMAKE_C_STANDARD 99)

if(Fortran)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
        set(Fortran_SOURCE src/gtoint-fortran.f90)
        set(Fortran_MODULE ${PROJECT_BINARY_DIR}/gtoint.mod)
    else()
        message(FATAL_ERROR "Fortran compiler not found")
    endif()
endif()

function(add_common_compile_options target)
    target_compile_options(
        ${target} PRIVATE
        $<$<COMPILE_LANGUAGE:C>:
            $<$<C_COMPILER_ID:GNU>:
                -Wall -Wextra -Wno-unused-parameter -Wno-misleading-indentation -pedantic -Werror
                $<$<CONFIG:Release>:-ffast-math>
            >
            $<$<C_COMPILER_ID:Clang>:
                -Wall -Wextra -Wno-unused-parameter -Wno-gnu-folding-constant -Wno-gnu-zero-variadic-macro-arguments -pedantic -Werror
                $<$<CONFIG:Release>:-ffast-math>
            >
            $<$<AND:$<C_COMPILER_ID:Intel>,$<NOT:$<PLATFORM_ID:Windows>>>:
                -Wall -Wextra -Wno-unused-parameter -pedantic -Werror
                $<$<CONFIG:Release>:-fp-model fast>
            >
            $<$<C_COMPILER_ID:MSVC>:
                /W4 /WX /wd4100 /wd4706
                $<$<CONFIG:Release>:/fp:fast>
            >
            $<$<AND:$<C_COMPILER_ID:Intel>,$<PLATFORM_ID:Windows>>:
                /W4 /WX
                $<$<CONFIG:Release>:/fp:fast>
            >
        >
        $<$<COMPILE_LANGUAGE:Fortran>:
            $<$<Fortran_COMPILER_ID:GNU>:
                -Wall -Wextra -Wno-unused-dummy-argument -Wno-compare-reals -pedantic -Werror
            >
            $<$<Fortran_COMPILER_ID:Flang>:
                -Wall -Wextra -pedantic -Werror
            >
            $<$<AND:$<Fortran_COMPILER_ID:Intel>,$<NOT:$<PLATFORM_ID:Windows>>>:
                -warn all,errors
            >
            $<$<AND:$<Fortran_COMPILER_ID:Intel>,$<PLATFORM_ID:Windows>>:
                /warn:all,errors
            >
        >
    )
endfunction()

add_library(
    gtoint
    src/array.c
    src/stack.c
    src/cache.c
    src/cache-oi.c
    src/cache-kei.c
    src/cache-nai.c
    src/cache-mmi.c
    src/cache-eri.c
    src/cache-ecp0.c
    src/cache-ecp1.c
    src/cache-ecp2.c
    src/cache-ecp2-w.c
    src/function.c
    src/function-boys.c
    src/function-bessel.c
    src/spherical.c
    src/integral-oi.c
    src/integral-kei.c
    src/integral-nai.c
    src/integral-mmi.c
    src/integral-eri.c
    src/integral-ecp.c
    src/integral-ecp-1-0.c
    src/integral-ecp-1-1.c
    src/integral-ecp-1-2.c
    src/integral-ecp-2.c
    src/integral-ecp-2-spherical.c
    src/integral-ecp-2-angular.c
    src/integral-ecp-2-radial.c
    src/integrator.c
    src/basis.c
    src/ecp.c
    ${Fortran_SOURCE}
)
add_common_compile_options(gtoint)

install(TARGETS gtoint)
install(FILES src/gtoint.h ${Fortran_MODULE} TYPE INCLUDE)

enable_testing()

add_subdirectory(tests)
