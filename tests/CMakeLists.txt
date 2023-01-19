function(set_test id)
    cmake_parse_arguments(SET_TEST "" "EXT" "" ${ARGN})
    if (NOT DEFINED SET_TEST_EXT)
        set(SET_TEST_EXT "c")
    endif()
    add_executable(
        test-${id}
        test-${id}.${SET_TEST_EXT}
    )
    add_common_compile_options(test-${id})
    target_include_directories(
        test-${id} PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_BINARY_DIR}
    )
    target_link_libraries(
        test-${id} PRIVATE
        gtoint
    )
    if(
        (${CMAKE_C_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU") OR
        (${CMAKE_C_COMPILER_ID} STREQUAL "Clang") OR (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang") OR
        (((${CMAKE_C_COMPILER_ID} STREQUAL "Intel") OR (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")) AND NOT (${CMAKE_SYSTEM_NAME} STREQUAL "Windows"))
    ) # NOTE: Generator expressions in target_link_libraries do not work correctly as of cmake 3.20.0.
        target_link_libraries(
            test-${id} PRIVATE
            m
        )
    endif()
    add_test(
        NAME ${id}
        COMMAND $<TARGET_FILE:test-${id}>
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
    )
endfunction()

set_test(function-boys)
set_test(function-bessel)

set_test(spherical)

set_test(integral-oi)
set_test(integral-oi-cartesian)
set_test(integral-oi-spherical)

set_test(integral-kei)
set_test(integral-kei-cartesian)
set_test(integral-kei-spherical)

set_test(integral-nai)
set_test(integral-nai-cartesian)
set_test(integral-nai-spherical)

set_test(integral-mmi)

set_test(integral-ecp-1-0)
set_test(integral-ecp-1-1)
set_test(integral-ecp-1-2)
set_test(integral-ecp-2)
set_test(integral-ecp-1-0-cartesian)
set_test(integral-ecp-1-1-cartesian)
set_test(integral-ecp-1-2-cartesian)
set_test(integral-ecp-2-0-cartesian)
set_test(integral-ecp-2-1-cartesian)
set_test(integral-ecp-2-2-cartesian)
set_test(integral-ecp-cartesian)
set_test(integral-ecp-spherical)

set_test(integral-eri)
set_test(integral-eri-cartesian)
set_test(integral-eri-spherical)

if(Fortran)
    set_test(interface-fortran EXT f90)
endif()
