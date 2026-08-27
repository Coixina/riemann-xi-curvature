cmake_minimum_required(VERSION 3.15)

project(
    rh_hardy_certifier
    VERSION 0.1.0
    LANGUAGES C
)

# ============================================================
# Language standard
# ============================================================

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# ============================================================
# Options
# ============================================================

option(
    RH_BUILD_VERIFY_LARGE_CONSTANTS
    "Build the rigorous large-t constants verifier"
    ON
)


# ============================================================
# FLINT / Arb
# ============================================================

# FLINT 3.x contains Arb/Acb.
#
# On Debian the Arb headers commonly live in /usr/include/flint,
# so detect both the include directory and library explicitly.

find_path(
    FLINT_INCLUDE_DIR
    NAMES arb.h
    PATH_SUFFIXES flint
)

find_library(
    FLINT_LIBRARY
    NAMES flint
)

if(NOT FLINT_INCLUDE_DIR)
    message(
        FATAL_ERROR
        "FLINT/Arb headers were not found (arb.h)."
    )
endif()

if(NOT FLINT_LIBRARY)
    message(
        FATAL_ERROR
        "FLINT library was not found (libflint)."
    )
endif()

message(
    STATUS
    "FLINT include dir : ${FLINT_INCLUDE_DIR}"
)

message(
    STATUS
    "FLINT library     : ${FLINT_LIBRARY}"
)


add_library(
    rh_flint
    INTERFACE
)

target_include_directories(
    rh_flint
    INTERFACE
    ${FLINT_INCLUDE_DIR}
)

target_link_libraries(
    rh_flint
    INTERFACE
    ${FLINT_LIBRARY}
)


# ============================================================
# Compiler warnings
# ============================================================

add_library(
    rh_warnings
    INTERFACE
)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")

    target_compile_options(
        rh_warnings
        INTERFACE
        -Wall
        -Wextra
        -Wpedantic
    )

endif()


# ============================================================
# Common project options
# ============================================================

add_library(
    rh_project_options
    INTERFACE
)

target_include_directories(
    rh_project_options
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
    rh_project_options
    INTERFACE
    rh_flint
    rh_warnings
)


# ============================================================
# Hardy finite-certifier core
# ============================================================

#
# IMPORTANT:
#
# Only production modules belonging to the final Hardy
# certifier are listed here.
#
# Modules are added only after their definitive Hardy-only
# implementation and tests have been completed.
#

set(
    RH_CORE_SOURCES

    src/interval_geometry.c
    src/lambda_model.c
    src/hardy_jet.c
    src/bilinear_certifier.c
    src/interval_engine.c
    src/certificate_writer.c
)


add_library(
    rh_core
    STATIC
    ${RH_CORE_SOURCES}
)

target_link_libraries(
    rh_core
    PUBLIC
    rh_project_options
)

# ============================================================
# Helper functions
# ============================================================

function(
    add_rh_c_program
    target_name
    source_file
)

    add_executable(
        ${target_name}
        ${source_file}
    )

    target_link_libraries(
        ${target_name}
        PRIVATE
        rh_project_options
    )

endfunction()


function(
    add_rh_core_program
    target_name
    source_file
)

    add_executable(
        ${target_name}
        ${source_file}
    )

    target_link_libraries(
        ${target_name}
        PRIVATE
        rh_core
    )

endfunction()


# ============================================================
# Large-t constants verifier
# ============================================================

if(RH_BUILD_VERIFY_LARGE_CONSTANTS)

    add_rh_c_program(
        verify_large_t_constants
        verify/verify_large_t_constants.c
    )

endif()

# ============================================================
# Production finite certifier
# ============================================================

add_rh_core_program(
    rh_certifier
    certifier/rh_certifier.c
)

# ============================================================
# Finite-certificate verifier
# ============================================================

add_rh_core_program(
    verify_certificate
    certifier/verify_certificate.c
)

# ============================================================
# CTest
# ============================================================

include(CTest)

if(BUILD_TESTING)

    add_rh_core_program(
        test_interval_geometry
        tests/test_interval_geometry.c
    )

    add_test(
        NAME interval_geometry
        COMMAND test_interval_geometry
    )


    add_rh_core_program(
        test_lambda_model
        tests/test_lambda_model.c
    )

    add_test(
        NAME lambda_model
        COMMAND test_lambda_model
    )

    add_rh_core_program(
        test_hardy_jet
        tests/test_hardy_jet.c
    )

    add_test(
        NAME hardy_jet
        COMMAND test_hardy_jet
    )

    add_rh_core_program(
        test_bilinear_certifier
        tests/test_bilinear_certifier.c
    )

    add_test(
        NAME bilinear_certifier
        COMMAND test_bilinear_certifier
    )

    add_rh_core_program(
        test_interval_engine
        tests/test_interval_engine.c
    )

    add_test(
        NAME interval_engine
        COMMAND test_interval_engine
    )

        add_rh_core_program(
        test_verifier_refinement
        tests/test_verifier_refinement.c
    )

    add_test(
        NAME verifier_refinement
        COMMAND test_verifier_refinement
    )


    # ========================================================
    # Certificate writer / verifier integration fixtures
    # ========================================================

    set(
        RH_CERTIFICATE_TEST_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/certificate_test_data"
    )

    file(
        MAKE_DIRECTORY
        "${RH_CERTIFICATE_TEST_DIR}"
    )


    add_rh_core_program(
        make_certificate_fixtures
        tests/make_certificate_fixtures.c
    )


    add_test(
        NAME certificate_fixture_generation
        COMMAND
            make_certificate_fixtures
            "${RH_CERTIFICATE_TEST_DIR}"
    )

    set_tests_properties(
        certificate_fixture_generation
        PROPERTIES
        FIXTURES_SETUP rh_certificate_fixture
    )


    add_test(
        NAME certificate_valid
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/valid_certificate.txt"
    )

    set_tests_properties(
        certificate_valid
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
    )

        add_test(
        NAME certificate_valid_320
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/valid_certificate.txt"
            320
    )

    set_tests_properties(
        certificate_valid_320
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
    )


    add_test(
        NAME certificate_valid_384
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/valid_certificate.txt"
            384
    )

    set_tests_properties(
        certificate_valid_384
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
    )


    add_test(
        NAME certificate_valid_512
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/valid_certificate.txt"
            512
    )

    set_tests_properties(
        certificate_valid_512
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
    )

        add_test(
        NAME certificate_512_insufficient_refinement
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/valid_certificate.txt"
            512
            0
    )

    set_tests_properties(
        certificate_512_insufficient_refinement
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
        WILL_FAIL TRUE
    )

    add_test(
        NAME certificate_bad_gap
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/bad_gap.txt"
    )

    set_tests_properties(
        certificate_bad_gap
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
        WILL_FAIL TRUE
    )


    add_test(
        NAME certificate_bad_overlap
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/bad_overlap.txt"
    )

    set_tests_properties(
        certificate_bad_overlap
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
        WILL_FAIL TRUE
    )


    add_test(
        NAME certificate_bad_duplicate
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/bad_duplicate.txt"
    )

    set_tests_properties(
        certificate_bad_duplicate
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
        WILL_FAIL TRUE
    )


    add_test(
        NAME certificate_bad_count
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/bad_count.txt"
    )

    set_tests_properties(
        certificate_bad_count
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
        WILL_FAIL TRUE
    )


    add_test(
        NAME certificate_precision_too_low
        COMMAND
            verify_certificate
            "${RH_CERTIFICATE_TEST_DIR}/valid_certificate.txt"
            128
    )

    set_tests_properties(
        certificate_precision_too_low
        PROPERTIES
        FIXTURES_REQUIRED rh_certificate_fixture
        WILL_FAIL TRUE
    )

    if(TARGET verify_large_t_constants)

        add_test(
            NAME large_t_constants
            COMMAND verify_large_t_constants 256
        )

    endif()

endif()


# ============================================================
# Configuration summary
# ============================================================

message(STATUS "")
message(STATUS "RH Hardy certifier configuration")
message(STATUS "--------------------------------")
message(STATUS "C compiler                  : ${CMAKE_C_COMPILER}")
message(STATUS "Build tests                 : ${BUILD_TESTING}")
message(STATUS "Build large-t verifier      : ${RH_BUILD_VERIFY_LARGE_CONSTANTS}")
message(STATUS "Production core             : Hardy-only")
message(STATUS "")