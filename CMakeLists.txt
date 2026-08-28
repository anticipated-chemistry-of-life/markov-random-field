cmake_minimum_required(VERSION 3.30)

project(
  acol
  VERSION 0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Nothing here uses C++20 modules. Without this, CMake + Ninja run a module
# dependency scan on every .cpp and fail on toolchains that ship no
# clang-scan-deps (the conda-forge clang, for one).
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Build type" FORCE)
endif()

# =============================================================================
# Options
#
# Which sources of information to compile in? These are independent: any
# combination is valid except LOTUS=OFF and SIMPLE_DATA=OFF, which Types.h
# rejects with a static_assert (Y would not be informed by anything).
#
# The justfile drives these with single letters: l -> LOTUS, s -> SIMPLE_DATA, m
# -> USE_MS_DATA (see `just --list`).
# =============================================================================

option(LOTUS "Compile in the LOTUS data source" ON)
option(SIMPLE_DATA "Compile in the simple error model data source" ON)
option(USE_MS_DATA "Compile in the mass spec data source" OFF)
option(BUILD_TESTING "Configure the acol_unitTests target" ON)

# Which implementation backs the field and the internal state. The concepts in
# src/storages/storage_concepts.h are what make this a choice at all; the
# aliases in src/storages/storage_backend.h are where it is spent.
#
# Two binaries built from these same sources, differing only in this value, must
# produce byte-identical output from the same seed. `just parity` is what checks
# it (tests/backend_parity/); `dense` is the obviously-correct reference and
# `sparse` is the path to larger runs.
set(ACOL_STORAGE_BACKEND
    "sparse"
    CACHE STRING "Storage backend for the field and the internal state")
set_property(CACHE ACOL_STORAGE_BACKEND PROPERTY STRINGS dense sparse)

if(NOT ACOL_STORAGE_BACKEND MATCHES "^(dense|sparse)$")
  message(
    FATAL_ERROR
      "Unknown ACOL_STORAGE_BACKEND '${ACOL_STORAGE_BACKEND}'. The backends are 'dense' and 'sparse'."
  )
endif()

string(TOUPPER "${ACOL_STORAGE_BACKEND}" ACOL_STORAGE_BACKEND_MACRO)

if(NOT LOTUS AND NOT SIMPLE_DATA)
  message(
    FATAL_ERROR
      "No source of data was compiled in: Y could not be informed by anything. "
      "Configure with at least one of -DLOTUS=ON or -DSIMPLE_DATA=ON (or 'just build l')."
  )
endif()

message(
  STATUS
    "acol data sources: LOTUS=${LOTUS} SIMPLE_DATA=${SIMPLE_DATA} USE_MS_DATA=${USE_MS_DATA}"
)
message(STATUS "acol storage backend: ${ACOL_STORAGE_BACKEND}")

# =============================================================================
# Dependencies
# =============================================================================

include(FetchContent)

# SOURCE_DIR points into the source tree (not the build tree) so that the debug
# and release build directories share a single checkout.
#
# SYSTEM makes the fetched targets hand their include directories to us as
# -isystem. Without it the strict warning set below is applied to the headers of
# armadillo, fmt and friends -- which xmake never did, because it pulls those in
# as packages -- and -Wformat=2 -Werror fails the build inside armadillo_bits.
FetchContent_Declare(
  coretools
  GIT_REPOSITORY https://bitbucket.org/wegmannlab/coretools.git
  GIT_TAG master
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/coretools" SYSTEM)

FetchContent_Declare(
  stattools
  GIT_REPOSITORY https://bitbucket.org/wegmannlab/stattools.git
  GIT_TAG master
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/stattools" SYSTEM)

FetchContent_MakeAvailable(coretools stattools)

# Header-only helpers used by src/ntfy. Prefer an installed copy (the conda
# environment ships nlohmann_json), otherwise fetch them.
find_package(nlohmann_json 3.11 QUIET)
if(NOT nlohmann_json_FOUND)
  FetchContent_Declare(
    nlohmann_json
    URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
        SYSTEM)
  FetchContent_MakeAvailable(nlohmann_json)
endif()

find_package(stduuid QUIET)
if(NOT stduuid_FOUND)
  set(UUID_BUILD_TESTS
      OFF
      CACHE BOOL "" FORCE)
  set(UUID_SYSTEM_GENERATOR
      OFF
      CACHE BOOL "" FORCE)
  FetchContent_Declare(
    stduuid
    GIT_REPOSITORY https://github.com/mariusbancila/stduuid.git
    GIT_TAG v1.2.3
    SYSTEM)
  FetchContent_MakeAvailable(stduuid)
endif()

# armadillo, fast_float, fmt, zlib and OpenMP are deliberately not searched for
# here: coretools links them PUBLIC, so they reach us through the coretools
# target, exactly as they did through xmake's `{public = true}` packages.
# Looking armadillo up again is not merely redundant -- coretools builds its own
# vendored copy from coretools/arma, so find_package would put a second,
# unrelated armadillo on the link line while the headers come from the vendored
# one.

# src/ntfy/TNtfyNotifier.h builds httplib with CPPHTTPLIB_OPENSSL_SUPPORT.
find_package(OpenSSL REQUIRED)

# =============================================================================
# Sources
# =============================================================================

file(GLOB_RECURSE ACOL_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE ACOL_HEADERS CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
file(GLOB_RECURSE ACOL_TEST_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")

# -Wreturn-local-addr is GCC-only; clang spells it -Wreturn-stack-address and
# would reject the unknown name under -Werror.
set(ACOL_WARNINGS
    -Wall
    -Wextra
    -Werror
    -Wpedantic
    -Wuninitialized
    -Warray-bounds
    -Wnull-dereference
    -Wdouble-promotion
    -Wformat=2
    -Wundef
    $<$<CXX_COMPILER_ID:GNU>:-Wreturn-local-addr>)

# Everything src/**.cpp needs, regardless of which binary it ends up in. The
# unit tests recompile src/ rather than linking a shared object library because
# they add CHECK_INTERVALS, which changes the coretools headers the sources are
# compiled against.
function(acol_add_common_sources target)
  target_sources(${target} PRIVATE ${ACOL_SOURCES} ${ACOL_HEADERS})
  target_include_directories(${target}
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
  target_link_libraries(
    ${target} PRIVATE coretools stattools OpenSSL::SSL OpenSSL::Crypto
                      nlohmann_json::nlohmann_json stduuid)
  target_compile_features(${target} PRIVATE cxx_std_20)
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)

  target_compile_definitions(
    ${target} PRIVATE ACOL_STORAGE_${ACOL_STORAGE_BACKEND_MACRO})

  if(LOTUS)
    target_compile_definitions(${target} PRIVATE USE_LOTUS)
  endif()
  if(SIMPLE_DATA)
    target_compile_definitions(${target} PRIVATE USE_SIMPLE_ERROR_MODEL)
  endif()
  if(USE_MS_DATA)
    target_compile_definitions(${target} PRIVATE USE_MS_DATA)
  endif()

  # httplib reads the trust store from the macOS keychain when built against
  # OpenSSL.
  if(APPLE)
    target_link_libraries(${target} PRIVATE "-framework CoreFoundation"
                                            "-framework Security")
  endif()
endfunction()

# =============================================================================
# Main executable
# =============================================================================

add_executable(${PROJECT_NAME} main.cpp)
acol_add_common_sources(${PROJECT_NAME})
target_compile_definitions(${PROJECT_NAME} PRIVATE DEVTOOLS DEV_LOCATION)
target_compile_options(${PROJECT_NAME} PRIVATE ${ACOL_WARNINGS})

# =============================================================================
# Unit tests
#
# EXCLUDE_FROM_ALL: `just build` only builds acol, `just test` asks for this
# target explicitly.
# =============================================================================

if(BUILD_TESTING)
  enable_testing()

  set(gtest_force_shared_crt
      ON
      CACHE BOOL "" FORCE)
  set(gtest_disable_pthreads
      ON
      CACHE BOOL "" FORCE)
  set(INSTALL_GTEST
      OFF
      CACHE BOOL "" FORCE)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG main)
  FetchContent_MakeAvailable(googletest)

  add_executable("${PROJECT_NAME}_unitTests" EXCLUDE_FROM_ALL
                 ${ACOL_TEST_SOURCES})
  acol_add_common_sources("${PROJECT_NAME}_unitTests")
  target_include_directories("${PROJECT_NAME}_unitTests"
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
  target_link_libraries("${PROJECT_NAME}_unitTests" PRIVATE GTest::gtest_main
                                                            GTest::gmock_main)
  target_compile_definitions(
    "${PROJECT_NAME}_unitTests"
    PRIVATE CHECK_INTERVALS
            ACOL_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/test_data")
  target_compile_options("${PROJECT_NAME}_unitTests" PRIVATE -Wall -Wextra)

  add_test(NAME "${PROJECT_NAME}_unitTests" COMMAND "${PROJECT_NAME}_unitTests")
endif()
