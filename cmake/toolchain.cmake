# Pick the compiler the machine can build with, and the prefix it looks for libraries in.
#
# The presets name this file, so every configure goes through it: the pixi tasks, `pixi run
# parity`, and a bare `cmake --preset debug`.
#
# The environment's compiler is the first choice, because the binary links the environment's
# libraries. CC and CXX win over it when they are set.
#
# On macOS that first choice needs a probe. The conda toolchain is pinned to the SDKs it knew
# about when it was built. A macOS SDK newer than that -- the weeks after an OS or Xcode upgrade,
# before conda-forge catches up -- breaks it two ways. Its linker fails to link even an empty
# program, sometimes only for C. Its libc++ headers fail to compile `<random>` at C++20, where
# they meet a newer `math.h` and stop at an undeclared INFINITY. Both are compiled here, and a
# failure falls back to the Xcode toolchain, which always matches its own SDK.

# Look in the environment before the machine. cmake searches CMAKE_PREFIX_PATH ahead of the
# platform's own prefixes, so a library the environment ships wins over one the machine happens
# to carry. Without this, find_package(OpenSSL) takes the distribution's libcrypto and hands it
# to the environment's linker, whose sysroot is older than the glibc that library was built
# against:
#
#     x86_64-conda-linux-gnu-ld: /usr/lib/x86_64-linux-gnu/libcrypto.so:
#     undefined reference to `dlopen@GLIBC_2.34'
#
# This runs before the early return below, because a reconfigure re-reads the file with the
# compilers already chosen and still has to search the same prefixes.
if(DEFINED ENV{CONDA_PREFIX})
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

# CMake re-reads this file for every language it detects and for every try_compile. The compilers
# are on the command line by then, so the work below happens once per build directory.
if(DEFINED CMAKE_C_COMPILER AND DEFINED CMAKE_CXX_COMPILER)
  return()
endif()

execute_process(
  COMMAND uname -s
  OUTPUT_VARIABLE acol_system
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

if(DEFINED ENV{CC} AND DEFINED ENV{CXX})
  set(acol_cc "$ENV{CC}")
  set(acol_cxx "$ENV{CXX}")
elseif(DEFINED ENV{CONDA_PREFIX})
  if(acol_system STREQUAL "Darwin")
    set(acol_cc "$ENV{CONDA_PREFIX}/bin/clang")
    set(acol_cxx "$ENV{CONDA_PREFIX}/bin/clang++")
  else()
    set(acol_cc "$ENV{CONDA_PREFIX}/bin/gcc")
    set(acol_cxx "$ENV{CONDA_PREFIX}/bin/g++")
  endif()
else()
  # No environment to take a compiler from: leave the choice to cmake.
  return()
endif()

if(NOT EXISTS "${acol_cc}" OR NOT EXISTS "${acol_cxx}")
  return()
endif()

set(CMAKE_C_COMPILER "${acol_cc}")
set(CMAKE_CXX_COMPILER "${acol_cxx}")

if(NOT acol_system STREQUAL "Darwin")
  return()
endif()

find_program(acol_xcrun xcrun)
if(NOT acol_xcrun)
  return()
endif()

set(acol_probe_dir "${CMAKE_BINARY_DIR}/toolchain_probe")
file(WRITE "${acol_probe_dir}/probe.c" "int main(){return 0;}\n")
file(WRITE "${acol_probe_dir}/probe.cpp" "#include <random>\nint main(){return 0;}\n")
execute_process(
  COMMAND "${acol_cc}" probe.c -o probe_c
  WORKING_DIRECTORY "${acol_probe_dir}"
  RESULT_VARIABLE acol_c_probe
  OUTPUT_QUIET
  ERROR_QUIET
)
execute_process(
  COMMAND "${acol_cxx}" -std=c++20 probe.cpp -o probe_cxx
  WORKING_DIRECTORY "${acol_probe_dir}"
  RESULT_VARIABLE acol_cxx_probe
  OUTPUT_QUIET
  ERROR_QUIET
)
file(REMOVE_RECURSE "${acol_probe_dir}")

if(acol_c_probe EQUAL 0 AND acol_cxx_probe EQUAL 0)
  return()
endif()

execute_process(
  COMMAND "${acol_xcrun}" -f clang
  OUTPUT_VARIABLE acol_xcode_cc
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE acol_xcode_found
  ERROR_QUIET
)
if(NOT acol_xcode_found EQUAL 0)
  return()
endif()
execute_process(
  COMMAND "${acol_xcrun}" -f clang++
  OUTPUT_VARIABLE acol_xcode_cxx
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
execute_process(
  COMMAND "${acol_xcrun}" --show-sdk-path
  OUTPUT_VARIABLE acol_xcode_sdk
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

set(CMAKE_C_COMPILER "${acol_xcode_cc}")
set(CMAKE_CXX_COMPILER "${acol_xcode_cxx}")
set(CMAKE_OSX_SYSROOT "${acol_xcode_sdk}")
message(STATUS "acol: the environment's clang cannot build against this SDK; using ${acol_xcode_cxx}")
