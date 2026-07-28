cmake_minimum_required(VERSION 3.30)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Project
project(acol LANGUAGES CXX)

# Executable
add_executable(${PROJECT_NAME} main.cpp)

file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS src/* tests/integrationtests/*)

# Download coretools and stattools
include(FetchContent)
FetchContent_Declare(coretools
 GIT_REPOSITORY https://bitbucket.org/wegmannlab/coretools.git
 GIT_TAG master
 SOURCE_DIR "${CMAKE_SOURCE_DIR}/coretools_cmake"
)

FetchContent_Declare(stattools
 GIT_REPOSITORY https://bitbucket.org/wegmannlab/stattools.git
 GIT_TAG master
 SOURCE_DIR "${CMAKE_SOURCE_DIR}/stattools_cmake"
)

FetchContent_MakeAvailable(coretools stattools)

find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif ()


# Which sources of information to compile in? These are independent: any combination is valid
# except all-off, which Types.h rejects with a static_assert. Defaults mirror xmake.lua.
option(LOTUS "Compile in the LOTUS data source" ON)
if (LOTUS)
  message("Using LOTUS data (run 'cmake -DLOTUS=OFF ..' to change)")
  # add compile definition to source files (everything behind -D will be replaced with bool)
  add_definitions(-DUSE_LOTUS)
else ()
  message("Not using LOTUS data (run 'cmake -DLOTUS=ON ..' to change)")
endif ()

option(SIMPLE_DATA "Compile in the simple error model data source" OFF)
if (SIMPLE_DATA)
  message("Using simple error model data (run 'cmake -DSIMPLE_DATA=OFF ..' to change)")
  add_definitions(-DUSE_SIMPLE_ERROR_MODEL)
else ()
  message("Not using simple error model data (run 'cmake -DSIMPLE_DATA=ON ..' to change)")
endif ()


target_sources(${PROJECT_NAME} PRIVATE ${CORE_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} PRIVATE coretools stattools ${ZLIB_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
target_compile_definitions(${PROJECT_NAME} PRIVATE DEVTOOLS DEV_LOCATION)

# Test
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_executable("${PROJECT_NAME}_unitTests" EXCLUDE_FROM_ALL)

  file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS tests/*)

  include(FetchContent)
  FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG main
    )
  FetchContent_MakeAvailable(googletest)

  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  option(gtest_disable_pthreads "Disable uses of pthreads in gtest." ON)

  target_sources("${PROJECT_NAME}_unitTests" PRIVATE ${CORE_SOURCES} ${TEST_SOURCES})
  target_include_directories("${PROJECT_NAME}_unitTests" PRIVATE tests src)
  target_link_libraries("${PROJECT_NAME}_unitTests" PRIVATE coretools stattools gtest_main gmock_main)
  target_compile_definitions("${PROJECT_NAME}_unitTests" PRIVATE CHECK_INTERVALS)
  target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
  set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
  target_compile_options("${PROJECT_NAME}_unitTests" PRIVATE -Wall -Wextra)
endif()
