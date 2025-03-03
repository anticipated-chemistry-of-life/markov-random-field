cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

# Project
project(acol LANGUAGES CXX)

# Executable
add_executable(${PROJECT_NAME} main.cpp)

file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS src/* tests/integrationtests/*)

# Download coretools and stattools
include(FetchContent)
FetchContent_Declare(coretools
        GIT_REPOSITORY https://bitbucket.org/wegmannlab/coretools.git
        GIT_TAG madleina/interface_stattools
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/coretools"
)

FetchContent_Declare(stattools
        GIT_REPOSITORY https://bitbucket.org/wegmannlab/stattools.git
        GIT_TAG madleina/interface
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/stattools"
)

FetchContent_Declare(genometools
        GIT_REPOSITORY https://bitbucket.org/wegmannlab/genometools.git
        GIT_TAG master
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/genometools"
)
FetchContent_MakeAvailable(coretools stattools genometools)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif ()

target_sources(${PROJECT_NAME} PRIVATE ${CORE_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} PRIVATE coretools stattools genometools ${ZLIB_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
target_compile_definitions(${PROJECT_NAME} PRIVATE DEVTOOLS DEV_LOCATION)

# Test
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
target_link_libraries("${PROJECT_NAME}_unitTests" PRIVATE coretools stattools genometools gtest_main gmock_main)
target_compile_definitions("${PROJECT_NAME}_unitTests" PRIVATE CHECK_INTERVALS)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options("${PROJECT_NAME}_unitTests" PRIVATE -Wall -Wextra)
