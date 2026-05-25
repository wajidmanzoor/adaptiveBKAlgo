cmake_minimum_required(VERSION 3.10)

project(bk_algorithm LANGUAGES CXX)

option(REORDERSIB_PROFILING "Enable ReorderSib profiling timers" OFF)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags per build type
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Include directories
include_directories(inc)

# Sources
add_executable(bk_algorithm
    main.cpp
    src/common.cpp
    src/graph.cpp
    src/helpers.cpp
)

if(REORDERSIB_PROFILING)
    target_compile_definitions(bk_algorithm PRIVATE PROFILING=1)
endif()
