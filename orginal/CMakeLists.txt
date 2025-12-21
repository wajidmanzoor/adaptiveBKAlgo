cmake_minimum_required(VERSION 3.10)

project(bk_algorithm LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
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
