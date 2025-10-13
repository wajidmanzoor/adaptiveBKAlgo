cmake_minimum_required(VERSION 3.10)

# Project name
project(bk_algorithm)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2")

# Include directories
include_directories(inc)

# Source files
set(SOURCES
    main.cpp
    src/common.cpp
    src/graph.cpp
    src/helpers.cpp
)

# Create executable
add_executable(bk_algorithm ${SOURCES})