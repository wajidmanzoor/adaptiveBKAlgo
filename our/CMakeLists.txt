cmake_minimum_required(VERSION 3.10)

project(bk_algorithm LANGUAGES CXX)

option(REORDERSIB_PROFILING "Enable ReorderSib profiling timers" OFF)
set(PURE_HITSET_VARIANT "128" CACHE STRING
    "Pure optimized hitting-set coverage: 128, 256, or dynamic")
set_property(CACHE PURE_HITSET_VARIANT PROPERTY STRINGS 128 256 dynamic)

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

set(BK_CORE_SOURCES
    src/common.cpp
    src/fast_list_bk.cpp
    src/fast_local_bitset.cpp
    src/fast_plex3.cpp
    src/graph.cpp
    src/helpers.cpp
    src/rmce_reduction.cpp
)

add_library(bk_core STATIC ${BK_CORE_SOURCES})
if(PURE_HITSET_VARIANT STREQUAL "dynamic")
    target_compile_definitions(bk_core PUBLIC PURE_HITSET_DYNAMIC=1)
elseif(PURE_HITSET_VARIANT STREQUAL "128")
    target_compile_definitions(bk_core PUBLIC PURE_HITSET_CAPACITY=128)
elseif(PURE_HITSET_VARIANT STREQUAL "256")
    target_compile_definitions(bk_core PUBLIC PURE_HITSET_CAPACITY=256)
else()
    message(FATAL_ERROR
        "PURE_HITSET_VARIANT must be one of: 128, 256, dynamic")
endif()
add_executable(bk_algorithm main.cpp)
target_link_libraries(bk_algorithm PRIVATE bk_core)

if(REORDERSIB_PROFILING)
    target_compile_definitions(bk_core PRIVATE PROFILING=1)
endif()
