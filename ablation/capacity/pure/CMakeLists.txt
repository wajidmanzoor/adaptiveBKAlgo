cmake_minimum_required(VERSION 3.16)

project(pure_mce LANGUAGES CXX)

set(PURE_HITSET_VARIANT "128" CACHE STRING
    "Seed-mask capacity: 64, 128, 512, 1024, or dynamic")
set_property(CACHE PURE_HITSET_VARIANT PROPERTY STRINGS
    64 128 512 1024 dynamic)

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

set(PURE_CORE_SOURCES
    src/fast_plex3.cpp
    src/graph.cpp
    src/helpers.cpp
)

add_library(pure_core STATIC ${PURE_CORE_SOURCES})
target_include_directories(pure_core PUBLIC inc)
if(PURE_HITSET_VARIANT STREQUAL "dynamic")
    target_compile_definitions(pure_core PUBLIC PURE_HITSET_DYNAMIC=1)
elseif(PURE_HITSET_VARIANT MATCHES "^(64|128|512|1024)$")
    target_compile_definitions(pure_core PUBLIC
        PURE_HITSET_CAPACITY=${PURE_HITSET_VARIANT})
else()
    message(FATAL_ERROR
        "PURE_HITSET_VARIANT must be 64, 128, 512, 1024, or dynamic")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pure_core PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(pure_mce main.cpp)
target_link_libraries(pure_mce PRIVATE pure_core)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pure_mce PRIVATE -Wall -Wextra -Wpedantic)
endif()
