cmake_minimum_required(VERSION 3.10)

project(bk_algorithm LANGUAGES CXX)

option(REORDERSIB_PROFILING "Enable ReorderSib profiling timers" OFF)
option(REORDERSIB_SOLVER_TRACE "Enable per-sibling solver trace logging" OFF)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# HBBMC is built with OpenMP when available and also uses AVX on x86. Mirror
# those opportunities here so our comparison is not handicapped by weaker
# compile flags.
include(CheckCXXCompilerFlag)
find_package(OpenMP)
check_cxx_compiler_flag("-mavx" REORDERSIB_COMPILER_SUPPORTS_MAVX)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" REORDERSIB_ARCH_LOWER)
set(REORDERSIB_X86 FALSE)
if(REORDERSIB_ARCH_LOWER MATCHES "^(x86_64|amd64|i[3-6]86)$")
    set(REORDERSIB_X86 TRUE)
endif()

# Compiler flags per build type
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Include directories
include_directories(inc)

# Sources
add_executable(bk_algorithm
    main.cpp
    src/fast_plex3.cpp
    src/graph.cpp
    src/helpers.cpp
)

if(REORDERSIB_PROFILING)
    target_compile_definitions(bk_algorithm PRIVATE PROFILING=1)
endif()

if(REORDERSIB_SOLVER_TRACE)
    target_compile_definitions(bk_algorithm PRIVATE SOLVER_TRACE=1)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(bk_algorithm PRIVATE OpenMP::OpenMP_CXX)
    if(REORDERSIB_COMPILER_SUPPORTS_MAVX AND REORDERSIB_X86)
        target_compile_options(bk_algorithm PRIVATE $<$<CONFIG:Release>:-mavx>)
    endif()
endif()
