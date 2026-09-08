cmake_minimum_required(VERSION 3.16)

project(pure_mce LANGUAGES CXX)

option(PURE_DISABLE_PRUNING_NORMALIZATION "Disable normalization" OFF)
option(PURE_DISABLE_PRUNING_SUBSUMPTION "Disable subsumption" OFF)
option(PURE_DISABLE_PRUNING_UNIT "Disable unit propagation" OFF)
option(PURE_DISABLE_PRUNING_USEFULNESS "Disable usefulness filtering" OFF)
option(PURE_DISABLE_PRUNING_ANTICHAIN "Disable antichain maintenance" OFF)
option(PURE_DISABLE_PRUNING_FAIL_FIRST "Disable fail-first pruning" OFF)
option(PURE_DISABLE_PRUNING_ZERO_COVERAGE "Disable zero-coverage pruning" OFF)

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
foreach(rule IN ITEMS NORMALIZATION SUBSUMPTION UNIT USEFULNESS ANTICHAIN
                      FAIL_FIRST ZERO_COVERAGE)
    if(PURE_DISABLE_PRUNING_${rule})
        target_compile_definitions(pure_core PUBLIC
            PURE_DISABLE_PRUNING_${rule}=1)
    endif()
endforeach()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pure_core PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(pure_mce main.cpp)
target_link_libraries(pure_mce PRIVATE pure_core)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pure_mce PRIVATE -Wall -Wextra -Wpedantic)
endif()
