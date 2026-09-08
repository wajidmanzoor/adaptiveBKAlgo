cmake_minimum_required(VERSION 3.16)

project(adaptive_bk_algorithm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(reorder_core STATIC
    src/fast_plex3.cpp
    src/graph.cpp
    src/helpers.cpp
)
target_include_directories(reorder_core PUBLIC inc)

add_executable(adaptive_bk main.cpp)
target_link_libraries(adaptive_bk PRIVATE reorder_core)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(reorder_core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(adaptive_bk PRIVATE -Wall -Wextra -Wpedantic)
endif()
