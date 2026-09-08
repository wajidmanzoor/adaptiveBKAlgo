cmake_minimum_required(VERSION 3.16)

project(final_memory_sampler LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(rss_sampler rss_sampler.cpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(rss_sampler PRIVATE -Wall -Wextra -Wpedantic)
endif()
