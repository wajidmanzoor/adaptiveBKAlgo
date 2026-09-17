cmake_minimum_required(VERSION 3.16)

project(adaptive_bk_algorithm LANGUAGES CXX)

include(CheckIPOSupported)
set(PURE_HITSET_VARIANT "128" CACHE STRING
    "Seed-mask capacity: 64, 128, 512, 1024, or dynamic")
set_property(CACHE PURE_HITSET_VARIANT PROPERTY STRINGS
    64 128 512 1024 dynamic)

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
if(PURE_HITSET_VARIANT STREQUAL "dynamic")
    target_compile_definitions(reorder_core PUBLIC PURE_HITSET_DYNAMIC=1)
elseif(PURE_HITSET_VARIANT MATCHES "^(64|128|512|1024)$")
    target_compile_definitions(reorder_core PUBLIC
        PURE_HITSET_CAPACITY=${PURE_HITSET_VARIANT})
else()
    message(FATAL_ERROR
        "PURE_HITSET_VARIANT must be 64, 128, 512, 1024, or dynamic")
endif()

add_executable(adaptive_bk main.cpp)
target_link_libraries(adaptive_bk PRIVATE reorder_core)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(reorder_core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(adaptive_bk PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(reorder_core PRIVATE
        $<$<CONFIG:Release>:-march=native>)
    target_compile_options(adaptive_bk PRIVATE
        $<$<CONFIG:Release>:-march=native>)
endif()

check_ipo_supported(RESULT REORDER_IPO_SUPPORTED OUTPUT REORDER_IPO_ERROR)
if(REORDER_IPO_SUPPORTED)
    set_property(TARGET reorder_core adaptive_bk PROPERTY
        INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(ccr_exhaustive tests/ccr_exhaustive.cpp)
    target_link_libraries(ccr_exhaustive PRIVATE reorder_core)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(ccr_exhaustive PRIVATE
            -Wall -Wextra -Wpedantic)
        target_compile_options(ccr_exhaustive PRIVATE
            $<$<CONFIG:Release>:-march=native>)
    endif()
    if(REORDER_IPO_SUPPORTED)
        set_property(TARGET ccr_exhaustive PROPERTY
            INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
    add_test(
        NAME ccr_exhaustive
        COMMAND ccr_exhaustive
            "${CMAKE_CURRENT_BINARY_DIR}/ccr_exhaustive.graph")
    set_tests_properties(ccr_exhaustive PROPERTIES TIMEOUT 180)
endif()
