cmake_minimum_required(VERSION 3.16)

project(adaptive_bk_algorithm LANGUAGES CXX)

include(CheckIPOSupported)
option(PURE_DISABLE_PRUNING_NORMALIZATION "Disable normalization" OFF)
option(PURE_DISABLE_PRUNING_SUBSUMPTION "Disable subsumption" OFF)
option(PURE_DISABLE_PRUNING_UNIT "Disable unit propagation" OFF)
option(PURE_DISABLE_PRUNING_USEFULNESS "Disable usefulness filtering" OFF)
option(PURE_DISABLE_PRUNING_ANTICHAIN "Disable antichain maintenance" OFF)
option(PURE_DISABLE_PRUNING_FAIL_FIRST "Disable fail-first pruning" OFF)
option(PURE_DISABLE_PRUNING_ZERO_COVERAGE "Disable zero-coverage pruning" OFF)

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
foreach(rule IN ITEMS NORMALIZATION SUBSUMPTION UNIT USEFULNESS ANTICHAIN
                      FAIL_FIRST ZERO_COVERAGE)
    if(PURE_DISABLE_PRUNING_${rule})
        target_compile_definitions(reorder_core PUBLIC
            PURE_PRUNING_${rule}=0)
    else()
        target_compile_definitions(reorder_core PUBLIC
            PURE_PRUNING_${rule}=1)
    endif()
endforeach()

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
