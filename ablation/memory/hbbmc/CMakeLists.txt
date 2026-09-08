cmake_minimum_required(VERSION 3.16)

project(hbbmc_faithful_reimplementation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(hbbmc_faithful_lib
    src/graph.cpp
    src/truss_order.cpp
    src/early_termination.cpp
    src/enumerator.cpp
    src/reduction.cpp
)

target_include_directories(hbbmc_faithful_lib PUBLIC inc)
target_compile_options(hbbmc_faithful_lib PRIVATE -Wall -Wextra -Wpedantic)

add_executable(hbbmc_faithful src/main.cpp)
target_link_libraries(hbbmc_faithful PRIVATE hbbmc_faithful_lib)
target_compile_options(hbbmc_faithful PRIVATE -Wall -Wextra -Wpedantic)

include(CTest)
if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_exhaustive.cpp")
    add_executable(hbbmc_faithful_tests tests/test_exhaustive.cpp)
    target_link_libraries(hbbmc_faithful_tests PRIVATE hbbmc_faithful_lib)
    target_compile_options(hbbmc_faithful_tests PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_definitions(hbbmc_faithful_tests PRIVATE
        HBBMC_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME exhaustive_differential COMMAND hbbmc_faithful_tests)
endif()
