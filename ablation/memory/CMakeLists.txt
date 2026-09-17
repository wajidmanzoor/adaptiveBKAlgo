cmake_minimum_required(VERSION 3.16)

project(final_memory_sampler LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(rss_sampler rss_sampler.cpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(rss_sampler PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(tomita_stream_adapter tomita_stream_adapter.cpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(tomita_stream_adapter PRIVATE
        -O3 -Wall -Wextra -Wpedantic)
endif()

set(TOMITA_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../compare/Tomita/src")
add_executable(tomita_retained
    "${TOMITA_SOURCE_ROOT}/LinkedList.c"
    "${TOMITA_SOURCE_ROOT}/MemoryManager.c"
    "${TOMITA_SOURCE_ROOT}/adjlist_algorithm.c"
    "${TOMITA_SOURCE_ROOT}/misc.c"
    "${TOMITA_SOURCE_ROOT}/adjlist.c"
)
target_include_directories(tomita_retained PRIVATE "${TOMITA_SOURCE_ROOT}")
target_compile_definitions(tomita_retained PRIVATE
    RETURN_CLIQUES_ONE_BY_ONE
    MINIMUM_CLIQUE_SIZE=3
    SMARTLENGTH
    ALLOW_ALLOC_ZERO_BYTES
)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(tomita_retained PRIVATE -O3 -fgnu89-inline)
endif()

set(CCRMCE_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../compare/CCRMCE/src")
add_executable(ccrmce_retained
    "${CCRMCE_SOURCE_ROOT}/graph/graph.cpp"
    "${CCRMCE_SOURCE_ROOT}/MClique/mcecrrv3.cpp"
    "${CCRMCE_SOURCE_ROOT}/run/maximalCliqueLargeClique.cpp"
)
target_include_directories(ccrmce_retained PRIVATE
    "${CCRMCE_SOURCE_ROOT}/graph"
    "${CCRMCE_SOURCE_ROOT}/MClique"
    "${CCRMCE_SOURCE_ROOT}/tools"
)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(ccrmce_retained PRIVATE
        -O3 -march=native -Wall -Wextra)
endif()
