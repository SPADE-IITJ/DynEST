cmake_minimum_required(VERSION 3.16)

project(DynEST LANGUAGES CXX)

option(DYNEST_NATIVE "Optimize binaries for the build host CPU" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP REQUIRED)

function(dynest_target target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)

    if(MSVC)
        target_compile_options(${target} PRIVATE /O2 /W4)
    else()
        target_compile_options(${target} PRIVATE -O3 -Wall -Wextra -Wpedantic)
        if(DYNEST_NATIVE)
            target_compile_options(${target} PRIVATE -march=native)
        endif()
    endif()
endfunction()

dynest_target(dynest src/d_escent.cpp)
dynest_target(dynest_benchmark src/benchmark.cpp)
dynest_target(static_eist src/static_eist_construction.cpp)
dynest_target(test_st src/test_st.cpp)
