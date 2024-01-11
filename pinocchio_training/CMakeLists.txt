cmake_minimum_required(VERSION 3.20)
project(pinocchio_training)
set(CMAKE_PREFIX_PATH "path to source or build dir")


add_executable(compile compile.cpp)

list(APPEND compile_targes compile)

foreach(target IN LISTS compile_targes)
    target_compile_options(${target} PRIVATE -Wall -Wextra -O2 -march=native -mtune=native)
    target_compile_features(${target} PRIVATE cxx_std_20)
endforeach()
