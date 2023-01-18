cmake_minimum_required(VERSION 3.10)
project(example)
include(FetchContent)

function(enable_ninja_output_coloring target)
    if(UNIX AND CMAKE_GENERATOR STREQUAL "Ninja")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            target_compile_options(${target} PRIVATE -fcolor-diagnostics)
            target_compile_options(${target} PRIVATE -fcolor-diagnostics)
        endif()
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            target_compile_options(${target} PRIVATE -fdiagnostics-color)
            target_compile_options(${target} PRIVATE -fdiagnostics-color)
        endif()
    endif()
endfunction(enable_ninja_output_coloring)

function(compile_target target)
    enable_ninja_output_coloring(${target})
    target_include_directories(${target} PUBLIC Eigen3::Eigen OsqpEigen::OsqpEigen ${osqp_INCLUDE_DIR})
    target_link_libraries(${target} Eigen3::Eigen OsqpEigen::OsqpEigen)
endfunction(compile_target)

# FetchContent_Declare(
#     osqp-cpp
#     GIT_REPOSITORY https://github.com/google/osqp-cpp.git
#     GIT_TAG        origin/master
#   )
#   FetchContent_MakeAvailable(osqp-cpp)
#   message("${osqp_INCLUDE_DIR}")

set(EIGEN_NO_DEBUG ON)
find_package(Eigen3 REQUIRED)
# find_package(GiNaC REQUIRED)
find_package(OsqpEigen REQUIRED)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/show.plt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(qp_train qp_train.cpp)
# add_executable(osqp_cpp osqp_cpp.cpp)
add_executable(osqp_sample osqp_sample.cpp)
add_executable(sample_revise sample_revise.cpp)
add_executable(noconstraint noconstraint.cpp)

compile_target(sample_revise)
compile_target(noconstraint)

enable_ninja_output_coloring(qp_train)
enable_ninja_output_coloring(osqp_sample)
target_include_directories(qp_train PUBLIC Eigen3::Eigen OsqpEigen::OsqpEigen)
target_link_libraries(qp_train Eigen3::Eigen OsqpEigen::OsqpEigen)
target_link_libraries(osqp_sample Eigen3::Eigen OsqpEigen::OsqpEigen)
install(TARGETS qp_train DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
