cmake_minimum_required(VERSION 3.10)
project(eigen_train)

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
    target_include_directories(${target} PUBLIC Eigen3::Eigen)
    target_link_libraries(${target} Eigen3::Eigen)
endfunction(compile_target)


set(EIGEN_NO_DEBUG ON)
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/show.plt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(train main.cpp)
compile_target(train)