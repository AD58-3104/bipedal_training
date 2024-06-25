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

include(FetchContent) # If not included already

FetchContent_Declare(sciplot
  GIT_REPOSITORY https://github.com/sciplot/sciplot.git
  GIT_TAG master)

FetchContent_GetProperties(sciplot)
if(NOT sciplot_content_POPULATED)
  FetchContent_Populate(sciplot)
endif()
file(GLOB_RECURSE SCPLOT_SOURCES ${sciplot_SOURCE_DIR}/sciplot/*)


function(compile_target target)
    enable_ninja_output_coloring(${target})
    target_include_directories(${target} PUBLIC Eigen3::Eigen ${sciplot_SOURCE_DIR})
    target_link_libraries(${target} Eigen3::Eigen)
    target_sources(${target} PRIVATE ${SCPLOT_SOURCES})
endfunction(compile_target)


set(EIGEN_NO_DEBUG ON)
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/show.plt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(main main.cpp)
compile_target(main)