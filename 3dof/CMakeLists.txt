cmake_minimum_required(VERSION 3.14)
project(3dof_kalman_filter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(3dof_kalman_filter
  3dof_kalman_filter.cc
)
target_include_directories(3dof_kalman_filter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(3dof_kalman_filter_main
  3dof_kalman_filter_main.cc
)
target_link_libraries(3dof_kalman_filter_main PRIVATE 3dof_kalman_filter)
