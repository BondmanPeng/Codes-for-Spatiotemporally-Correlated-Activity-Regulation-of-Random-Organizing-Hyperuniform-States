cmake_minimum_required(VERSION 3.15)
project(modulate_IG CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

add_executable(randorg_binary_idealA main_binary_idealA.cpp)
add_executable(compute_msd_B compute_msd_B.cpp)
add_executable(compute_msd_A compute_msd_A.cpp)
add_executable(compute_activity_B compute_activity_B.cpp)

if(EXISTS "${CMAKE_SOURCE_DIR}/compute_sq_species.cpp")
  add_executable(compute_sq_species compute_sq_species.cpp)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(compute_sq_species PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

enable_testing()
add_executable(test_disk tests/test_disk.cpp)
target_include_directories(test_disk PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_options(test_disk PRIVATE -UNDEBUG)
add_test(NAME exact_disk_selfcheck COMMAND test_disk)
