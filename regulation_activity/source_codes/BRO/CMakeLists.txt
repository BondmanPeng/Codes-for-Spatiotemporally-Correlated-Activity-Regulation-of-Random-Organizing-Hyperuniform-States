cmake_minimum_required(VERSION 3.16)
project(modulated_BRO CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

add_executable(randorg_binary_randA main_binary_randA.cpp)

if(EXISTS "${CMAKE_SOURCE_DIR}/compute_sq_species.cpp")
  add_executable(compute_sq_species compute_sq_species.cpp)
  find_package(OpenMP REQUIRED)
  target_link_libraries(compute_sq_species PRIVATE OpenMP::OpenMP_CXX)
endif()

enable_testing()
add_executable(test_binary tests/test_binary.cpp)
target_include_directories(test_binary PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_options(test_binary PRIVATE -UNDEBUG)
add_test(NAME binary_selfcheck COMMAND test_binary)
