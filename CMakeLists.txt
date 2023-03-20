cmake_minimum_required(VERSION 3.14...3.26)

get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release default")
endif()

project(stlbfgs LANGUAGES CXX)

enable_testing()

option(STLBFGS_UNIT_TESTS "STLBFGS_UNIT_TESTS" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)

find_package(OpenMP COMPONENTS CXX)

if(MSVC)
  add_compile_options(/W4)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|Intel")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(RELATIVE_BIN_DIR bin/)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)

add_library(stlbfgs linesearch.h stlbfgs.h linesearch.cpp stlbfgs.cpp)
target_link_libraries(stlbfgs PRIVATE $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_include_directories(stlbfgs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET stlbfgs PROPERTY SOVERSION 1.0)

if (STLBFGS_UNIT_TESTS)
  add_subdirectory(tests)
endif()

add_executable(stlbfgs-helloworld helloworld.cpp)
target_link_libraries(stlbfgs-helloworld ${CMAKE_DL_LIBS} stlbfgs $<$<NOT:$<BOOL:${WIN32}>>:m>)

file(GENERATE OUTPUT .gitignore CONTENT "*")
