cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-stlbfgs test-stlbfgs.cpp)
target_link_libraries(test-stlbfgs PRIVATE ${CMAKE_DL_LIBS} stlbfgs Catch2::Catch2)

add_executable(test-linesearch test-linesearch.cpp)
target_link_libraries(test-linesearch PRIVATE ${CMAKE_DL_LIBS} stlbfgs Catch2::Catch2)

