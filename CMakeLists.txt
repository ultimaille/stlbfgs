cmake_minimum_required(VERSION 3.19)
project(stlbfgs LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES linesearch.cpp stlbfgs.cpp linesearch.h stlbfgs.h)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(stlbfgs 
  ${SOURCES} ${HEADERS})
target_include_directories(stlbfgs PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(stlbfgs OpenMP::OpenMP_CXX)
endif()
