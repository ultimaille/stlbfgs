add_executable(helloworld helloworld.cpp)
target_link_libraries(helloworld ${CMAKE_DL_LIBS} stlbfgs)

IF (NOT WIN32)
  target_link_libraries(helloworld m)
ENDIF()
