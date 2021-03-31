add_executable(untangle untangle.cpp)
target_link_libraries(untangle ${CMAKE_DL_LIBS} stlbfgs)

IF (NOT WIN32)
  target_link_libraries(untangle m)
ENDIF()
