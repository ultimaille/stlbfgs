Include(FetchContent)

set(FETCHCONTENT_QUIET false)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.3.1
  GIT_SHALLOW true
  TLS_VERIFY ON
  INACTIVITY_TIMEOUT 60
  )

FetchContent_MakeAvailable(Catch2)

FILE(GLOB SRCTEST *.cpp)

add_executable(stlbfgs-test-all ${SRCTEST})

target_link_libraries(stlbfgs-test-all PRIVATE ${CMAKE_DL_LIBS} stlbfgs Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(stlbfgs-test-all)
