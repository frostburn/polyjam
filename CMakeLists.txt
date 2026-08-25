cmake_minimum_required(VERSION 3.20)
project(polyjam LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(polyjam
  src/main.cpp
  src/geometry.cpp
  src/search.cpp
)

target_include_directories(polyjam PRIVATE src)

if (MSVC)
  target_compile_options(polyjam PRIVATE /W4 /permissive-)
else()
  target_compile_options(polyjam PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

find_package(Threads REQUIRED)
target_link_libraries(polyjam PRIVATE Threads::Threads)

include(CTest)
if (BUILD_TESTING)
  add_executable(polyjam_smoke tests/smoke.cpp src/geometry.cpp)
  target_include_directories(polyjam_smoke PRIVATE src)
  if (MSVC)
    target_compile_options(polyjam_smoke PRIVATE /W4 /permissive-)
  else()
    target_compile_options(polyjam_smoke PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME polyjam_geometry_smoke COMMAND polyjam_smoke)

  add_executable(polyjam_search_test tests/search.cpp src/geometry.cpp src/search.cpp)
  target_include_directories(polyjam_search_test PRIVATE src)
  target_link_libraries(polyjam_search_test PRIVATE Threads::Threads)
  if (MSVC)
    target_compile_options(polyjam_search_test PRIVATE /W4 /permissive-)
  else()
    target_compile_options(polyjam_search_test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME polyjam_lattice_search COMMAND polyjam_search_test)
endif()
