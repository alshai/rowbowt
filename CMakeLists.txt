cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project (rowbowt)
enable_testing()
set(CMAKE_CXX_FLAGS "--std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/sdsl-lite/include)
include_directories(${PROJECT_SOURCE_DIR}/pfbwt-f/include)

add_subdirectory(pfbwt-f)
add_subdirectory(${PROJECT_SOURCE_DIR}/tests)

message("Building in ${CMAKE_BUILD_TYPE} mode")

add_executable(rb_build src/rb_build.cpp)
add_executable(rb_align src/rb_align.cpp)
add_executable(rb_markers src/rb_markers.cpp)
add_executable(rb_locs src/rb_markers_tsa.cpp)
add_executable(build_midx src/build_midx.cpp)
TARGET_LINK_LIBRARIES(rb_align z)
TARGET_LINK_LIBRARIES(rb_markers z)
TARGET_LINK_LIBRARIES(rb_locs z)

install(TARGETS rb_build DESTINATION bin)
install(TARGETS rb_align DESTINATION bin)
install(TARGETS rb_markers DESTINATION bin)
install(TARGETS rb_locs DESTINATION bin)
install(TARGETS build_midx DESTINATION bin)
