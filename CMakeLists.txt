cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project (rowbowt)

find_library(LIB_HTS hts)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/9c047902ac8843cd6f972ae9ff472f2baf790879.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


set(CMAKE_CXX_FLAGS "--std=c++17 -D__STDC_FORMAT_MACROS")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/sdsl-lite/include)
include_directories(${PROJECT_SOURCE_DIR}/pfbwt-f/include)
include_directories(${PROJECT_SOURCE_DIR}/faster-minuter/include)

add_subdirectory(pfbwt-f)
# add_subdirectory(${PROJECT_SOURCE_DIR}/tests)

message("Building in ${CMAKE_BUILD_TYPE} mode")

enable_testing()
add_executable(rb_build src/rb_build.cpp)
add_executable(rb_align src/rb_align.cpp)
add_executable(rb_markers src/rb_markers.cpp)
add_executable(rb_locs src/rb_markers_tsa.cpp)
add_executable(build_midx src/build_midx.cpp)
add_executable(rb_tests tests/rb_tests.cpp)
target_link_libraries(rb_align z)
target_link_libraries(rb_markers z pthread)
target_link_libraries(rb_locs z)
target_link_libraries(rb_tests z gtest_main)

include(GoogleTest)
gtest_discover_tests(rb_tests)
include(CPack)

install(TARGETS rb_build DESTINATION bin)
install(TARGETS rb_align DESTINATION bin)
install(TARGETS rb_markers DESTINATION bin)
install(TARGETS rb_locs DESTINATION bin)
install(TARGETS build_midx DESTINATION bin)