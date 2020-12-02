cmake_minimum_required(VERSION 3.16)
project(opna_2)
set(CMAKE_CXX_STANDARD 17)
set(Boost_USE_STATIC_LIBS ON)

include_directories(libs/catch)
include_directories(libs/prettyprint)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

set(FORT_ENABLE_TESTING OFF CACHE INTERNAL "")
add_subdirectory(libs/libfort)

add_executable(opna_2 main.cpp)
target_link_libraries(opna_2 PRIVATE ${Boost_LIBRARIES} fort)

add_executable(opna_2_tests main.cpp)
target_link_libraries(opna_2_tests PRIVATE ${Boost_LIBRARIES} fort)
target_compile_definitions(opna_2_tests PRIVATE TESTING_ENABLED)
