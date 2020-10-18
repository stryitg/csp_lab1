project(cv-labs)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(lab ${SOURCES})
target_include_directories(lab PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/)
