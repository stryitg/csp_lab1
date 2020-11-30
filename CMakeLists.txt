project(cv-labs)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(lab1 ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/lab1/main.cpp)
target_include_directories(lab1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/)

find_package(OpenCV COMPONENTS core highgui imgcodecs imgproc REQUIRED)

add_executable(lab2 ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/lab2/main.cpp)
target_link_directories(lab2 PRIVATE ${OpenCV_LIBRARY_DIR})
target_link_libraries(lab2 PRIVATE ${OpenCV_LIBRARIES})
target_include_directories(lab2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/ ${OpenCV_INCLUDE_DIRS})
