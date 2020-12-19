project(cv-labs)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

find_package(OpenCV COMPONENTS core highgui imgcodecs imgproc REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

macro(create_executable name)
    add_executable(${name} ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/${name}/main.cpp)
    target_link_directories(${name} PRIVATE ${OpenCV_LIBRARY_DIR})
    target_link_libraries(${name} PRIVATE Eigen3::Eigen ${OpenCV_LIBRARIES})
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/ ${OpenCV_INCLUDE_DIRS})
endmacro()

create_executable(lab1)
create_executable(lab2)
create_executable(lab3)