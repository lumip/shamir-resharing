cmake_minimum_required(VERSION 3.0.2)
project(shamir-resharing CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
find_package(GMP REQUIRED)

include_directories(${GMP_INCLUDES})

add_executable(resharing main.cpp)
target_link_libraries(resharing ${GMP_LIBRARIES})

