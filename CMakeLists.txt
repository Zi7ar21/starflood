cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(src)

project(starflood VERSION 1.0)

add_executable(starflood src/image.cpp src/main.cpp)

target_include_directories(starflood PUBLIC src)
