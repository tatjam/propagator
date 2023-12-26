cmake_minimum_required(VERSION 3.27)
set(CMAKE_CXX_STANDARD 17)

project(propagador)

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(propagador ${SOURCES})
include_directories(src)
