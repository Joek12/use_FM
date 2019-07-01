cmake_minimum_required(VERSION 3.14)
project(use_FM)

set(CMAKE_CXX_STANDARD 14)

add_executable(use_FM src/main.cpp src/fastaReader.cpp src/fastaReader.h)