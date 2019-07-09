cmake_minimum_required(VERSION 3.14)
project(use_FM)

set(CMAKE_CXX_STANDARD 14)

include_directories(/Users/josephkang/Downloads/boost_1_68_0/boost/python.hpp)

add_executable(use_FM src/main.cpp src/fastaReader.cpp src/fastaReader.h src/spammodule.cpp src/spammodule.h src/FMCommands.h src/pbar.h src/pbar.cpp)