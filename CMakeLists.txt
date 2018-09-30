cmake_minimum_required(VERSION 3.10)
project(Software_ORAM_new)

set(CMAKE_CXX_STANDARD 11)

add_executable(Software_ORAM_new main.cpp ORAM.cpp ORAM.h Util.cpp Util.h)