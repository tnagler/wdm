cmake_minimum_required(VERSION 3.13)

project(wdm_example LANGUAGES CXX)

find_package(wdm CONFIG REQUIRED)

add_executable(wdm_example main.cpp)
target_link_libraries(wdm_example PRIVATE wdm)
