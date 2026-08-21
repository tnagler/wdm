cmake_minimum_required(VERSION 3.13)

project(wdm_consumer LANGUAGES CXX)

find_package(wdm REQUIRED)

add_executable(wdm_consumer main.cpp)
target_link_libraries(wdm_consumer PRIVATE wdm)
