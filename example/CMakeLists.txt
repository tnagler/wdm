cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project (Example)

# Setting default folders
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# C++ compile flags
if (NOT WIN32)
  set(CMAKE_CXX_FLAGS "-std=gnu++11 -Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type -O2 -DNDEBUG")
endif()

# Find wdm package
find_package(wdm                     REQUIRED)

# Include header files
include_directories(${wdm_INCLUDE_DIRS})

# Add main executable
add_executable(main main.cpp)
