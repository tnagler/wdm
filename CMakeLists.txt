cmake_minimum_required(VERSION 3.1)
project(wdm)

set(VERSION "0.1.0")

# request c++11 without gnu extension for the whole project and enable more warnings
set(CMAKE_CXX_STANDARD 11)

add_library(wdm INTERFACE)

target_include_directories(wdm INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    )

include(CMakePackageConfigHelpers)
set(WDM_CMAKE_DIR "lib/cmake/wdm" CACHE STRING
  "Installation directory for cmake files, relative to ${CMAKE_INSTALL_PREFIX}.")
set(version_config "${PROJECT_BINARY_DIR}/wdm-config-version.cmake")
set(project_config "${PROJECT_BINARY_DIR}/wdm-config.cmake")
set(targets_export_name wdm-targets)

# Generate the version, config and target files into the build directory.
write_basic_package_version_file(
    ${version_config}
    VERSION ${VERSION}
    COMPATIBILITY AnyNewerVersion)
configure_package_config_file(
        ${PROJECT_SOURCE_DIR}/wdm-config.cmake.in
        ${project_config}
    INSTALL_DESTINATION ${WDM_CMAKE_DIR})
export(TARGETS wdm NAMESPACE wdm::
    FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)


# Install version, config and target files.
install(
    FILES  ${version_config} ${project_config}
DESTINATION ${WDM_CMAKE_DIR})
install(EXPORT ${targets_export_name} DESTINATION ${WDM_CMAKE_DIR})

# Install the header file and export the target
install(TARGETS wdm EXPORT ${targets_export_name} DESTINATION lib)

file(GLOB_RECURSE main_hpp ${PROJECT_SOURCE_DIR}/include/wdm.hpp)
file(GLOB_RECURSE impl_hpp ${PROJECT_SOURCE_DIR}/include/wdm/*.hpp)
install(FILES ${main_hpp} DESTINATION include)
install(FILES ${impl_hpp} DESTINATION include/wdm)
