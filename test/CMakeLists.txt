add_executable(test_wdm test.cpp)
target_link_libraries(test_wdm wdm)

find_package(Eigen3 QUIET NO_MODULE)
if(Eigen3_FOUND)
    add_executable(test_wdm_eigen test_eigen.cpp)
    target_link_libraries(test_wdm_eigen wdm Eigen3::Eigen)
endif()
