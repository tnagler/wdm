function(add_wdm_test name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} wdm)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_wdm_test(test_cxi test_cxi.cpp)
add_wdm_test(test_estimators test_estimators.cpp)
add_wdm_test(test_input test_input.cpp)
add_wdm_test(test_inference test_inference.cpp)
add_wdm_test(test_ranks test_ranks.cpp)

set(wdm_headers
    wdm.hpp
    wdm/bbeta.hpp
    wdm/cxi.hpp
    wdm/hoeffd.hpp
    wdm/ktau.hpp
    wdm/methods.hpp
    wdm/nan_handling.hpp
    wdm/prho.hpp
    wdm/random.hpp
    wdm/ranks.hpp
    wdm/srho.hpp
    wdm/utils.hpp
)

foreach(header IN LISTS wdm_headers)
    string(REPLACE "/" "_" header_target ${header})
    string(REPLACE "." "_" header_target ${header_target})
    add_executable(test_header_${header_target} test_headers.cpp)
    target_link_libraries(test_header_${header_target} wdm)
    target_compile_definitions(test_header_${header_target}
        PRIVATE WDM_HEADER="${header}")
    add_test(NAME test_header_${header_target}
        COMMAND test_header_${header_target})
endforeach()

find_package(Eigen3 QUIET NO_MODULE)
if(Eigen3_FOUND)
    add_executable(test_header_wdm_eigen_hpp test_headers.cpp)
    target_link_libraries(test_header_wdm_eigen_hpp wdm Eigen3::Eigen)
    target_compile_definitions(test_header_wdm_eigen_hpp
        PRIVATE WDM_HEADER="wdm/eigen.hpp")
    add_test(NAME test_header_wdm_eigen_hpp COMMAND test_header_wdm_eigen_hpp)

    add_executable(test_wdm_eigen test_eigen.cpp)
    target_link_libraries(test_wdm_eigen wdm Eigen3::Eigen)
    add_test(NAME test_wdm_eigen COMMAND test_wdm_eigen)
endif()

set(install_prefix ${PROJECT_BINARY_DIR}/test-install)
set(consumer_build ${PROJECT_BINARY_DIR}/test-consumer)
add_test(NAME test_install
    COMMAND ${CMAKE_COMMAND} --install ${PROJECT_BINARY_DIR}
            --prefix ${install_prefix} --config $<CONFIG>)
set_tests_properties(test_install PROPERTIES FIXTURES_SETUP wdm_installed)

add_test(NAME test_consumer_configure
    COMMAND ${CMAKE_COMMAND}
            -S ${CMAKE_CURRENT_SOURCE_DIR}/install
            -B ${consumer_build}
            -DCMAKE_PREFIX_PATH=${install_prefix})
set_tests_properties(test_consumer_configure PROPERTIES
    FIXTURES_REQUIRED wdm_installed
    FIXTURES_SETUP wdm_consumer_configured)

add_test(NAME test_consumer_build
    COMMAND ${CMAKE_COMMAND} --build ${consumer_build} --config $<CONFIG>)
set_tests_properties(test_consumer_build PROPERTIES
    FIXTURES_REQUIRED wdm_consumer_configured)
