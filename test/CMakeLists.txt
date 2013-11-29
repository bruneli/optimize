find_library(optimizelib NAMES optimize HINTS ${CMAKE_SOURCE_DIR}/lib 
    ${CMAKE_SOURCE_DIR}/lib/Debug ${CMAKE_SOURCE_DIR}/lib/Release})
find_package(Boost COMPONENTS unit_test_framework REQUIRED) 
if(Boost_FOUND AND optimizelib)
    message("optimize library is ${optimizelib}")
    include_directories(${Boost_INCLUDE_DIRS})
    add_definitions (-DBOOST_TEST_DYN_LINK)

    add_executable(test_func_var test_func_var.cpp)
    target_link_libraries(test_func_var ${Boost_LIBRARIES} ${optimizelib})
    set_target_properties(test_func_var
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
    set_target_properties(test_func_var
        PROPERTIES OUTPUT_NAME test_func_var${BUILD_POSTFIX}
    )

    add_executable(test_line_search test_line_search.cpp)
    target_link_libraries(test_line_search ${Boost_LIBRARIES} ${optimizelib})
    set_target_properties(test_line_search
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
    set_target_properties(test_line_search
        PROPERTIES OUTPUT_NAME test_line_search${BUILD_POSTFIX}
    )

    add_executable(test_minimizer test_minimizer.cpp)
    target_link_libraries(test_minimizer ${Boost_LIBRARIES} ${optimizelib})
    set_target_properties(test_minimizer
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
    set_target_properties(test_minimizer
        PROPERTIES OUTPUT_NAME test_minimizer${BUILD_POSTFIX}
    )
endif()
