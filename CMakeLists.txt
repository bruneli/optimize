cmake_minimum_required (VERSION 2.6)
project (optimize)

LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)

enable_testing()
add_test(TestFuncVar ${CMAKE_SOURCE_DIR}/bin/test_func_var)
add_test(TestLineSearch ${CMAKE_SOURCE_DIR}/bin/test_line_search)
add_test(TestMinimizer ${CMAKE_SOURCE_DIR}/bin/test_minimizer)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
