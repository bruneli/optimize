cmake_minimum_required (VERSION 2.6)
project (optimize)

if(WIN32)
    # Define a custom Boost folder via commands
    set(BOOST_INCLUDEDIR "C:\\boost\\boost_1_54_0\\include\\boost-1_54")
    set(BOOST_LIBRARYDIR "C:\\boost\\boost_1_54_0\\lib")
    # Define a custom Python folder via commands
    set(PYTHON_LIBRARY "C:\\Anaconda\\python27.dll")
    set(PYTHON_INCLUDE_DIR "C:\\Anaconda\\include")
endif()

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
if(MSVC)
    if(NOT DEFINED BOOST_ALL_DYN_LINK)
        set(BOOST_ALL_DYN_LINK TRUE)
    endif()
    if(BOOST_ALL_DYN_LINK)
        message("Set special flags for MSCV")
        add_definitions( -DBOOST_ALL_NO_LIB )
        add_definitions( -DBOOST_ALL_DYN_LINK )
        add_definitions( -DBOOST_PYTHON_SOURCE )
        #add_definitions( -DBOOST_LIB_DIAGNOSTIC )
        add_definitions( -D_SCL_SECURE_NO_WARNINGS )
    endif()
endif()

set(INCLUDE_INSTALL_DIR "${CMAKE_SOURCE_DIR}/optimize")
set(LIBRARY_INSTALL_DIR "${CMAKE_SOURCE_DIR}/lib")
#LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)

enable_testing()
add_test(TestFuncVar ${CMAKE_SOURCE_DIR}/bin/test_func_var)
add_test(TestLineSearch ${CMAKE_SOURCE_DIR}/bin/test_line_search)
add_test(TestMinimizer ${CMAKE_SOURCE_DIR}/bin/test_minimizer)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
