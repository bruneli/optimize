set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF)
#set(BOOST_ROOT /usr/local/boost/boost_1_54_0)
#set(BOOST_INCLUDEDIR /usr/local/boost/boost_1_54_0/include)
#set(BOOST_LIBRARYDIR /usr/local/boost/boost_1_54_0/lib)
find_package(Boost COMPONENTS thread REQUIRED) 
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(ex_find_mle ex_find_mle.cpp)
    target_link_libraries(ex_find_mle optimize)
    set_target_properties(ex_find_mle
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
    set_target_properties(ex_find_mle
        PROPERTIES OUTPUT_NAME ex_find_mle${BUILD_POSTFIX}
    )

    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(ex_mt ex_mt.cpp)
    target_link_libraries(ex_mt ${Boost_LIBRARIES} optimize)
    set_target_properties(ex_mt
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
    set_target_properties(ex_mt
        PROPERTIES OUTPUT_NAME ex_mt${BUILD_POSTFIX}
    )
endif()
