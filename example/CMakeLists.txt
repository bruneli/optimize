find_library(optimizelib NAMES optimize HINTS ${CMAKE_SOURCE_DIR}/lib 
    ${CMAKE_SOURCE_DIR}/lib/Debug ${CMAKE_SOURCE_DIR}/lib/Release})
find_package(Boost COMPONENTS thread date_time system REQUIRED)
if(Boost_FOUND AND optimizelib)
    message("optimize library is ${optimizelib}")
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(ex_find_mle ex_find_mle.cpp)
    target_link_libraries(ex_find_mle ${optimizelib})
    set_target_properties(ex_find_mle
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
    set_target_properties(ex_find_mle
        PROPERTIES OUTPUT_NAME ex_find_mle${BUILD_POSTFIX}
    )

    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(ex_mt ex_mt.cpp)
    target_link_libraries(ex_mt ${Boost_LIBRARIES} ${optimizelib})
    set_target_properties(ex_mt
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
    set_target_properties(ex_mt
        PROPERTIES OUTPUT_NAME ex_mt${BUILD_POSTFIX}
    )
endif()
