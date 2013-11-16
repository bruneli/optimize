set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS thread python REQUIRED)
message("Include dirs of boost: " ${Boost_INCLUDE_DIRS} )
message("Libs of boost: " ${Boost_LIBRARIES} )

find_package(PythonLibs REQUIRED)
message("Include dirs of Python: " ${PYTHON_INCLUDE_PATH} )
message("Libs of Python: " ${PYTHON_LIBRARIES} )

if(Boost_FOUND AND PYTHONLIBS_FOUND)
    include_directories(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
    add_library(optimizelib SHARED variable.cpp function.cpp minimizer.cpp
                minimizer_nm.cpp minimizer_bfgs.cpp
                line_search.cpp scalar_search.cpp)
    target_link_libraries(optimizelib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
    set_target_properties(optimizelib 
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
    )
    set_target_properties(optimizelib
        PROPERTIES OUTPUT_NAME optimize${BUILD_POSTFIX}
    )
else()
    message(FATAL_ERROR "Unable to find Boost.")
endif()
