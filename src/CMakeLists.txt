find_package(Boost COMPONENTS thread date_time system python)
message("Include dirs of boost: " ${Boost_INCLUDE_DIRS} )
message("Libs of boost: " ${Boost_LIBRARIES} )

if(WIN32)
    set(configuration STATIC)
else()
    set(configuration SHARED)
    find_package(PythonLibs REQUIRED)
    message("Include dirs of Python: " ${PYTHON_INCLUDE_PATH} )
    message("Libs of Python: " ${PYTHON_LIBRARIES} )
endif()
message("Library configuration: " ${configuration} )

if(Boost_FOUND)
    include_directories(
        ${Boost_INCLUDE_DIRS} 
        ${PYTHON_INCLUDE_PATH}
        ${CMAKE_CURRENT_BINARY_DIR}
    )
    if (PYTHONLIBS_FOUND AND Boost_PYTHON_FOUND)
        message("Compile liboptimize with python bindings.")
        add_library(optimizelib ${configuration}
                    variable.cpp function.cpp minimizer.cpp
                    minimizer_nm.cpp minimizer_bfgs.cpp
                    line_search.cpp scalar_search.cpp
                    pybindings.cpp)
        target_link_libraries(optimizelib ${Boost_LIBRARIES}
                              ${PYTHON_LIBRARIES})
    else()
        message("Compile liboptimize without python bindings.")
        add_library(optimizelib ${configuration}
                    variable.cpp function.cpp minimizer.cpp
                    minimizer_nm.cpp minimizer_bfgs.cpp
                    line_search.cpp scalar_search.cpp
                    )
        target_link_libraries(optimizelib ${Boost_LIBRARIES})
    endif()
    message("library install dir ${LIBRARY_INSTALL_DIR}")
    #install(TARGETS optimizelib DESTINATION ${LIBRARY_INSTALL_DIR})
    #install(FILES optimize_export.hpp DESTINATION ${INCLUDE_INSTALL_DIR})
    set_target_properties(optimizelib 
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_INSTALL_DIR}
    )
    if(WIN32)
        set_target_properties(optimizelib 
            PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${LIBRARY_INSTALL_DIR}
        )
        set_target_properties(optimizelib 
            PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${LIBRARY_INSTALL_DIR}
        )
    endif()
    set_target_properties(optimizelib
        PROPERTIES OUTPUT_NAME optimize${BUILD_POSTFIX}
    )
else()
    message(FATAL_ERROR "Unable to find Boost.")
endif()
