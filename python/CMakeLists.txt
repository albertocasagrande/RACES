find_package(Python COMPONENTS Interpreter Development)
find_package(Boost COMPONENTS system python)

set(COMPILE_PYTHON_BINDINGS TRUE)

if(NOT ${Boost_python_FOUND})
message("Boost.Python is missing: Python bindings will NOT be compiled")
set(COMPILE_PYTHON_BINDINGS FALSE)
endif()

if(NOT ${Boost_system_FOUND})
message("Boost.System is missing: Python bindings will NOT be compiled")
set(COMPILE_PYTHON_BINDINGS FALSE)
endif()

if(NOT ${PYTHONLIBS_FOUND})
message("Python libraries are missing: Python bindings will NOT be compiled")
set(COMPILE_PYTHON_BINDINGS FALSE)
endif()

if(${COMPILE_PYTHON_BINDINGS})
    message("Python bindings will be compiled")

    # include directories
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
    include_directories(Python::Module)
    include_directories(${Boost_INCLUDE_DIRS})

    # create the lib
    add_library(pyRACES MODULE
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/bindings.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/genotype.cpp  
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/epigenetic_rates.cpp  
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation_wrapper.cpp)

    set_target_properties(pyRACES PROPERTIES OUTPUT_NAME RACES PREFIX "")

    # link
    target_link_libraries(pyRACES ${Boost_LIBRARIES} Python::Module libRACES)

    find_program(PYTHON_FULL_PATH Python::Interpreter)

    if(EXISTS ${PYTHON_FULL_PATH})
        message("Python bindings will tested")

        enable_testing()
        add_test (NAME python-tests
            COMMAND ${PYTHON_FULL_PATH} -m pytest # Or just COMMAND pytest
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/
            )
        set_tests_properties(python-tests
            PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}:$ENV{PYTHONPATH}")
    endif()
endif()