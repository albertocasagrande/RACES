find_package(Boost COMPONENTS unit_test_framework)

if(${Boost_unit_test_framework_FOUND})

    message("C++ code will tested")

    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(TEST_LIBRARIES ${Boost_LIBRARIES} libRACES)

    set(RACES_TESTS archive)
    foreach(TEST ${RACES_TESTS})
        ADD_EXECUTABLE( test_${TEST} ${TEST}.cpp)
        target_link_libraries(test_${TEST} ${TEST_LIBRARIES})
        ADD_TEST(${TEST} test_${TEST})
    endforeach()
endif()