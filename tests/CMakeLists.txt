find_package(Boost COMPONENTS unit_test_framework)

if(${Boost_unit_test_framework_FOUND})

    message("C++ code will tested")

    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(TEST_LIBRARIES ${Boost_LIBRARIES} libRACES)

    set(RACES_TESTS archive fragment fasta_reader context_positions)

    set(ASSET_DIR "${CMAKE_CURRENT_SOURCE_DIR}/assets")
    set(SBS_EXAMPLE "${ASSET_DIR}/v3.3_SBS3_PROFILE.txt")
    if(NOT EXISTS ${SBS_EXAMPLE})
        message("Downloading SBS profile example from sanger...")
        file(DOWNLOAD
            https://cancer.sanger.ac.uk/signatures/documents/1944/v3.3_SBS3_PROFILE.txt
            ${SBS_EXAMPLE}
            STATUS DOWNLOAD_RESULT)

        list(GET DOWNLOAD_RESULT 1 DOWNLOAD_RESULT_MESSAGE)
        if(${DOWNLOAD_RESULT_MESSAGE} STREQUAL "\"No error\"")
            add_executable( test_snv_signature snv_signature.cpp)
            target_link_libraries(test_snv_signature ${TEST_LIBRARIES})
            target_compile_definitions(test_snv_signature PUBLIC SBS_EXAMPLE="${SBS_EXAMPLE}")
        endif()
    else()
        add_executable( test_snv_signature snv_signature.cpp)
        target_link_libraries(test_snv_signature ${TEST_LIBRARIES})
        target_compile_definitions(test_snv_signature PUBLIC SBS_EXAMPLE="${SBS_EXAMPLE}")
    endif()

    foreach(TEST ${RACES_TESTS})
        ADD_EXECUTABLE( test_${TEST} ${TEST}.cpp)
        target_link_libraries(test_${TEST} ${TEST_LIBRARIES})
        ADD_TEST(${TEST} test_${TEST})
    endforeach()

    target_compile_definitions(test_fasta_reader PUBLIC FASTA_FILE="${ASSET_DIR}/fasta_test.fasta")
    target_compile_definitions(test_context_positions PUBLIC FASTA_FILE="${ASSET_DIR}/context_positions_test.fasta")
endif()