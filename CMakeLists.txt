cmake_minimum_required(VERSION 3.11)

project(races VERSION 0.1.0
              DESCRIPTION "RACES is an Advanced Cancer Evolution Simulation"
              HOMEPAGE_URL https://github.com/albertocasagrande/RACES
              LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -Wpedantic)

find_package( Boost COMPONENTS program_options REQUIRED )

find_package(SDL2)
find_package(SDL2_image)
find_package(SDL2_ttf)

set(COMPILE_SDL2_CODE TRUE)
if(${SDL2_ttf_FOUND})
    message("SDL2_ttf found")
elseif()
    message("SDL2_ttf NOT found!")
    set(COMPILE_SDL2_CODE FALSE)
endif()

if(${SDL2_image_FOUND})
    message("SDL2_image found")
elseif()
    message("SDL2_image NOT found!")
    set(COMPILE_SDL2_CODE FALSE)
endif()

if(NOT EXISTS ${CMAKE_SOURCE_DIR}/include/indicators.hpp)
    message("Downloading \"p-ranav/indicators\" library...")
    file(DOWNLOAD
        https://raw.githubusercontent.com/p-ranav/indicators/master/single_include/indicators/indicators.hpp
        ${CMAKE_SOURCE_DIR}/include/indicators.hpp
        STATUS DOWNLOAD_RESULT)
    list(GET DOWNLOAD_RESULT 1 DOWNLOAD_RESULT_MESSAGE)
    if(${DOWNLOAD_RESULT_MESSAGE} STREQUAL "\"No error\"")
        add_compile_definitions(WITH_INDICATORS)
        message("done")
        message("Using \"p-ranav/indicators\"")
    else()
        message("Download error: ${DOWNLOAD_RESULT_MESSAGE}")
        message("Avoiding \"p-ranav/indicators\"")
    endif()
else()
    add_compile_definitions(WITH_INDICATORS)
    message("Using \"p-ranav/indicators\"")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)

set(RACES_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/basic_io.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/archive.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/driver_genotype.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/species.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/position.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/tissue.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell_event.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/timed_mutations.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/binary_logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/statistics.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/palette.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/plot_2D.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/progress_bar.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/context.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/snv_signature.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/snv.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/fragment.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/genomic_region.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/genomic_position.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/fasta_reader.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/fasta_utils.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/genome_mutations.cpp)

set(RACES_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/basic_io.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/archive.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/driver_genotype.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/species.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/position.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/simulation.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/tissue.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/cell.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/cell_event.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/timed_mutations.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/logger.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/binary_logger.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/statistics.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/palette.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/plot_2D.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/progress_bar.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/phyloXML.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/phylogenetic_forest.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/palette.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sampler.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/statistics.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/tissue_plotter.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/time.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/simulation.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/context.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/context_index.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/snv_signature.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/snv.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/fragment.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/genomic_region.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/genomic_position.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/fasta_reader.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/fasta_utils.hpp)

set(RACES_LIBRARIES ${Boost_LIBRARIES})

if(${COMPILE_SDL2_CODE})
    add_compile_definitions(WITH_SDL2)

    include_directories(${SDL2_INCLUDE_DIRS} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)

    set(RACES_SOURCES ${RACES_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/SDL_plot.cpp)
    set(RACES_HEADER ${RACES_HEADER} ${CMAKE_CURRENT_SOURCE_DIR}/include/SDL_plot.hpp)

    set(RACES_LIBRARIES ${RACES_LIBRARIES} ${SDL2_LIBRARIES} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)
endif()

add_library(libRACES SHARED ${RACES_SOURCES})
        
target_link_libraries(libRACES ${RACES_LIBRARIES})
set_target_properties(libRACES PROPERTIES VERSION ${PROJECT_VERSION} OUTPUT_NAME RACES)
target_sources(libRACES PRIVATE ${RACES_SOURCES}
                        PUBLIC FILE_SET HEADERS
                        BASE_DIRS ${PROJECT_SOURCE_DIR}
                        FILES ${RACES_HEADER})

install(TARGETS libRACES FILE_SET HEADERS)

set(SIM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/races_sim.cpp)

set(PHYLO_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/races_phylo.cpp
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/phyloXML.cpp
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/sampler.cpp
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/phylogenetic_forest.cpp)


set(BUILD_INDEX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/races_build_index.cpp)

add_executable(races_sim ${SIM_SOURCES})
target_include_directories(races_sim PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(races_sim libRACES)

add_executable(races_phylo ${PHYLO_SOURCES})
target_include_directories(races_phylo PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(races_phylo libRACES)

add_executable(races_build_index ${BUILD_INDEX_SOURCES})
target_include_directories(races_build_index PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(races_build_index libRACES)

add_subdirectory(python)

enable_testing()

add_subdirectory(tests)
