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

message("Downloading \"p-ranav/indicators\" library...")
file(DOWNLOAD
    https://raw.githubusercontent.com/p-ranav/indicators/master/single_include/indicators/indicators.hpp
    ${CMAKE_SOURCE_DIR}/include/indicators.hpp
    STATUS DOWNLOAD_RESULT)
list(GET DOWNLOAD_RESULT 1 DOWNLOAD_RESULT_MESSAGE)
if(${DOWNLOAD_RESULT_MESSAGE} STREQUAL "\"No error\"")
    add_compile_definitions(WITH_INDICATORS)
    message("done")
    message("Using p-ranav/indicators")
else()
    message("Download error: ${DOWNLOAD_RESULT_MESSAGE}")
    message("Avoiding p-ranav/indicators")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)

set(RACES_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/archive.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/driver_genotype.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/species.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/position.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/tissue.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell_event.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/somatic_mutations.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/binary_logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/statistics.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/palette.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/plot_2D.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/progress_bar.cpp)

set(RACES_LIBRARIES ${Boost_LIBRARIES})

if(${COMPILE_SDL2_CODE})
    add_compile_definitions(WITH_SDL2)

    include_directories(${SDL2_INCLUDE_DIRS} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)

    set(RACES_SOURCES ${RACES_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/SDL_plot.cpp)

    set(RACES_LIBRARIES ${RACES_LIBRARIES} ${SDL2_LIBRARIES} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)
endif()

add_library(libRACES SHARED ${RACES_SOURCES})
target_link_libraries(libRACES ${RACES_LIBRARIES})
set_target_properties(libRACES PROPERTIES VERSION ${PROJECT_VERSION} OUTPUT_NAME RACES)

set(SIM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/races_sim.cpp)

set(PHYLO_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/races_phylo.cpp
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/phyloXML.cpp
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/sampler.cpp
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/phylogenetic_forest.cpp)


add_executable(races_sim ${SIM_SOURCES})
target_include_directories(races_sim PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(races_sim libRACES)

add_executable(races_phylo ${PHYLO_SOURCES})
target_include_directories(races_phylo PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(races_phylo libRACES)

add_subdirectory(python)

enable_testing()

add_subdirectory(tests)
