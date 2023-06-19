cmake_minimum_required(VERSION 3.11)

project (races CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -Wpedantic)

find_package(SDL2 REQUIRED)
find_package(SDL2_image REQUIRED)
find_package(SDL2_ttf REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include ${SDL2_INCLUDE_DIRS} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)

set(SIM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/driver_genotype.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutation_graphs.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/species.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/position.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/tissue.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/SDL_plot.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/plot_2D.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/driver_palette.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/simulator_main.cpp)

add_executable(races_sim ${SIM_SOURCES})

target_link_libraries(races_sim ${SDL2_LIBRARIES} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)
