cmake_minimum_required(VERSION 3.11)

project (races CXX)

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

include_directories(${CMAKE_SOURCE_DIR}/include)

set(SIM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/driver_genotype.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutation_graphs.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/species.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/position.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/tissue.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/plot_2D.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/driver_palette.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/simulator_main.cpp)

if(${COMPILE_SDL2_CODE})
    add_compile_definitions(WITH_SDL2)

    include_directories(${SDL2_INCLUDE_DIRS} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)

    set(SIM_SOURCES ${SIM_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/SDL_plot.cpp)
endif()

add_executable(races_sim ${SIM_SOURCES})

target_include_directories( races_sim PRIVATE ${Boost_INCLUDE_DIR})

if(${COMPILE_SDL2_CODE})
    target_link_libraries(races_sim ${SDL2_LIBRARIES} ${Boost_LIBRARIES} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)
else()
    target_link_libraries(races_sim ${SDL2_LIBRARIES} ${Boost_LIBRARIES})
endif()
