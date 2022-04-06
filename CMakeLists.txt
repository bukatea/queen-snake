cmake_minimum_required(VERSION 3.12)
project(QueenSnake)

add_definitions("-std=c++11")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set(CMAKE_CXX_FLAGS "-O3 -Wall")

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB QueenSnake_SRC CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(QueenSnake main.cpp ${QueenSnake_SRC})
target_link_libraries(QueenSnake ${CONAN_LIBS})