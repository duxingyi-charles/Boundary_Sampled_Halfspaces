cmake_minimum_required(VERSION 3.15)
project(PSI VERSION 0.1.0 LANGUAGES CXX)

SET(CMAKE_CXX_STANDARD 14)

include(FetchContent)
include(cmake/Eigen3.cmake)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST} main.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_compile_options(${PROJECT_NAME} PRIVATE "-Wall")

