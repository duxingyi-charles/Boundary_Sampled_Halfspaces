cmake_minimum_required(VERSION 3.15)
project(PSI VERSION 0.1.0 LANGUAGES CXX)

SET(CMAKE_CXX_STANDARD 14)

include(FetchContent)
include(cmake/Eigen3.cmake)
include(cmake/cli11.cmake)
include(cmake/nlohmann-json.cmake)

option(LIBIGL_WITH_CGAL "Build CGAL support in libigl" On)
include(cmake/libigl.cmake)


aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
add_library(psi STATIC ${SRC_LIST})
target_include_directories(psi PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(psi PUBLIC Eigen3::Eigen)
target_compile_options(psi PUBLIC "-Wall")

add_executable(psi_cli apps/psi_cli.cpp)
target_link_libraries(psi_cli PRIVATE psi CLI11::CLI11 nlohmann::json)

add_library(cell_partition STATIC
    ${PROJECT_SOURCE_DIR}/apps/PyMesh/CellPartition.cpp)
target_link_libraries(cell_partition PUBLIC igl::core igl::cgal)

add_executable(mesh_arrangement apps/mesh_arrangement.cpp)
target_link_libraries(mesh_arrangement psi CLI11::CLI11 nlohmann::json igl::core
    cell_partition)
