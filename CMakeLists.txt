cmake_minimum_required(VERSION 3.15)
project(PSI VERSION 0.1.0 LANGUAGES CXX C)

SET(CMAKE_CXX_STANDARD 14)

include(FetchContent)
include(cmake/Eigen3.cmake)
include(cmake/cli11.cmake)
include(cmake/nlohmann-json.cmake)
include(cmake/fast_arrangement.cmake)
include(cmake/nlopt.cmake)
include(cmake/armadillo.cmake)

option(LIBIGL_WITH_CGAL "Build CGAL support in libigl" On)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Build with UI" On)
option(LIBIGL_WITH_OPENGL "" On)
option(LIBIGL_WITH_OPENGL_GLFW "" On)
include(cmake/libigl.cmake)

# boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS graph)
message(STATUS "Boost version: ${Boost_VERSION}")

# This is needed if your Boost version is newer than your CMake version
# or if you have an old version of CMake (<3.5)
if(NOT TARGET Boost::graph)
    add_library(Boost::graph IMPORTED INTERFACE)
    set_property(TARGET Boost::graph PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::graph PROPERTY
            INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
#add_library(psi STATIC ${SRC_LIST})
#target_include_directories(psi PUBLIC ${PROJECT_SOURCE_DIR}/src)
#
#target_include_directories(psi PUBLIC ${PROJECT_SOURCE_DIR}/apps)

#target_link_libraries(psi PUBLIC Eigen3::Eigen Boost::graph)
#target_compile_options(psi PUBLIC "-Wall")

# VIPSS

include_directories(${PROJECT_SOURCE_DIR}/src/vipss/surfacer)
#aux_source_directory(. MAIN)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/vipss VIPSS_SRC_LIST)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/vipss/surfacer VIPSS_SURFACER_LIST)

#add_executable(${PROJECT_NAME} ${SRC_LIST} ${MAIN} ${SURFACER_LIST})

add_library(vipss STATIC
        ${VIPSS_SRC_LIST}
        ${VIPSS_SURFACER_LIST})
target_link_libraries(vipss PUBLIC armadillo::armadillo nlopt::nlopt)
set_target_properties(vipss PROPERTIES LINKER_LANGUAGE CXX)


add_library(arrangement STATIC
    ${PROJECT_SOURCE_DIR}/src/PyMesh/Arrangement.cpp
    ${PROJECT_SOURCE_DIR}/src/PyMesh/FastArrangement.cpp
    ${PROJECT_SOURCE_DIR}/src/PyMesh/MeshArrangement.cpp)
target_link_libraries(arrangement PUBLIC igl::core igl::cgal
    fast_arrangement::fast_arrangement)
#target_compile_options(arrangement PUBLIC -O0)

add_library(psi STATIC ${SRC_LIST})
target_include_directories(psi PUBLIC ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/apps)
target_link_libraries(psi PUBLIC Eigen3::Eigen Boost::graph nlohmann::json igl::core arrangement vipss)
target_compile_options(psi PUBLIC "-Wall")
target_compile_definitions(psi PUBLIC "-DLIBIGL_NUM_THREADS=2")

add_executable(mesh_arrangement apps/mesh_arrangement.cpp)
target_link_libraries(mesh_arrangement psi CLI11::CLI11 nlohmann::json igl::core
    arrangement)

add_executable(psi_cli apps/psi_cli.cpp)
target_link_libraries(psi_cli PRIVATE psi CLI11::CLI11)

add_executable(psi_cli_add_new_implicit apps/psi_cli_add_new_implicit.cpp)
target_link_libraries(psi_cli_add_new_implicit PRIVATE psi CLI11::CLI11)

add_executable(psi_cli_sample_reduction apps/psi_cli_sample_reduction.cpp)
target_link_libraries(psi_cli_sample_reduction PRIVATE psi CLI11::CLI11)

#add_executable(psi_cli_connect apps/psi_cli_search_for_connected.cpp)
#target_link_libraries(psi_cli_connect PRIVATE psi CLI11::CLI11)

add_executable(psi_ui apps/psi_ui.cpp apps/UI.h)
target_link_libraries(psi_ui psi CLI11::CLI11 nlohmann::json igl::core igl_opengl_glfw_imgui
    arrangement)

add_executable(vipss_fitting apps/vipss_fitting.cpp)
target_link_libraries(vipss_fitting PRIVATE psi CLI11::CLI11)

add_executable(efficient_ransac_test apps/efficient_ransac_test.cpp)
target_link_libraries(efficient_ransac_test PRIVATE CLI11::CLI11 igl::cgal)

add_executable(efficient_ransac_fitting apps/efficient_ransac_fitting.cpp)
target_link_libraries(efficient_ransac_fitting PRIVATE CLI11::CLI11 igl::cgal)

add_executable(mixed_fitting apps/mixed_fitting.cpp)
target_link_libraries(mixed_fitting PRIVATE psi CLI11::CLI11 igl::cgal)

add_executable(export_sampled_implicits apps/export_sampled_implicits.cpp)
target_link_libraries(export_sampled_implicits PRIVATE psi CLI11::CLI11)
