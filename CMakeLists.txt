cmake_minimum_required(VERSION 3.15)
project(PSI VERSION 0.1.0 LANGUAGES CXX)

SET(CMAKE_CXX_STANDARD 11)

aux_source_directory(. MAIN)
aux_source_directory(./src SRC_LIST)

add_executable(${PROJECT_NAME} ${SRC_LIST} ${MAIN} src/Hermite_RBF_sImplicit.cpp src/Hermite_RBF_sImplicit.h src/Grid.h)

# Eigen
set(EIGEN_INCLUDE_DIR /usr/local/include CACHE PATH "Eigen3 include directory")
target_include_directories(${PROJECT_NAME} PUBLIC "${EIGEN_INCLUDE_DIR}")

target_include_directories(${PROJECT_NAME} PUBLIC ${SRC_LIST})

