cmake_minimum_required(VERSION 3.27)

project(NMIF)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(EIGEN_DIR "/usr/local/include/eigen3")
INCLUDE_DIRECTORIES ( "${EIGEN_DIR}" )

set(CMAKE_CXX_STANDARD 20)

add_executable(NMIF main.cpp
        MeshClass.cpp
        MeshClass.h
        linAlg.cpp
        linAlg.h
)
