cmake_minimum_required(VERSION 3.12)

project(SphereMeshes)

add_executable(SphereMeshes 
    SphereMeshes.cpp
    DirectionalWidth.cpp
    SQEM.h
    DirectionalWidth.h
)

target_include_directories(SphereMeshes PUBLIC extern/igl
extern/eigen3
)