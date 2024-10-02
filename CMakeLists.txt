cmake_minimum_required(VERSION 3.12)

project(SphereMeshes)

add_executable(SphereMeshes sphere-meshes.cpp
SQEM.h
)

target_include_directories(SphereMeshes PUBLIC extern/igl
extern/eigen3
)