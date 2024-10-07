cmake_minimum_required(VERSION 3.12)

project(SphereMeshes)

add_executable(SphereMeshes 
    main.cpp
    SphereMeshes.cpp
    DirectionalWidth.cpp
    SQEM_bundle.cpp
    SQEM.h
    SphereMeshes.h
    DirectionalWidth.h
    SQEM_bundle.h
)

add_executable(test_normal
    tests/test_normal.cpp
    SphereMeshes.cpp
    DirectionalWidth.cpp
    SQEM_bundle.cpp
    SQEM.h
    SphereMeshes.h
    DirectionalWidth.h
    SQEM_bundle.h
)
target_include_directories(SphereMeshes PUBLIC extern/igl
    extern/eigen3
)
target_include_directories(test_normal PUBLIC extern/igl
    extern/eigen3
)