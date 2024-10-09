cmake_minimum_required(VERSION 3.12)

project(SphereMeshes)
find_package(ZLIB REQUIRED)
set(ENABLE_STATIC ON)
add_subdirectory("extern/cnpy")
add_executable(SphereMeshes 
    main.cpp
    SphereMeshes.cpp
    DirectionalWidth.cpp
    SQEM.h
    SphereMeshes.h
    DirectionalWidth.h
)

add_executable(test_normal
    tests/test_normal.cpp
    SphereMeshes.cpp
    DirectionalWidth.cpp
    SQEM.h
    SphereMeshes.h
    DirectionalWidth.h
)
target_include_directories(SphereMeshes PUBLIC extern/libigl
    extern/eigen3
    extern/cnpy
    ${ZLIB_INCLUDE_DIRS}
)
target_include_directories(test_normal PUBLIC extern/libigl
    extern/eigen3
)
target_link_libraries(SphereMeshes cnpy)