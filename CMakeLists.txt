cmake_minimum_required (VERSION 2.8)
project (slic)

# The version number.
set (slic_VERSION_MAJOR 0)
set (slic_VERSION_MINOR 0)

# Debug||Release
set(CMAKE_BUILD_TYPE Release)

# Color
set(CMAKE_COLOR_MAKEFILE ON)

# Source
set(SRC_LISTS slic.cpp main.cpp)

# Add Eecutable
add_executable(slic ${SRC_LISTS})

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# GDAL Support
FIND_PACKAGE(GDAL)
IF(GDAL_FOUND)
        INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
        TARGET_LINK_LIBRARIES(slic ${GDAL_LIBRARY})
ELSE(GDAL_FOUND)
        MESSAGE(FATEL_ERROR "GDAL library not found!")
ENDIF(GDAL_FOUND)
