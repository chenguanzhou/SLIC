cmake_minimum_required (VERSION 2.8)
project (SLIC)

# The version number.
set (ChangeDetection_VERSION_MAJOR 0)
set (ChangeDetection_VERSION_MINOR 0)

# Color
set(CMAKE_COLOR_MAKEFILE ON)

# Source
set(SRC_LISTS slic.cpp)

# Add Eecutable
add_executable(slic ${SRC_LISTS})

# GDAL Support
FIND_PACKAGE(GDAL)
IF(GDAL_FOUND)
        INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
        TARGET_LINK_LIBRARIES(slic ${GDAL_LIBRARY})
ELSE(GDAL_FOUND)
        MESSAGE(FATEL_ERROR "GDAL library not found!")
ENDIF(GDAL_FOUND)