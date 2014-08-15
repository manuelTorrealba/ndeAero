cmake_minimum_required (VERSION 2.6)
project (NDEAERO)

#add_subdirectory(src)
include_directories(${NDEAERO_SOURCE_DIR}src)

set(SOURCE_FILES src/AerodynamicBody2D.cpp
                 src/Airfoil.cpp
                 src/Interpolator.cpp
                 src/Panel.cpp
                 src/testVector.cpp)
 
add_executable (tests ${SOURCE_FILES}) 
