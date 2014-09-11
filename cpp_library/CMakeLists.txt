

set(INCLUDE_DIRECTORY include)
include_directories(${INCLUDE_DIRECTORY})

set(SOURCE_FILES AerodynamicBody2D.cpp
                 Airfoil.cpp
                 Interpolator.cpp
                 MathUtils.cpp
                 Panel.cpp
		 Wing.cpp)

set(HEADER_FILES ${INCLUDE_DIRECTORY}/AerodynamicBody2D.hpp
                 ${INCLUDE_DIRECTORY}/Airfoil.hpp
                 ${INCLUDE_DIRECTORY}/Interpolator.hpp
                 ${INCLUDE_DIRECTORY}/MathUtils.hpp
		 ${INCLUDE_DIRECTORY}/Matrix.hpp
                 ${INCLUDE_DIRECTORY}/Panel.hpp
		 ${INCLUDE_DIRECTORY}/Vector.hpp
		 ${INCLUDE_DIRECTORY}/Wing.hpp
		 ${INCLUDE_DIRECTORY}/SmartPtr.hpp)

add_library(ndeAero STATIC ${SOURCE_FILES} ${HEADER_FILES}) 
