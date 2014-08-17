
set(SOURCE_FILES AerodynamicBody2D.cpp
                 Airfoil.cpp
                 Interpolator.cpp
                 Panel.cpp
                 testVector.cpp)

set(HEADER_FILES AerodynamicBody2D.hpp
                 Airfoil.hpp
                 Interpolator.hpp
		 Matrix.hpp
                 Panel.hpp
		 Vector.hpp
		 SmartPtr.hpp)

add_executable(tests ${SOURCE_FILES} ${HEADER_FILES}) 
