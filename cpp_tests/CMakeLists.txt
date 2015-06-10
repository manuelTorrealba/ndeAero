
include_directories(${NDEAERO_SOURCE_DIR}/cpp_library/include)
link_directories (${NDEAERO_SOURCE_DIR}/cpp_library)

set(SOURCE_FILES main.cpp
					  TestNACA4Digits.cpp
					  TestODESolver.cpp
					  TestThinAirfoil.cpp
					  TestThinAirfoilTheory.cpp)

set(HEADER_FILES Test.hpp
					  TestODESolver.hpp
                 TestThinAirfoil.hpp)

add_executable(tests ${SOURCE_FILES})

target_link_libraries(tests ndeAero)
