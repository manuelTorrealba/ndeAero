
include_directories(${NDEAERO_SOURCE_DIR}/cpp_library/include)
link_directories (${NDEAERO_SOURCE_DIR}/cpp_library)

set(SOURCE_FILES main.cpp
					  TestThinAirfoil.cpp)

add_executable(tests ${SOURCE_FILES})

target_link_libraries(tests ndeAero)
