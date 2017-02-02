cmake_minimum_required (VERSION 2.6)
project (WIsH)


find_package(OpenMP)
if (OPENMP_FOUND)
	option(OPENMP "OpenMP support" ON)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	option(OPENMP "OpenMP support" OFF)
endif()



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
file(GLOB SOURCES "*.cpp")
configure_file ("main.h.in" "main.h")
add_executable(WIsH ${SOURCES})


