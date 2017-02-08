cmake_minimum_required (VERSION 2.6)
project (WIsH)

set (WIsH_VERSION_MAJOR 1)
set (WIsH_VERSION_MINOR 0)

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
configure_file ("README.md.in" "README.md")


add_executable(WIsH ${SOURCES})


