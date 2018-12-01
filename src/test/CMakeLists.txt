cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0015 NEW)

# Project name
project (test)

# Check the C++14 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++1y" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64 -m64 -fopenmp")

# Add include directories
include_directories(include ../../lib /opt/intel/mkl/include/ ../../ext/matio/src ../krr)
link_directories(../../lib)
	
add_executable(error error.cc)
add_executable(analyzer analyzer.cc)
add_executable(sine sineplot.cc)
add_executable(unit unit.cc)
add_executable(krr main.cc)
SET (CMAKE_C_COMPILER "/usr/bin/gcc")
SET (CMAKE_CXX_COMPILER "/usr/bin/g++")

#link libraries
target_link_libraries(analyzer mkl_intel_ilp64 mkl_core mkl_gnu_thread pthread m dl matio bz2 liukrr)
target_link_libraries(krr mkl_intel_ilp64 mkl_core mkl_gnu_thread pthread m dl matio bz2 liukrr)
target_link_libraries(error mkl_intel_ilp64 mkl_core mkl_gnu_thread pthread m dl matio bz2 liukrr)
target_link_libraries(sine mkl_intel_ilp64 mkl_core mkl_gnu_thread pthread m dl matio bz2 liukrr)