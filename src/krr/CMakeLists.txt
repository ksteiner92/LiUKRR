cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0015 NEW)

# Project name
project (liukrr)

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

# Configure Doxygen
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})
if(BUILD_DOCUMENTATION)
    set(DOT_PATH, ${DOXYGEN_DOT_PATH})
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64 -m64 -fopenmp")

# Add include directories
include_directories(include ../../lib /opt/intel/mkl/include/ ../../ext/matio/src)
link_directories(../../lib)

#add_executable(krr main.cc machine.cc descriptor.cc dataset.cc utiles.cc progargparser.cc)
add_library(liukrr SHARED machine.cc descriptor.cc dataset.cc utiles.cc progargparser.cc)
SET (CMAKE_C_COMPILER "/usr/bin/gcc")
SET (CMAKE_CXX_COMPILER "/usr/bin/g++")

#link libraries
target_link_libraries(liukrr mkl_intel_ilp64 mkl_core mkl_gnu_thread pthread m dl matio bz2)