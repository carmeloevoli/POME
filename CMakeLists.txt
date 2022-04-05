cmake_minimum_required(VERSION 3.2)

file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION VERSION)
set(LANGS CXX)
project(POME VERSION ${VERSION} LANGUAGES ${LANGS})

set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")
set(CMAKE_CXX_FLAGS "-Wall -O3 -ftree-vectorize -std=c++14")
set(POME_EXTRA_INCLUDES)
set(POME_EXTRA_LIBRARIES)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL " ")

include_directories("${CMAKE_SOURCE_DIR}")
include_directories(include/ )
include_directories(external/plog/include/)
include_directories(external/json/include/)

set(POME_SRCS 
     src/model/escape.cpp
     src/model/losses.cpp
     src/model/spectrum.cpp
     src/common.cpp
     src/modelState.cpp
     src/timer.cpp
     src/tridiag.cpp
     src/utilities.cpp
)

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
     set(HAVE_LIBGSL 1)
     list(APPEND POME_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
     list(APPEND POME_EXTRA_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
     message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif(GSL_FOUND)

# copy data files in build
file(GLOB POME_DATA "data/*")
file(COPY ${POME_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# make library
include_directories(include ${POME_EXTRA_INCLUDES})
add_library(POME_LIB ${POME_SRCS} ${POME_EXTRA_INCLUDES})

# make executables
#add_executable(pome src/main.cpp)
#target_link_libraries (pome POME_LIB ${POME_EXTRA_LIBRARIES})

add_executable(print_sourcespectrum examples/printSourceSpectrum.cpp)
target_link_libraries (print_sourcespectrum POME_LIB ${POME_EXTRA_LIBRARIES})

add_executable(print_timescales examples/printTimescales.cpp)
target_link_libraries (print_timescales POME_LIB ${POME_EXTRA_LIBRARIES})

add_executable(evolve_electrons examples/evolvingElectrons.cpp)
target_link_libraries (evolve_electrons POME_LIB ${POME_EXTRA_LIBRARIES})

add_executable(evolve_electrons_cn examples/evolvingElectronsCN.cpp)
target_link_libraries (evolve_electrons_cn POME_LIB ${POME_EXTRA_LIBRARIES})

add_executable(evolve_pwn examples/evolvingPwn.cpp)
target_link_libraries (evolve_pwn POME_LIB ${POME_EXTRA_LIBRARIES})



# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "POME_INCLUDES ........ = ${POME_EXTRA_INCLUDES}")
message (STATUS "POME_LIBRARIES ....... = ${POME_EXTRA_LIBRARIES}")
