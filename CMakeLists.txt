project(DAS CXX)

cmake_minimum_required(VERSION 2.6)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-O3 --std=c++11" )
set( DAS_EXTRA_INCLUDES )
set( DAS_EXTRA_LIBRARIES )

set( DAS_SRCS 
     src/main.cpp
     src/galaxy.cpp
     src/tridiag.cpp
)

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
  set(HAVE_LIBGSL 1)
  list(APPEND DAS_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
  list(APPEND DAS_EXTRA_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
  message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif(GSL_FOUND)

include_directories (include ${DAS_EXTRA_INCLUDES})
add_executable(DAS ${DAS_SRCS} )
target_link_libraries (DAS ${DAS_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
