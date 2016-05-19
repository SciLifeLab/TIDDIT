

cmake_minimum_required (VERSION 2.6)
project(TIDDIT)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_definitions( -Wno-deprecated )

# set our library and executable destination dirs
set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )

# define compiler flags for all code
set( CMAKE_BUILD_TYPE Release )

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/lib/")
include_directories("${PROJECT_SOURCE_DIR}/lib/bamtools/src")


# source code
file(GLOB TIDDIT_FILES
    ${PROJECT_SOURCE_DIR}/src/TIDDIT.cpp
    ${PROJECT_SOURCE_DIR}/src/data_structures/Translocation.cpp
	${PROJECT_SOURCE_DIR}/src/data_structures/findTranslocationsOnTheFly.cpp
    ${PROJECT_SOURCE_DIR}/src/common.h
    ${PROJECT_SOURCE_DIR}/src/data_structures/CoverageModule.cpp
)


add_subdirectory(lib)


# TIDDIT executable
add_executable(TIDDIT ${TIDDIT_FILES})
target_link_libraries(TIDDIT ${ZLIB_LIBRARIES})
target_link_libraries(TIDDIT  BamTools)
