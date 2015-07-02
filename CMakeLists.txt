

cmake_minimum_required (VERSION 2.6)
project(FindTranslocations)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_definitions( -Wno-deprecated )

find_package(Boost COMPONENTS  program_options system filesystem REQUIRED)

# set our library and executable destination dirs
set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )

# define compiler flags for all code
set( CMAKE_BUILD_TYPE Release )

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/lib/")
include_directories("${PROJECT_SOURCE_DIR}/lib/bamtools/src")
include_directories( ${Boost_INCLUDE_DIRS} )



# sorgenti da compilare
file(GLOB FIND_TRANS_FILES
    ${PROJECT_SOURCE_DIR}/src/FindTranslocations.cpp
    ${PROJECT_SOURCE_DIR}/src/data_structures/Translocation.cpp
	${PROJECT_SOURCE_DIR}/src/data_structures/ExtractModule.cpp
	${PROJECT_SOURCE_DIR}/src/data_structures/findTranslocationsOnTheFly.cpp
    ${PROJECT_SOURCE_DIR}/src/common.h
)


add_subdirectory(lib)


# Find_translocations executable
add_executable(FindTranslocations ${FIND_TRANS_FILES})
target_link_libraries(FindTranslocations ${ZLIB_LIBRARIES})
target_link_libraries(FindTranslocations  BamTools)
target_link_libraries(FindTranslocations ${Boost_LIBRARIES})


