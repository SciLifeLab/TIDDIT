project( BamTools )

# BamTools version information
set (BamTools_VERSION_MAJOR 2)
set (BamTools_VERSION_MINOR 0)
set (BamTools_VERSION_BUILD 5)

# add our includes root path
include_directories( bamtools/src )

# list subdirectories to build in
add_subdirectory( bamtools/src/api )

