cmake_minimum_required(VERSION 3.8)

project(CentroidFold CXX C)

add_definitions(-DPACKAGE_VERSION="0.0.17")

find_package(PkgConfig REQUIRED)
pkg_check_modules(VIENNARNA REQUIRED RNAlib2)
add_definitions(-DHAVE_LIBRNA)
add_definitions(-DHAVE_VIENNA18)
add_definitions(-DHAVE_VIENNA20)
add_definitions(-DHAVE_NAVIEW_XY_COORDINATES)
add_definitions(-DHAVE_MKSTEMP)

find_package(Boost 1.34 REQUIRED
              COMPONENTS program_options random system)

add_subdirectory(src)