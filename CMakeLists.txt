# Open Data Dector project
#
# (c) 2021 CERN for the benefit of the ODD project
#
# Mozilla Public License Version 2.0

# Set up the project.
cmake_minimum_required( VERSION 3.14 )
project( OpenDataDetector VERSION 0.1 LANGUAGES CXX )

# Standard CMake include(s).
include( GNUInstallDirs )

# Find the DD4Hep package with corresponding library
find_package(DD4hep 1.11 REQUIRED CONFIG COMPONENTS DDCore DDDetectors)

# Switch the Acts extension on
option(ODD_ACTS_EXTENSION "Build within an ACTS context" OFF)

# the `<project_name>_VERSION` variables set by `setup(... VERSION ...)` have
# only local scope, i.e. they are not accessible her for dependencies added
# via `add_subdirectory`. this overrides the `project(...)` funcion for
# sub-projects such that the resulting `<project_name>_VERSION` has
# global scope and is accessible within the main project later on.
cmake_policy(SET CMP0048 NEW)

# Add the main factory 
add_subdirectory(factory)
