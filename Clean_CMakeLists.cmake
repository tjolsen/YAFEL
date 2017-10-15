cmake_minimum_required(VERSION 3.5)
project(yafel VERSION 0.1 LANGUAGES CXX)

# For installing properly
include(GNUInstallDirs)

# Configuring and including dependencies
include("${CMAKE_SOURCE_DIR}/cmake/yafel_dependencies.cmake")