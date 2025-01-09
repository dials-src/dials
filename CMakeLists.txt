cmake_minimum_required(VERSION 3.20...3.30 FATAL_ERROR)

project(dials)

# Add the included modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

# General cmake environment configuration
include(SetDefaultBuildRelWithDebInfo) # Default builds to release with debug info
include(CoverageBuildConfiguration) # Custom module to make turning coverage on easy
include(AlwaysColourCompilation) # Always show coloured compiler output
include(CheckCXXSymbolExists)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Generate compile_commands.json
set(CMAKE_CXX_STANDARD 14)

message(STATUS "CMAKE CXX COMPILER ID - ${CMAKE_CXX_COMPILER_ID}")

find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(CCTBX COMPONENTS cctbx scitbx ccp4io annlib REQUIRED)
find_package(DXTBX REQUIRED)
find_package(OpenMP)
# Handle msgpack.. between versions 3 and 6, split into separate packages.
# While we still support both, look for both names
find_package(msgpack-cxx REQUIRED NAMES msgpack-cxx msgpack)
if (TARGET msgpackc AND NOT TARGET msgpack-cxx)
    # If the old version of the library, let us refer to the new name
    add_library(msgpack-cxx ALIAS msgpackc)
endif()

# Find the boost::python library for this version of python
set(Boost_USE_STATIC_LIBS OFF) # This is the default everywhere except Windows
find_package(Boost COMPONENTS thread "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}" REQUIRED)

# Create Boost::python alias so we don't need to carry the python version around
if(NOT TARGET Boost::python)
    add_library(Boost::python INTERFACE IMPORTED)
    set_target_properties(Boost::python PROPERTIES INTERFACE_LINK_LIBRARIES Python::Module)
endif()

# Put the libraries into lib/ so that we can run this in-place in a TBX install
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/lib")

# We depend on including header sources from above this current directory.
# This currently means this _must_ be in a "dials" folder. This is unstable,
# but will be dealt with when we move to src/ layout.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

link_libraries(DXTBX::DXTBX)
# Make boost a bit quieter - all these warnings come from inside boost.python
add_compile_definitions(BOOST_BIND_GLOBAL_PLACEHOLDERS BOOST_ALLOW_DEPRECATED_HEADERS)

add_subdirectory(src/dials)
add_subdirectory(tests)

# Install all of the DIALS extension modules
install(
    TARGETS
        dials_algorithms_background_ext
        dials_algorithms_background_glm_ext
        dials_algorithms_background_gmodel_ext
        dials_algorithms_background_modeller_ext
        dials_algorithms_background_simple_ext
        dials_algorithms_centroid_simple_ext
        dials_algorithms_filter_ext
        dials_algorithms_image_centroid_ext
        dials_algorithms_image_connected_components_ext
        dials_algorithms_image_fill_holes_ext
        dials_algorithms_image_filter_ext
        dials_algorithms_image_threshold_ext
        dials_algorithms_image_distortion_ext
        dials_algorithms_indexing_ext
        dials_algorithms_integration_bayes_ext
        dials_algorithms_integration_ext
        dials_algorithms_integration_fit_ext
        dials_algorithms_integration_integrator_ext
        dials_algorithms_integration_kapton_ext
        dials_algorithms_integration_parallel_integrator_ext
        dials_algorithms_integration_sum_ext
        dials_algorithms_polygon_clip_ext
        dials_algorithms_polygon_ext
        dials_algorithms_polygon_spatial_interpolation_ext
        dials_algorithms_profile_model_ellipsoid_ext
        dials_algorithms_profile_model_gaussian_rs_ext
        dials_algorithms_profile_model_gaussian_rs_transform_ext
        dials_algorithms_profile_model_modeller_ext
        dials_algorithms_shoebox_ext
        dials_algorithms_simulation_ext
        dials_algorithms_spot_finding_ext
        dials_algorithms_spot_prediction_ext
        dials_algorithms_statistics_ext
        dials_array_family_flex_ext
        dials_model_data_ext
        dials_pychef_ext
        dials_refinement_helpers_ext
        dials_scaling_ext
        dials_tof_scaling_ext
        dials_util_ext
        dials_util_streambuf_test_ext
        dials_viewer_ext
        recviewer_ext
    DESTINATION "${Python_SITEARCH}"
)
