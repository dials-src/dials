
if (CMAKE_CXX_COMPILER_ID STREQUAL Clang OR CMAKE_CXX_COMPILER_ID STREQUAL AppleClang)
    set(CMAKE_CXX_FLAGS_COVERAGE            "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-instr-generate -fcoverage-mapping")
    set(GCC_DEBUG_FLAGS "-g -Wall")

    set(CMAKE_CXX_FLAGS_COVERAGE "${GCC_DEBUG_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS_COVERAGE "${GCC_DEBUG_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_COVERAGE} ${CMAKE_C_FLAGS_COVERAGE}")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    set(CMAKE_CXX_FLAGS_COVERAGE            "${CMAKE_CXX_FLAGS_DEBUG} --coverage -fprofile-abs-path")
    set(CMAKE_EXE_LINKER_FLAGS_COVERAGE     "${CMAKE_EXE_LINKER_FLAGS_COVERAGE} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE  "${CMAKE_SHARED_LINKER_FLAGS_COVERAGE} --coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_COVERAGE}")
else()
    # It's an error if we requested coverage, otherwise ignore
    string( TOUPPER "${CMAKE_BUILD_TYPE}" _build_type )
    if (_build_type STREQUAL COVERAGE)
        message(FATAL_ERROR "Did not know how to run Coverage build for ${CMAKE_CXX_COMPILER_ID}")
    endif()
    message(VERBOSE "Did not know how to run Coverage build for ${CMAKE_CXX_COMPILER_ID}")
    return()
endif()

mark_as_advanced(
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE
)

# Add a set of strings to the build_type metadata description, if not present
function(_add_build_types_to_cache_strings)
    get_property(_config_list CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS)
    foreach(build_type IN ITEMS "${ARGN}")
        if (NOT "${build_type}" IN_LIST _config_list)
            list(APPEND _config_list "${build_type}")
        endif()
    endforeach()
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${_config_list}")
endfunction()

# Add the "Coverage" build type to GUI build options
_add_build_types_to_cache_strings("Coverage")
