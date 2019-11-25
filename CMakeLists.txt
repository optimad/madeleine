#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project(madeleine)

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
set(ENABLE_MPI 0 CACHE BOOL "If set, the program is compiled with MPI support")
set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")
set(BUILD_SHARED_LIBS 0 CACHE BOOL "Build Shared Libraries")

#------------------------------------------------------------------------------------#
# Internal variables
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
# Internal dependencies
#------------------------------------------------------------------------------------#
include_directories("src")

#------------------------------------------------------------------------------------#
# MPI
#------------------------------------------------------------------------------------#
if (ENABLE_MPI)
    find_package(MPI)

    include_directories("${MPI_INCLUDE_PATH}")

    if(MPI_C_COMPILE_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    endif()

    if(MPI_CXX_COMPILE_FLAGS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    endif()

    if(MPI_CXX_LINK_FLAGS)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    endif()
endif()

#------------------------------------------------------------------------------------#
# Compiler settings
#------------------------------------------------------------------------------------#
set (MADELEINE_DEFINITIONS "")

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)

set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (ENABLE_MPI)
    list (APPEND MADELEINE_DEFINITIONS "ENABLE_MPI=1")
else ()
    list (APPEND MADELEINE_DEFINITIONS "ENABLE_MPI=0")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if (ENABLE_WARNINGS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

if (NOT ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_LESS "2.8.12"))
    add_compile_options("-std=c++11")
else ()
    add_definitions("-std=c++11")
endif ()

# Define a preprocessor macro to recognize debug builds
IF(CMAKE_BUILD_TYPE_LOWER MATCHES "debug" OR CMAKE_BUILD_TYPE_LOWER MATCHES "debinfo")
    list (APPEND MADELEINE_DEFINITIONS "ENABLE_DEBUG=1")
else ()
    list (APPEND MADELEINE_DEFINITIONS "ENABLE_DEBUG=0")
endif ()

IF(NOT CMAKE_BUILD_TYPE_LOWER MATCHES "debug")
    list (APPEND MADELEINE_DEFINITIONS "NDEBUG")
endif ()

# Define an alias for building with scalasca
if (ENABLE_MPI)
    SET(C_FLAGS_INSTRUMENT   "-instrument mpicxx")
    SET(CXX_FLAGS_INSTRUMENT "-instrument mpic")
else ()
    SET(C_FLAGS_INSTRUMENT   "")
    SET(CXX_FLAGS_INSTRUMENT "")
endif ()

if (CMAKE_BUILD_TYPE_LOWER MATCHES "scalasca")
    file(WRITE scalasca_c_compiler
"#!/bin/bash
scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
"
    )

    file(WRITE scalasca_cxx_compiler
"#!/bin/bash
scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
"
    )

    file(INSTALL scalasca_cxx_compiler DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
    file(INSTALL scalasca_c_compiler   DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
    file(REMOVE  scalasca_cxx_compiler)
    file(REMOVE  scalasca_c_compiler)

    set(CMAKE_CXX_COMPILER "${PROJECT_BINARY_DIR}/scalasca_cxx_compiler")
    set(CMAKE_C_COMPILER   "${PROJECT_BINARY_DIR}/scalasca_c_compiler")
endif ()

# Add the definitions
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${MADELEINE_DEFINITIONS})

#------------------------------------------------------------------------------------#
# Customized build types
#------------------------------------------------------------------------------------#

# GNU Profiling
SET(CMAKE_CXX_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the C++ compiler during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_GNUPROFILING)

SET(CMAKE_C_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the C compiler during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_C_FLAGS_GNUPROFILING)

SET(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used for linking binaries during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING)

SET(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the shared libraries linker during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING)

# Scalasca Profiling
SET(CMAKE_CXX_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
    "Flags used by the C++ compiler during Scalasca profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_SCALASCAPROFILING)

SET(CMAKE_C_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
    "Flags used by the C compiler during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_C_FLAGS_SCALASCAPROFILING)

SET(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
    "Flags used for linking binaries during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING)

SET(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
    "Flags used by the shared libraries linker during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING)

# Set build type variable
if(NOT CMAKE_BUILD_TYPE)
    # Set default build type to Debug
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel GNUProfiling ScalascaProfiling."
        FORCE)

    # Set the possible values of build type for the GUI
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo" "GNUProfiling" "ScalascaProfiling")
endif()

#------------------------------------------------------------------------------------#
# Version
#------------------------------------------------------------------------------------#
set (VERSION_HEADER_FILENAME "${PROJECT_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}_version.hpp")
string(TOUPPER "${CMAKE_PROJECT_NAME}_VERSION" VERSION_DEFINE_NAME)

file(READ "${VERSION_HEADER_FILENAME}" header)

string(REGEX REPLACE ".*#[ \t]*define[ \t]*${VERSION_DEFINE_NAME}[ \t]*\"([^\n]*)\".*" "\\1" match "${header}")
set(MADELEINE_VERSION "${match}")

STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${MADELEINE_VERSION}")
set(MADELEINE_MAJOR_VERSION "${match}")

STRING(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${MADELEINE_VERSION}")
set(MADELEINE_MINOR_VERSION "${match}")

STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)(-[0-9A-Za-z-]+)?" "\\1" match "${MADELEINE_VERSION}")
set(MADELEINE_PATCH_VERSION "${match}")

STRING(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" match "${MADELEINE_VERSION}")
if (NOT match STREQUAL "")
    STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" "\\1" match "${MADELEINE_VERSION}")
    set(MADELEINE_TAG_VERSION "${match}")
else ()
    set(MADELEINE_TAG_VERSION "")
endif ()

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#

set(MADELEINE_EXTERNAL_LIBRARIES "")
set(MADELEINE_EXTERNAL_INCLUDE_DIRECTORIES "")

###     BITPIT      ###############################################
find_package(BITPIT REQUIRED)
include(${BITPIT_USE_FILE})
list (APPEND MADELEINE_EXTERNAL_LIBRARIES "${BITPIT_LIBRARIES}")

###      LAPACKE   ################################################
set(LAPACKE_STATIC 0 CACHE BOOL "Require LAPACKE static libraries")
mark_as_advanced(LAPACKE_STATIC)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/LAPACKE/cmake")
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
set(LAPACKE_PRIVATE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/external/LAPACKE/include")
unset(LAPACKEXT_DIR CACHE)
list (APPEND MADELEINE_EXTERNAL_LIBRARIES "${LAPACKE_LIBRARIES}" "${LAPACK_LIBRARIES}")
list (APPEND MADELEINE_EXTERNAL_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIRS}" "${LAPACK_INCLUDE_DIRS}")
list (APPEND MADELEINE_EXTERNAL_INCLUDE_DIRECTORIES "${LAPACKE_PRIVATE_INCLUDE_DIR}" "${LAPACK_INCLUDE_DIRS}")

###     XML2      #################################################
find_package(LibXml2 REQUIRED)
list (APPEND MADELEINE_EXTERNAL_LIBRARIES "${LIBXML2_LIBRARIES}")
list (APPEND MADELEINE_EXTERNAL_INCLUDE_DIRECTORIES "${LIBXML2_INCLUDE_DIR}")

###     MPI      #################################################
if (ENABLE_MPI)
    list (APPEND MADELEINE_EXTERNAL_LIBRARIES "${MPI_CXX_LIBRARIES}")
endif()

###     METIS/PARMETIS    ########################################
if (ENABLE_MPI)
    # ---- parmetis/metis search ---
    set(PARMETIS_DIR "/usr/parmetis-4.0.3/" CACHE PATH "path to parMETIS installation")
    set(METIS_DIR "/usr/metis-5.1.0" CACHE PATH "path to METIS installation")

    list(APPEND MADELEINE_EXTERNAL_INCLUDE_DIRECTORIES "${METIS_DIR}/include")
    list(APPEND MADELEINE_EXTERNAL_INCLUDE_DIRECTORIES "${PARMETIS_DIR}/include")

    find_library(METISLIB NAMES "libmetis.so" "libmetis.a" HINTS "${METIS_DIR}/lib")
    if(NOT METISLIB)
        unset(METISLIB CACHE)
        message(FATAL_ERROR "Cannot found any libmetis.so/libmetis.a in METIS lib installation path")
    endif()
    list(APPEND MADELEINE_EXTERNAL_LIBRARIES ${METISLIB})
    unset(METISLIB CACHE)

    find_library(PARMETISLIB NAMES "libparmetis.so" "libparmetis.a" HINTS "${PARMETIS_DIR}/lib")
    if(NOT PARMETISLIB)
        unset(PARMETISLIB CACHE)
        message(FATAL_ERROR "Cannot found any libparmetis.so/libparmetis.a in parMETIS lib installation path")
    endif()
    list(APPEND MADELEINE_EXTERNAL_LIBRARIES ${PARMETISLIB})
    unset(PARMETISLIB CACHE)
    
else()
    UNSET(METIS_DIR CACHE)
    UNSET(PARMETIS_DIR CACHE)
endif ()

include_directories(${MADELEINE_EXTERNAL_INCLUDE_DIRECTORIES})


#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Main program
ADD_SUBDIRECTORY(src)

# Docs
add_subdirectory(doc)

# Examples
add_subdirectory(examples)
