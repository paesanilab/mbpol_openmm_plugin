#---------------------------------------------------
# OpenMM MBPol Plugin
#
# Creates OpenMM MBPol plugin library, base name=OpenMMMBPol.
#----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
set (CMAKE_CXX_STANDARD 11)

# We need to know where OpenMM is installed so we can access the headers and libraries.
SET(OPENMM_DIR "/usr/local/openmm" CACHE PATH "Where OpenMM is installed")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include" "${OPENMM_DIR}/include/openmm/reference" "${OPENMM_DIR}/include/openmm//lepton/")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# Select where to install
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    IF(WIN32)
        SET(CMAKE_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install the plugin" FORCE)
    ELSE(WIN32)
        SET(CMAKE_INSTALL_PREFIX "/usr/local/openmm" CACHE PATH "Where to install the     plugin" FORCE)
    ENDIF(WIN32)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(MBPOL_PLUGIN_SOURCE_SUBDIRS openmmapi)

# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(MBPOL_LIBRARY_NAME OpenMMMBPol)
SET(MBPOL_MAJOR_VERSION 1)
SET(MBPOL_MINOR_VERSION 1)
SET(MBPOL_BUILD_VERSION 1)

SET(OPENMM_MAJOR_VERSION 6 CACHE INTEGER "OpenMM major version, e.g. 6 in 6.2.0")
SET(OPENMM_MINOR_VERSION 2 CACHE INTEGER "OpenMM minor version, e.g. 2 in 6.2.0")

ADD_DEFINITIONS(-DOPENMM_MAJOR_VERSION=${OPENMM_MAJOR_VERSION}
                -DOPENMM_MINOR_VERSION=${OPENMM_MINOR_VERSION})

SET(SHARED_MBPOL_TARGET ${MBPOL_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include" "openmmapi/include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${MBPOL_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Create the library.

ADD_LIBRARY(${SHARED_MBPOL_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES}       ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_MBPOL_TARGET} PROPERTIES COMPILE_FLAGS "-DMBPOL_BUILDING_SHARED_LIBRARY")
TARGET_LINK_LIBRARIES(${SHARED_MBPOL_TARGET} optimized OpenMM)
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_MBPOL_TARGET})

# install headers excluding headers in internal
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_OPENMM "openmmapi/include/openmm/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_OPENMM} DESTINATION include/openmm)
FILE(GLOB API_ONLY_INCLUDE_FILES_OPENMM_INTERNAL "openmmapi/include/openmm/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_OPENMM_INTERNAL} DESTINATION include/openmm/internal)

# Enable testing

ENABLE_TESTING()

# Build the implementations for different platforms

ADD_SUBDIRECTORY(platforms/reference)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# Optimized CPU platform

SET(MBPOL_BUILD_CPU_LIB ON CACHE BOOL "Build optimized CPU platform")
IF(MBPOL_BUILD_CPU_LIB)
   ADD_SUBDIRECTORY(platforms/cpu)
ENDIF(MBPOL_BUILD_CPU_LIB)

FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    SET(MBPOL_BUILD_CUDA_LIB ON CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND)
    SET(MBPOL_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND)
IF(MBPOL_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(MBPOL_BUILD_CUDA_LIB)

# Build the Python API


FIND_PROGRAM(PYTHON_EXECUTABLE python)
FIND_PROGRAM(SWIG_EXECUTABLE swig)
IF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(MBPOL_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
ELSE(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(MBPOL_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Build wrappers for Python")
ENDIF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
IF(MBPOL_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(MBPOL_BUILD_PYTHON_WRAPPERS)
