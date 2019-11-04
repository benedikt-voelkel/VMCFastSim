# @author Benedikt Volkel
# @brief  cmake setup for FastSim using VMC

set(MODULE_NAME "VMCFastSim")

# Minimum version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.0 FATAL_ERROR)

project(${MODULE_NAME})

# Include for building tests, this generates option BUILD_TESTING which is ON
# by default
include(CTest)

################################################################################
# Extend module file path
################################################################################
set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH}
)

################################################################################
# Basic configurations
################################################################################
# Sources for built are here.
set(CXX_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CXX_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(MACRO_SOURCE_DIR ${CMAKE_SOURCE_DIR}/macro)

# Install paths
set(INSTALL_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(INSTALL_BINARY_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/${MODULE_NAME})
#set(INSTALL_MACRO_DIR ${CMAKE_INSTALL_PREFIX}/macro)
# Collect for configure script
#set(_configure_shared_library_paths)
#set(_configure_static_library_paths)

################################################################################
# Some options to decide about backend installation
################################################################################
#option(WITH_LWTNN "Build with LWTNN backend" ON)
#option(WITH_XGBoost "Build with XGBoost backend" ON)

################################################################################
# Find required packages
################################################################################

########
# ROOT #
########

# NOTE Rely on $ENV{ROOTSYS} at the moment
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core RIO MathCore EG)
include(${ROOT_USE_FILE})
# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

#######
# VMC #
#######

# NOTE Rely on $ENV{ROOTSYS} at the moment
find_package(VMC REQUIRED)
# Include ROOT headers
include_directories(${VMC_INCLUDE_DIRS})


#########
# Other #
#########

# Nothing yet

################################################################################
# Set C++ standard
################################################################################

# NOTE So far only based on ROOT
#message(STATUS "Derive CXX standard from ROOT installation")
SET(CMAKE_CXX_STANDARD 11)
#execute_process(
#    COMMAND $ENV{ROOTSYS}/bin/root-config --cflags
#    OUTPUT_VARIABLE ROOT_C_FLAGS
#    OUTPUT_STRIP_TRAILING_WHITESPACE)
# Find the cxx standard
#string(FIND ${ROOT_C_FLAGS} "-std=" POSITION)
#if (${POSITION} GREATER -1)
#  string(SUBSTRING ${ROOT_C_FLAGS} ${POSITION} 11 ROOT_CXX_STD)
#  if(${ROOT_CXX_STD} STREQUAL "-std=c++1z " OR ${ROOT_CXX_STD} STREQUAL "-std=c++17 ")
#    SET(CMAKE_CXX_STANDARD 17)
#  elseif(${ROOT_CXX_STD} STREQUAL "-std=c++1y " OR ${ROOT_CXX_STD} STREQUAL "-std=c++14 ")
#    SET(CMAKE_CXX_STANDARD 14)
#  endif()
#endif()
#message(STATUS "Compile with ${CMAKE_CXX_STANDARD}")


################################################################################
# Compiler and linker flags
################################################################################
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2 -g")
set(CMAKE_CXX_FLAGS  "-O2")


################################################################################
# Base sources and headers
################################################################################
# HEADERS and LIBS sum everything required for building.
# NOTE So far everything is compiled into one lib.
set(SRCS)
set(HEADERS
   ${CXX_INCLUDE_DIR}/${MODULE_NAME}/FastSim.h
)

################################################################################
# ROOT LinkDef
################################################################################
SET(ROOT_DICT_LINKDEF_FILE "${CXX_SOURCE_DIR}/${MODULE_NAME}LinkDef.h")
message(STATUS "LinkDef at ${ROOT_DICT_LINKDEF_FILE}")
SET(ROOT_DICT_NAME "G__${MODULE_NAME}")

ROOT_GENERATE_DICTIONARY(${ROOT_DICT_NAME} ${HEADERS} LINKDEF ${ROOT_DICT_LINKDEF_FILE})
# Files produced by the dictionary generation
SET(ROOT_DICT_LIB_FILES
    "${PROJECT_BINARY_DIR}/lib${MODULE_NAME}_rdict.pcm"
    "${PROJECT_BINARY_DIR}/lib${MODULE_NAME}.rootmap"
)

################################################################################
# Used for any additional directories and libraries to be included and used
# for built
################################################################################
#set(INCLUDES)
#set(LIBS)
#set(SCRIPTS)


################################################################################
# Check for backends
################################################################################
#include(DeployBackends)


################################################################################
# Setup for library built and installation
################################################################################

include_directories(${CXX_INCLUDE_DIR})
set(LIBRARY_NAME ${MODULE_NAME})
add_library(${LIBRARY_NAME} SHARED ${SRCS} "${ROOT_DICT_NAME}.cxx")
target_link_libraries(${LIBRARY_NAME} ${ROOT_LIBRARIES} ${VMC_LIBRARIES})
list(APPEND _configure_shared_library_paths ${INSTALL_LIBRARY_DIR})

#add_executable(run ${CXX_SOURCE_DIR}/run.cxx)
#target_link_libraries(run ${LIBRARY_NAME})

################################################################################
# Configure the config and setup script
################################################################################
# This project
#set(SETUP_SCRIPT_IN ${CONFIG_SOURCE_DIR}/setup.in)
#set(SETUP_SCRIPT_OUT ${CMAKE_BINARY_DIR}/setup)
#include(BuildConfig)

################################################################################
# Build tests
################################################################################
#if(NOT BUILD_TESTING)
#  message(WARNING "No tests are built.")
#else()
#  add_subdirectory(test)
#endif(NOT BUILD_TESTS)

################################################################################
# Install the project
################################################################################
# Install headers
install(FILES ${HEADERS} DESTINATION ${INSTALL_INCLUDE_DIR})
# Install libraries
install(TARGETS ${LIBRARY_NAME} DESTINATION ${INSTALL_LIBRARY_DIR})
# Install binaries
#install(TARGETS run DESTINATION ${INSTALL_BINARY_DIR})
# Install the ROOT dictionary files
install(FILES ${ROOT_DICT_LIB_FILES} DESTINATION ${INSTALL_LIBRARY_DIR})
