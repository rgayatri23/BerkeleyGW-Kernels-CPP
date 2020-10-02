cmake_minimum_required (VERSION 3.12)

# don't allow in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(STATUS "Warning! Building from the source directory is not allow")
  message(STATUS "Remove 'CMakeCache.txt' and 'CMakeFiles' and build from a separate directory")
  message(ERROR "In-source build")
endif()

project(
    BerkeleyGW-kernels
    LANGUAGES C CXX)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS ON)

# creates an option and save the variable name
MACRO(GPP_ADD_OPTION _VAR _MSG _DEF)
    LIST(APPEND GPP_OPTIONS ${_VAR})
    OPTION(${_VAR} "${_MSG}" ${_DEF})
ENDMACRO()

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

GPP_ADD_OPTION(GPP_BUILD_FF "Build FF" OFF)
GPP_ADD_OPTION(GPP_BUILD_TESTING "Enable testing" OFF)
GPP_ADD_OPTION(GPP_USE_OPENMP "Enable OpenMP support" ON)
GPP_ADD_OPTION(GPP_USE_OPENMP_TARGET "Enable OpenMP target support" OFF)
GPP_ADD_OPTION(GPP_USE_TIMEMORY "Enable timemory support" OFF)
GPP_ADD_OPTION(GPP_USE_STD_COMPLEX "Enable std::complex" ON)

INCLUDE(GPP_Utilities)  # some macros
INCLUDE(GPP_Packages)   # creates interface libs and finds packages
                        # based on options

MESSAGE(STATUS "")
MESSAGE(STATUS "Configuration options:")
FOREACH(_OPT ${GPP_OPTIONS})
    # print ON/OFF
    MESSAGE(STATUS "    ${_OPT} = ${${_OPT}}")
    # set to zero or one for generator expressions
    IF(${_OPT})
        SET(${_OPT} 1)
    ELSE()
        SET(${_OPT} 0)
    ENDIF()
ENDFOREACH()
MESSAGE(STATUS "")

ADD_EXECUTABLE(gpp ${PROJECT_SOURCE_DIR}/GPP/gpp.cpp)

TARGET_LINK_LIBRARIES(gpp PRIVATE
    gpp-source
    gpp-compile-options
    gpp-std-complex
    $<IF:${GPP_USE_OPENMP},gpp-openmp,>
    $<IF:${GPP_USE_OPENMP_TARGET},gpp-openmp-target,>)

IF(GPP_BUILD_FF)
    ADD_EXECUTABLE(ff ${PROJECT_SOURCE_DIR}/FF/ff.cpp)

    TARGET_LINK_LIBRARIES(ff PRIVATE
        gpp-source
        gpp-compile-options
        gpp-std-complex
        $<IF:${GPP_USE_OPENMP},gpp-openmp,>
        $<IF:${GPP_USE_OPENMP_TARGET},gpp-openmp-target,>)
ENDIF()

IF(GPP_BUILD_TESTING)
    ADD_SUBDIRECTORY(test)
ENDIF()