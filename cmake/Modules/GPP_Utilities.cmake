# include guard
include_guard(DIRECTORY)

cmake_policy(PUSH)
cmake_policy(SET CMP0054 NEW)

include(CMakeDependentOption)
include(CMakeParseArguments)


#------------------------------------------------------------------------------#
# macro for creating a library target
#
FUNCTION(GPP_CREATE_EXECUTABLE)
    # for include dirs, compile flags, definitions, etc. --> use INTERFACE libs
    # and add them to "LINK_LIBRARIES"
    # list of arguments taking multiple values
    set(multival_args
        HEADERS SOURCES PROPERTIES LINK_LIBRARIES INSTALL_DESTINATION)

    # parse args
    cmake_parse_arguments(EXE
        "INSTALL;EXCLUDE_FROM_ALL"   # options
        "TARGET_NAME;"               # single value args
        "${multival_args}"           # multiple value args
        ${ARGN})

    set(_EXCLUDE)
    if(EXE_EXCLUDE_FROM_ALL)
        set(_EXCLUDE EXCLUDE_FROM_ALL)
    endif()
    # create library
    add_executable(${EXE_TARGET_NAME} ${_EXCLUDE} ${EXE_SOURCES} ${EXE_HEADERS})

    # link library
    target_link_libraries(${EXE_TARGET_NAME} ${EXE_LINK_LIBRARIES})

    # target properties
    if(NOT "${EXE_PROPERTIES}" STREQUAL "")
        set_target_properties(${EXE_TARGET_NAME} PROPERTIES ${EXE_PROPERTIES})
    endif()

    if(EXE_INSTALL AND NOT EXE_INSTALL_DESTINATION)
        set(EXE_INSTALL_DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()

    # Install the exe
    if(EXE_INSTALL_DESTINATION)
        install(TARGETS ${EXE_TARGET_NAME} DESTINATION ${EXE_INSTALL_DESTINATION})
    endif()
ENDFUNCTION()

#------------------------------------------------------------------------------#
# macro add_googletest()
#
# Adds a unit test and links against googletest. Additional arguments are linked
# against the test.
#
FUNCTION(GPP_GOOGLE_TEST TEST_NAME)
    set(_OPTS )
    if(NOT GPP_BUILD_TESTING)
        set(_OPTS EXCLUDE_FROM_ALL)
    endif()
    if(NOT TARGET gpp-test)
        add_custom_target(gpp-test
            COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target test
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Running tests...")
    endif()
    include(GoogleTest)
    # list of arguments taking multiple values
    set(multival_args SOURCES PROPERTIES DEFINITIONS LINK_LIBRARIES COMMAND OPTIONS ENVIRONMENT)
    # parse args
    cmake_parse_arguments(TEST "DISCOVER_TESTS;ADD_TESTS" "" "${multival_args}" ${ARGN})

    if(NOT TARGET google-test-debug-options)
        add_library(google-test-debug-options INTERFACE)
        target_compile_definitions(google-test-debug-options INTERFACE
            $<$<CONFIG:Debug>:DEBUG> GPP_TESTING)
    endif()
    list(APPEND TEST_LINK_LIBRARIES google-test-debug-options)

    GPP_CREATE_EXECUTABLE(${_OPTS}
        TARGET_NAME     ${TEST_NAME}
        OUTPUT_NAME     ${TEST_NAME}
        SOURCES         ${TEST_SOURCES}
        LINK_LIBRARIES  gpp-google-test ${TEST_LINK_LIBRARIES}
        PROPERTIES      "${TEST_PROPERTIES}")

    if(TEST_DEFINITIONS)
        target_compile_definitions(${TEST_NAME} PRIVATE ${TEST_DEFINITIONS})
    endif()

    # always add as a dependency
    add_dependencies(gpp-test ${TEST_NAME})

    if("${TEST_COMMAND}" STREQUAL "")
        set(TEST_COMMAND $<TARGET_FILE:${TEST_NAME}>)
    endif()

    if(TEST_DISCOVER_TESTS)
        GTEST_DISCOVER_TESTS(${TEST_NAME}
            ${TEST_OPTIONS})
    elseif(TEST_ADD_TESTS)
        GTEST_ADD_TESTS(TARGET ${TEST_NAME}
            ${TEST_OPTIONS})
    else()
        ADD_TEST(
            NAME                ${TEST_NAME}
            COMMAND             ${TEST_COMMAND}
            WORKING_DIRECTORY   ${CMAKE_CURRENT_LIST_DIR}
            ${TEST_OPTIONS})
        SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES ENVIRONMENT "${TEST_ENVIRONMENT}")
    endif()

ENDFUNCTION()

#----------------------------------------------------------------------------------------#
# macro gpp_checkout_git_submodule()
#
#   Run "git submodule update" if a file in a submodule does not exist
#
#   ARGS:
#       RECURSIVE (option) -- add "--recursive" flag
#       RELATIVE_PATH (one value) -- typically the relative path to submodule
#                                    from PROJECT_SOURCE_DIR
#       WORKING_DIRECTORY (one value) -- (default: PROJECT_SOURCE_DIR)
#       TEST_FILE (one value) -- file to check for (default: CMakeLists.txt)
#       ADDITIONAL_CMDS (many value) -- any addition commands to pass
#
FUNCTION(GPP_CHECKOUT_GIT_SUBMODULE)
    # parse args
    cmake_parse_arguments(
        CHECKOUT
        "RECURSIVE"
        "RELATIVE_PATH;WORKING_DIRECTORY;TEST_FILE;REPO_URL;REPO_BRANCH"
        "ADDITIONAL_CMDS"
        ${ARGN})

    if(NOT CHECKOUT_WORKING_DIRECTORY)
        set(CHECKOUT_WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    endif()

    if(NOT CHECKOUT_TEST_FILE)
        set(CHECKOUT_TEST_FILE "CMakeLists.txt")
    endif()

    # default assumption
    if(NOT CHECKOUT_REPO_BRANCH)
        set(CHECKOUT_REPO_BRANCH "master")
    endif()

    find_package(Git)
    set(_DIR "${CHECKOUT_WORKING_DIRECTORY}/${CHECKOUT_RELATIVE_PATH}")
    # ensure the (possibly empty) directory exists
    if(NOT EXISTS "${_DIR}")
        if(NOT CHECKOUT_REPO_URL)
            message(FATAL_ERROR "submodule directory does not exist")
        endif()
    endif()

    # if this file exists --> project has been checked out
    # if not exists --> not been checked out
    set(_TEST_FILE "${_DIR}/${CHECKOUT_TEST_FILE}")
    # assuming a .gitmodules file exists
    set(_SUBMODULE "${PROJECT_SOURCE_DIR}/.gitmodules")

    set(_TEST_FILE_EXISTS OFF)
    if(EXISTS "${_TEST_FILE}" AND NOT IS_DIRECTORY "${_TEST_FILE}")
        set(_TEST_FILE_EXISTS ON)
    endif()

    if(_TEST_FILE_EXISTS)
        return()
    endif()

    find_package(Git REQUIRED)

    set(_SUBMODULE_EXISTS OFF)
    if(EXISTS "${_SUBMODULE}" AND NOT IS_DIRECTORY "${_SUBMODULE}")
        set(_SUBMODULE_EXISTS ON)
    endif()

    set(_HAS_REPO_URL OFF)
    if(NOT "${CHECKOUT_REPO_URL}" STREQUAL "")
        set(_HAS_REPO_URL ON)
    endif()

    # if the module has not been checked out
    if(NOT _TEST_FILE_EXISTS AND _SUBMODULE_EXISTS)
        # perform the checkout
        execute_process(
            COMMAND
                ${GIT_EXECUTABLE} submodule update --init ${_RECURSE}
                    ${CHECKOUT_ADDITIONAL_CMDS} ${CHECKOUT_RELATIVE_PATH}
            WORKING_DIRECTORY
                ${CHECKOUT_WORKING_DIRECTORY}
            RESULT_VARIABLE RET)

        # check the return code
        if(RET GREATER 0)
            set(_CMD "${GIT_EXECUTABLE} submodule update --init ${_RECURSE}
                ${CHECKOUT_ADDITIONAL_CMDS} ${CHECKOUT_RELATIVE_PATH}")
            message(STATUS "function(gpp_checkout_git_submodule) failed.")
            message(FATAL_ERROR "Command: \"${_CMD}\"")
        else()
            set(_TEST_FILE_EXISTS ON)
        endif()
    endif()

    if(NOT _TEST_FILE_EXISTS AND _HAS_REPO_URL)
        message(STATUS "Checking out '${CHECKOUT_REPO_URL}' @ '${CHECKOUT_REPO_BRANCH}'...")

        # remove the existing directory
        if(EXISTS "${_DIR}")
            execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${_DIR})
        endif()

        # perform the checkout
        execute_process(
            COMMAND
                ${GIT_EXECUTABLE} clone -b ${CHECKOUT_REPO_BRANCH}
                    ${CHECKOUT_ADDITIONAL_CMDS}
                    ${CHECKOUT_REPO_URL} ${CHECKOUT_RELATIVE_PATH}
            WORKING_DIRECTORY
                ${CHECKOUT_WORKING_DIRECTORY}
            RESULT_VARIABLE RET)

        # perform the submodule update
        if(CHECKOUT_RECURSIVE AND EXISTS "${_DIR}" AND IS_DIRECTORY "${_DIR}")
            execute_process(
                COMMAND
                    ${GIT_EXECUTABLE} submodule update --init ${_RECURSE}
                WORKING_DIRECTORY
                    ${_DIR}
                RESULT_VARIABLE RET)
        endif()

        # check the return code
        if(RET GREATER 0)
            set(_CMD "${GIT_EXECUTABLE} clone -b ${CHECKOUT_REPO_BRANCH}
                ${CHECKOUT_ADDITIONAL_CMDS} ${CHECKOUT_REPO_URL} ${CHECKOUT_RELATIVE_PATH}")
            message(STATUS "function(gpp_checkout_git_submodule) failed.")
            message(FATAL_ERROR "Command: \"${_CMD}\"")
        else()
            set(_TEST_FILE_EXISTS ON)
        endif()
    endif()

    if(NOT EXISTS "${_TEST_FILE}" OR NOT _TEST_FILE_EXISTS)
        message(FATAL_ERROR "Error checking out submodule: '${CHECKOUT_RELATIVE_PATH}' to '${_DIR}'")
    endif()

ENDFUNCTION()


#----------------------------------------------------------------------------------------#
# macro to add an interface lib
#
FUNCTION(GPP_ADD_INTERFACE_LIBRARY _TARGET)
    add_library(${_TARGET} INTERFACE)
    add_library(${PROJECT_NAME}::${_TARGET} ALIAS ${_TARGET})
ENDFUNCTION()


#----------------------------------------------------------------------------------------#
# macro to build a library of type: shared, static, object
#
macro(GPP_BUILD_LIBRARY)

    # options
    set(_options    PIC NO_CACHE_LIST)
    # single-value
    set(_onevalue   TYPE
                    OUTPUT_NAME
                    TARGET_NAME
                    OUTPUT_DIR
                    LANGUAGE
                    LINKER_LANGUAGE)
    # multi-value
    set(_multival   SOURCES
                    LINK_LIBRARIES
                    COMPILE_DEFINITIONS
                    INCLUDE_DIRECTORIES
                    C_COMPILE_OPTIONS
                    CXX_COMPILE_OPTIONS
                    CUDA_COMPILE_OPTIONS
                    LINK_OPTIONS
                    EXTRA_PROPERTIES)

    cmake_parse_arguments(
        LIBRARY "${_options}" "${_onevalue}" "${_multival}" ${ARGN})

    if("${LIBRARY_LANGUAGE}" STREQUAL "")
        set(LIBRARY_LANGUAGE CXX)
    endif()

    if("${LIBRARY_LINKER_LANGUAGE}" STREQUAL "")
        set(LIBRARY_LINKER_LANGUAGE CXX)
    endif()

    if("${LIBRARY_OUTPUT_DIR}" STREQUAL "")
        set(LIBRARY_OUTPUT_DIR ${PROJECT_BINARY_DIR})
    endif()

    if(NOT WIN32 AND NOT XCODE)
        list(APPEND LIBRARY_EXTRA_PROPERTIES
            VERSION                     ${PROJECT_VERSION}
            SOVERSION                   ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
    endif()

    if(NOT WIN32)
        set(LIB_PREFIX )
        list(APPEND LIBRARY_EXTRA_PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY    ${LIBRARY_OUTPUT_DIR}
            ARCHIVE_OUTPUT_DIRECTORY    ${LIBRARY_OUTPUT_DIR}
            RUNTIME_OUTPUT_DIRECTORY    ${LIBRARY_OUTPUT_DIR})
    else()
        set(LIB_PREFIX lib)
    endif()

    # add the library or sources
    if(NOT TARGET ${LIBRARY_TARGET_NAME})
        add_library(${LIBRARY_TARGET_NAME} ${LIBRARY_TYPE} ${LIBRARY_SOURCES})
        add_library(${PROJECT_NAME}::${LIBRARY_TARGET_NAME} ALIAS ${LIBRARY_TARGET_NAME})
    else()
        target_sources(${LIBRARY_TARGET_NAME} PRIVATE ${LIBRARY_SOURCES})
    endif()

    # append include directories
    target_include_directories(${LIBRARY_TARGET_NAME}
        PUBLIC ${LIBRARY_INCLUDE_DIRECTORIES})

    # compile definitions
    target_compile_definitions(${LIBRARY_TARGET_NAME}
        PUBLIC ${LIBRARY_COMPILE_DEFINITIONS})

    # compile flags
    target_compile_options(${LIBRARY_TARGET_NAME}
        PRIVATE
            $<$<COMPILE_LANGUAGE:C>:${LIBRARY_C_COMPILE_OPTIONS}>
            $<$<COMPILE_LANGUAGE:CXX>:${LIBRARY_CXX_COMPILE_OPTIONS}>)

    # cuda flags
    get_property(LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
    if(CMAKE_CUDA_COMPILER AND "CUDA" IN_LIST LANGUAGES)
        target_compile_options(${LIBRARY_TARGET_NAME}
            PRIVATE
                $<$<COMPILE_LANGUAGE:CUDA>:${LIBRARY_CUDA_COMPILE_OPTIONS}>)
    endif()

    # link options
    if(NOT CMAKE_VERSION VERSION_LESS 3.13)
        target_link_options(${LIBRARY_TARGET_NAME} PUBLIC ${LIBRARY_LINK_OPTIONS})
    elseif(NOT "${LIBRARY_LINK_OPTIONS}" STREQUAL "")
        list(APPEND LIBRARY_EXTRA_PROPERTIES LINK_OPTIONS ${LIBRARY_LINK_OPTIONS})
    endif()

    # link libraries
    target_link_libraries(${LIBRARY_TARGET_NAME}
        PUBLIC ${LIBRARY_LINK_LIBRARIES})

    # other properties
    set_target_properties(
        ${LIBRARY_TARGET_NAME}      PROPERTIES
        OUTPUT_NAME                 ${LIB_PREFIX}${LIBRARY_OUTPUT_NAME}
        LANGUAGE                    ${LIBRARY_LANGUAGE}
        LINKER_LANGUAGE             ${LIBRARY_LINKER_LANGUAGE}
        POSITION_INDEPENDENT_CODE   ${LIBRARY_PIC}
        ${LIBRARY_EXTRA_PROPERTIES})

    if(NOT LIBRARY_NO_CACHE_LIST)
        # add to cached list of compiled libraries
        set(COMPILED_TYPES "SHARED" "STATIC" "MODULE")
        if("${LIBRARY_TYPE}" IN_LIST COMPILED_TYPES)
            cache_list(APPEND ${PROJECT_NAME_UC}_COMPILED_LIBRARIES ${LIBRARY_TARGET_NAME})
        endif()
    endif()
    unset(COMPILED_TYPES)

endmacro()


#----------------------------------------------------------------------------------------#
# require variable
#
function(GPP_CHECK_REQUIRED VAR)
    if(NOT DEFINED ${VAR} OR "${${VAR}}" STREQUAL "")
        message(FATAL_ERROR "Variable '${VAR}' must be defined and not empty")
    endif()
endfunction()


#----------------------------------------------------------------------------------------#

cmake_policy(POP)
