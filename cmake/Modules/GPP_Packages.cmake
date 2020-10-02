# include guard
include_guard(DIRECTORY)

enable_testing()

##########################################################################################
#
#                       External Packages are found here
#
##########################################################################################

gpp_add_interface_library(gpp-source
    "Provides the source code and the include directories without any flags")
gpp_add_interface_library(gpp-compile-options
    "Provides general compiler options")
gpp_add_interface_library(gpp-google-test
    "Provides general compiler options")
gpp_add_interface_library(gpp-timemory
    "Provides timemory")
gpp_add_interface_library(gpp-std-complex
    "Provies OpenMP target")

if(GPP_USE_OPENMP)
    gpp_add_interface_library(gpp-openmp
        "Provies OpenMP")
    find_package(OpenMP REQUIRED)
endif()

if(GPP_USE_OPENMP_TARGET)
    gpp_add_interface_library(gpp-openmp-target
        "Provies OpenMP Target")
    find_package(OpenMP REQUIRED)
    find_package(CUDA REQUIRED)
endif()

if(GPP_USE_STD_COMPLEX)
    target_compile_definitions(gpp-std-complex INTERFACE std_complex)
endif()

#----------------------------------------------------------------------------------------#
#
#                            Source
#
#----------------------------------------------------------------------------------------#

target_include_directories(gpp-source INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ComplexClass>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/arrayMD>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/GPP>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/common>
)

file(GLOB GPP_GPP_SOURCE
    ${PROJECT_SOURCE_DIR}/ComplexClass/*.h
    ${PROJECT_SOURCE_DIR}/arrayMD/*.h
    ${PROJECT_SOURCE_DIR}/common/*.h
)

target_sources(gpp-source INTERFACE ${GPP_SOURCE})

target_compile_options(gpp-compile-options INTERFACE
    -W
    -Wall
    -Wno-unused-variable
    -Wno-unused-parameter
    -Wno-sign-compare
    -Wno-unknown-pragmas)

#----------------------------------------------------------------------------------------#
#
#                           Google Test
#
#----------------------------------------------------------------------------------------#

gpp_checkout_git_submodule(RECURSIVE
    RELATIVE_PATH external/google-test
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    REPO_URL https://github.com/google/googletest.git
    REPO_BRANCH master)

# add google-test
set(INSTALL_GTEST OFF CACHE BOOL "Install gtest")
set(BUILD_GMOCK ON CACHE BOOL "Build gmock")
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON CACHE BOOL "Enable MACOS_RPATH on targets to suppress warnings")
    mark_as_advanced(CMAKE_MACOSX_RPATH)
endif()
add_subdirectory(${PROJECT_SOURCE_DIR}/external/google-test)
target_link_libraries(gpp-google-test INTERFACE gtest gmock gtest_main)
target_include_directories(gpp-google-test SYSTEM INTERFACE
    ${PROJECT_SOURCE_DIR}/google-test/googletest/include
    ${PROJECT_SOURCE_DIR}/google-test/googlemock/include)


#----------------------------------------------------------------------------------------#
#
#                               Timemory
#
#----------------------------------------------------------------------------------------#

if(GPP_USE_TIMEMORY)
    find_package(timemory REQUIRED COMPONENTS headers arch cxx)
    target_link_libraries(gpp-timemory INTERFACE timemory)
endif()


#----------------------------------------------------------------------------------------#
#
#                               OpenMP
#
#----------------------------------------------------------------------------------------#

if(GPP_USE_OPENMP)
    target_compile_definitions(gpp-openmp INTERFACE GPP_USE_OPENMP)
    target_link_libraries(gpp-openmp INTERFACE OpenMP::OpenMP_CXX)
endif()

if(GPP_USE_OPENMP_TARGET)
    target_compile_definitions(gpp-openmp INTERFACE
        GPP_USE_OPENMP_TARGET
        OPENMP_TARGET)
    target_link_libraries(gpp-openmp-target INTERFACE OpenMP::OpenMP_CXX)

    if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_include_directories(gpp-openmp-target INTERFACE
            ${CUDA_TOOLKIT_ROOT_DIR}/include)
        target_compile_options(gpp-openmp-target INTERFACE
            -fopenmp-targets=nvptx64-nvidia-cuda
            --cuda-path=${CUDA_TOOLKIT_ROOT_DIR}
            -fopenmp-cuda-mode
            -ffast-math
            -ffp-contract=fast
            -fstrict-aliasing)
        target_link_options(gpp-openmp-target INTERFACE
            -fopenmp-targets=nvptx64-nvidia-cuda
            -fopenmp=libiomp5)
    elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "XL")
        target_compile_options(gpp-openmp-target INTERFACE
            -qsmp=omp:noauto
            -qoffload)
        target_link_options(gpp-openmp-target INTERFACE
            -qsmp=omp:noauto
            -qoffload)
    else()
        message(FATAL_ERROR "Invalid compiler for openmp target")
    endif()
endif()


