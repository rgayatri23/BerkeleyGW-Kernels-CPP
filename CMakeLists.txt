cmake_minimum_required (VERSION 3.12)
project(BerkeleyGW-kernels
  LANGUAGES C CXX
  )

SET(CMAKE_CXX_STANDARD 11)
# don't allow in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(STATUS "Warning! Building from the source directory is not allow")
  message(STATUS "Remove 'CMakeCache.txt' and 'CMakeFiles' and build from a separate directory")
  message(ERROR "In-source build")
endif()

SET(GPPTarget gpp)
FILE(GLOB gpp-sources
    ${CMAKE_SOURCE_DIR}/GPP/gpp.cpp
    ${CMAKE_SOURCE_DIR}/ComplexClass/CustomComplex.h
    ${CMAKE_SOURCE_DIR}/arrayMD/arrayMD.h
    ${CMAKE_SOURCE_DIR}/common/Defines.h
  )

MESSAGE(STATUS "CMAKE_SOURCE_DIR = " ${CMAKE_SOURCE_DIR})

SET(BUILD_FF OFF CACHE BOOL "BUILD_FF Option")
MESSAGE(STATUS "BUILDING FF " ${BUILD_FF})

ADD_EXECUTABLE(
    ${GPPTarget}
    ${gpp-sources}
  )

if(${OPENMP})
  find_package(OpenMP REQUIRED)
  target_link_libraries(${GPPTarget} PUBLIC OpenMP::OpenMP_CXX)
endif()

if(${OPENMP_TARGET})
  find_package(CUDA REQUIRED)
  find_package(OpenMP REQUIRED)
  target_compile_definitions(${GPPTarget} PRIVATE OPENMP_TARGET)
  target_link_libraries(${GPPTarget} PUBLIC OpenMP::OpenMP_CXX)

  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  target_compile_options(${GPPTarget} PRIVATE -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=${CUDA_TOOLKIT_ROOT_DIR} -I${CUDA_TOOLKIT_ROOT_DIR}/include -ffp-contract=fast)
  target_compile_options(${GPPTarget} PRIVATE -fopenmp-cuda-mode)
  target_compile_options(${GPPTarget} PRIVATE -ffast-math -fstrict-aliasing -Wall -Wno-unused-variable -Wall)
  target_link_options(${GPPTarget} PRIVATE -fopenmp-targets=nvptx64-nvidia-cuda -fopenmp=libiomp5)
  endif()

  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "XL")
  target_compile_options(${GPPTarget} PRIVATE -qsmp=omp:noauto -qoffload)
  target_link_options(${GPPTarget} PRIVATE -qsmp=omp:noauto -qoffload)
  endif()
  target_link_libraries(${GPPTarget} PUBLIC OpenMP::OpenMP_CXX)
endif()

if(${BUILD_FF})
  SET(FFTarget ff)
  FILE(GLOB ff-sources
      ${CMAKE_SOURCE_DIR}/FF/ff.cpp
      ${CMAKE_SOURCE_DIR}/ComplexClass/CustomComplex.h
      ${CMAKE_SOURCE_DIR}/arrayMD/*.h
    )

  ADD_EXECUTABLE(
      ${FFTarget}
      ${ff-sources}
    )

if(${OPENMP})
  target_link_libraries(${FFTarget} PUBLIC OpenMP::OpenMP_CXX)
endif()

if(${OPENMP_TARGET})
  target_compile_definitions(${FFTarget} PRIVATE OPENMP_TARGET)

  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(${FFTarget} PRIVATE -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=${CUDA_TOOLKIT_ROOT_DIR} -I${CUDA_TOOLKIT_ROOT_DIR}/include -ffp-contract=fast)
    target_compile_options(${FFTarget} PRIVATE -fopenmp-cuda-mode)
    target_compile_options(${FFTarget} PRIVATE -ffast-math -fstrict-aliasing -Wall -Wno-unused-variable -Wall)
    target_link_options(${FFTarget} PRIVATE -fopenmp-targets=nvptx64-nvidia-cuda -fopenmp=libiomp5)
  endif()

  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "XL")
    target_compile_options(${FFTarget} PRIVATE -qsmp=omp:noauto -qoffload)
    target_link_options(${FFTarget} PRIVATE -qsmp=omp:noauto -qoffload)
  endif()
  target_link_libraries(${FFTarget} PUBLIC OpenMP::OpenMP_CXX)
endif()

endif()
