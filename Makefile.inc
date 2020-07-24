NVCC:=nvcc
OPENMP:=n
OPENMP_TARGET:=n
OPENACC:=n
GPU:=n
GPU_OPENACC_UM:=n
GPU_CUDA_UM:=n
BENCHMARK:=n
CUDA_PATH?=$(CUDA_ROOT)
USE_OPENMP_SCHEDULE:=n
USE_INLINE_MATH:=n
TIMEMORY_PROFILE:=n

ifeq ($(findstring lang,$(COMP)),lang)
     CC=clang
     CXX=clang++
     FC=flang
     F77=flang
     COMP_CHOICE=clang
else ifeq ($(findstring xl,$(COMP)),xl)
     CC=xlc_r
     CXX=xlC_r
     FC=xlf90_r
     F77=xlf_r
     COMP_CHOICE=xl
else ifeq ($(findstring pg,$(COMP)),pg)
     CC=pgcc
     CXX=pgc++
     FC=pgf90
     F77=pgf77
     COMP_CHOICE=pgi
else ifeq ($(findstring cray,$(COMP)),cray)
     CC=CC
     CXX=CC
     COMP_CHOICE=cray
else ifeq ($(findstring ic,$(COMP)),ic)
     CC=icc
     CXX=icpc
     FC=ifort
     F77=ifort
     COMP_CHOICE=intel
else ifeq ($(findstring gc,$(COMP)),gc)
     CC=gcc
     CXX=g++
     FC=gfortran
     F77=gfortran
     COMP_CHOICE=gcc
else ifeq ($(findstring gf,$(COMP)),gf)
     CC=gcc
     CXX=g++
     FC=gfortran
     F77=gfortran
     COMP_CHOICE=gcc
else ifeq ($(findstring g+,$(COMP)),g+)
     CC=gcc
     CXX=g++
     FC=gfortran
     F77=gfortran
     COMP_CHOICE=gcc
else
     $(warning Warning COMP is not set)
endif

# We can use OpenMP target code on the host and the GPU
ifeq ($(OPENMP_TARGET),y)
  override CPPFLAGS += -DOPENMP_TARGET
endif

# We can use OpenMP target code on the host and the GPU
ifeq ($(BENCHMARK),y)
  override CPPFLAGS += -DBENCHMARK
endif

# Allow the user to override CPPFLAGS only, e.g.
# make CC=gcc CPPFLAGS=-DSTREAM_ARRAY_SIZE=100000
ifeq ($(GPU),y)
  ifeq ($(OPENMP_TARGET),y)
    # We can only use OpenMP target code for the GPU
    override CPPFLAGS += -DGPU
  endif
  ifeq ($(OPENACC),y)
    override CPPFLAGS += -DOPENACC -DGPU
  endif
  ifeq ($(GPU_OPENACC_UM),y)
    override CPPFLAGS += -DGPU_OPENACC_UM -DGPU
  endif
  ifeq ($(GPU_CUDA_UM),y)
    NVCCFLAGS = -ccbin=$(CXX)
    override CPPFLAGS += -DGPU_CUDA_UM -DGPU
  endif
endif

ifeq ($(USE_OPENMP_SCHEDULE),y)
  override CPPFLAGS += -DUSE_OPENMP_SCHEDULE
endif
ifeq ($(USE_INLINE_MATH),y)
  override CPPFLAGS += -DUSE_INLINE_MATH
endif

CDLDFLAGS=-lstdc++ -L$(dir $(shell which nvcc))../lib64 -lcudart

#
# IBM XL COMPILER - xlc_r
#
ifeq ($(COMP_CHOICE),xl)
  CFLAGS += -Ofast

  ifeq ($(OPENMP),y)
    CFLAGS += -qsmp=omp:noauto
  else ifeq ($(OPENMP_TARGET),y)
    CFLAGS += -qsmp=omp:noauto
  endif

  # Compiler flags for different platforms
  ifeq ($(LMOD_SYSTEM_NAME),summitdev)
    CFLAGS += -qarch=pwr8 -qtune=pwr8
  else ifeq ($(LMOD_SYSTEM_NAME),summit)
    CFLAGS += -qarch=pwr9 -qtune=pwr9
  endif

  # Explicitly turn off with -qnooffload. This seems to execute with 1 thread though
  ifeq ($(GPU),y)
    CFLAGS += -qoffload
    ifeq ($(LMOD_SYSTEM_NAME),summitdev)
       CFLAGS += -qtgtarch=sm_60
    else ifeq ($(LMOD_SYSTEM_NAME),summit)
       CFLAGS += -qtgtarch=sm_70
    endif
  endif

  ifeq ($(VERBOSE),y)
    CFLAGS += -qreport
    ifeq ($(GPU),y)
      CFLAGS += -Xptxas -v
    endif
  endif

  FCFLAGS := $(CFLAGS)
  CFLAGS += -fstrict-aliasing -Werror -Wall -Wno-unused-variable
  CXXFLAGS += $(CFLAGS) -std=c++11

#
# LLVM COMPILER - clang
#
else ifeq ($(COMP_CHOICE),clang)

  # -O3 and -Ofast cause a deadlock with Clang-7.0.0
  CFLAGS += -O2 -ffast-math -ffp-contract=fast -fstrict-aliasing -Wall -Wno-unused-variable
  LDFLAGS += -lm

  ifeq ($(OPENMP),y)
    CFLAGS += -fopenmp
  else ifeq ($(OPENMP_TARGET),y)
    CFLAGS += -fopenmp
  endif

  # Compiler flags for different platforms
  ifeq ($(LMOD_SYSTEM_NAME),summitdev)
    CFLAGS += --target="powerpc64le-unknown-linux-gnu"
  else ifeq ($(LMOD_SYSTEM_NAME),summit)
    CFLAGS += --target="powerpc64le-unknown-linux-gnu"
  endif

  ifeq ($(GPU),y)
    # clang-ykt option: -fopenmp-nonaliased-maps
    CFLAGS += -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=${CUDA_PATH} -ffp-contract=fast

    # Arch specific flags: llvm/tools/clang/test/Driver/openmp-offload-gpu.c
    ifeq ($(LMOD_SYSTEM_NAME),summitdev)
      CFLAGS += -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_60
    else ifeq ($(LMOD_SYSTEM_NAME),summit)
      CFLAGS += -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_70
    else
      # Macro defines used for bug https://bugs.llvm.org/show_bug.cgi?id=38464 - only needed on x86
      CFLAGS += -D__NO_MATH_INLINES -U__SSE2_MATH__ -U__SSE_MATH__
    endif
  endif

  ifeq ($(VERBOSE),y)
    # clang-ykt option: -Rpass-analysis (shows "Combined construct detected.")
    # -Rpass=.* -Rpass=inline
    CFLAGS += -Rpass-analysis
    ifeq ($(GPU),y)
      CFLAGS += -Xcuda-ptxas -v -Wno-unused-command-line-argument
    endif
  endif

  CXXFLAGS += $(CFLAGS) -std=c++11
  FCFLAGS += $(CFLAGS)
#
# Cray COMPILER - CC
#
else ifeq ($(COMP_CHOICE),cray)
  ifeq ($(OPENMP_TARGET),y)
      ifeq ($(GPU),y)
            CFLAGS += -hunroll1 -O2 -hvector2,scalar3 -h flex_mp=tolerant -hlist=m
            CFLAGS += -hnomessage=7275 -lcrypto
      endif
  endif
  CXXFLAGS += $(CFLAGS)

#
# PGI COMPILER - pgcc
#
else ifeq ($(COMP_CHOICE),pgi)

  CFLAGS += -fast
  ifeq ($(OPENMP),y)
    CFLAGS += -mp
  else ifeq ($(OPENMP_TARGET),y)
    CFLAGS += -mp
  else ifeq ($(OPENACC),y)
    CFLAGS += -acc #-Minfo=accel #-ta=tesla:managed
  endif

  ifeq ($(GPU),y)
    ifeq ($(OPENACC),y)
      tgtflags=tesla:fastmath

      ifeq ($(GPU_OPENACC_UM),y)
        tgtflags:=$(tgtflags),managed
      endif
      ifeq ($(VERBOSE),y)
        tgtflags:=$(tgtflags),ptxinfo
      endif
      ifeq ($(LMOD_SYSTEM_NAME),summitdev)
        tgtflags:=$(tgtflags),cc60,cuda9.0
      else ifeq ($(LMOD_SYSTEM_NAME),summit)
        tgtflags:=$(tgtflags),cc70,cuda9.2
      else ifeq ($(HOSTNAME),psgcluster)
        tgtflags:=$(tgtflags),cc70,cuda9.1
      else ifeq ($(SLURM_CLUSTER_NAME),escori)
        tgtflags:=$(tgtflags),cc70,cuda9.2
      endif

      CFLAGS += -ta=$(tgtflags)
    endif
  else
    ifeq ($(OPENACC),y)
      CFLAGS += -ta=multicore
    endif
  endif

  ifeq ($(VERBOSE),y)
    ifeq ($(OPENMP),y)
      CFLAGS += -Minfo=mp
    else ifeq ($(OPENACC),y)
      CFLAGS += -Minfo=accel
    endif
  endif

  CXXFLAGS += $(CFLAGS) -std=c++11
  FCFLAGS += $(CFLAGS)


#
# INTEL COMPILER - icc
#
else ifeq ($(COMP_CHOICE),intel)

  # KNL: -xMIC-AVX512
  CFLAGS += -Ofast -qopt-report=5
  ifeq ($(OPENMP),y)
    CFLAGS += -qopenmp
    ifeq ($(CRAY_CPU_TARGET),haswell)
        CFLAGS += -xcore-avx2
    else ifeq ($(CRAY_CPU_TARGET),mic-knl)
        CFLAGS += -xmic-avx512
    endif
  endif
  ifeq ($(OPENMP_TARGET),y)
    CFLAGS += -qopenmp-offload=host -qopt-report=5
      LDFLAGS += -qopenmp-offload=mic
  endif

  CXXFLAGS += $(CFLAGS) -std=c++11
  FCFLAGS += $(CFLAGS)

#
# GNU COMPILER - gcc
#
else ifeq ($(COMP_CHOICE),gcc)

  CFLAGS += -Ofast -fstrict-aliasing
  ifeq ($(OPENMP),y)
    CFLAGS += -fopenmp
  else ifeq ($(OPENMP_TARGET),y)
    CFLAGS += -fopenmp
  else ifeq ($(OPENACC),y)
    CFLAGS += -fopenacc
  endif

  ifeq ($(GPU),y)
    CFLAGS += -foffload=nvptx-none="-Ofast -lm -misa=sm_35"
  else
    ifeq ($(OPENMP_TARGET),y)
      CFLAGS += -foffload=disable
    endif
  endif


  # Compiler flags for different platforms
  ifeq ($(LMOD_SYSTEM_NAME),summitdev)
    CFLAGS += -mcpu=power8 -mtune=power8
  else ifeq ($(LMOD_SYSTEM_NAME),summit)
    CFLAGS += -mcpu=power9 -mtune=power9
  endif

  ifeq ($(VERBOSE),y)
    ifeq ($(GPU),y)
      CFLAGS += -fopt-info-all-omp -fopt-info-vec
    endif
  endif

  CXXFLAGS += $(CFLAGS) -std=c++11
  FCFLAGS += $(CFLAGS)
endif

#Set TiMemory
ifeq ($(TIMEMORY_PROFILE),y)
  CXXFLAGS += -I/project/projectdirs/m1759/timemory/corigpu/include
  CXXFLAGS += -DUSE_TIMEMORY
  CXXFLAGS += -I/opt/cray/pe/papi/5.6.0.6/include
  CXXFLAGS += -DTIMEMORY_USE_PAPI
endif
