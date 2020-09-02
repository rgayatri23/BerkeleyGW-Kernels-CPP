KOKKOS_PATH = $(HOME)/Kokkos/kokkos
CXXFLAGS = -O3

DEVICE=openmp

ifeq ($(DEVICE),openmp)
CXX=CC
CXXFLAGS += -qopenmp -xcore_avx2
KOKKOS_INC=$(KOKKOS_PATH)/install_OMP_icc/include
KOKKOS_LIB=$(KOKKOS_PATH)/install_OMP_icc/lib64
KOKKOS_DEVICES="OpenMP"
LDFLAGS = -qopenmp
else ifeq ($(DEVICE),cuda)
CXX=$(KOKKOS_PATH)/install_cuda/bin/nvcc_wrapper
CXXFLAGS += --expt-extended-lambda
KOKKOS_INC=$(KOKKOS_PATH)/install_cuda/include
KOKKOS_LIB=$(KOKKOS_PATH)/install_cuda/lib64
KOKKOS_DEVICES="Cuda"
LDFLAGS =
endif

gpp_obj = gpp.o

CXXFLAGS += -I$(CUDA_ROOT)/include -I$(KOKKOS_INC)
LDFLAGS += -L$(KOKKOS_LIB) -lkokkoscore

all: gpp
EXEC=gpp

gpp: $(gpp_obj)
	$(CXX) $(gpp_obj) -o $(EXEC).ex $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

gpp.cpp: ../common/Defines.h

clean:
	rm -f *.o *.lst $(EXEC).ex *.ptx *.cub
