include ../Makefile.inc


gpp_obj = gpp.o

CXXFLAGS += -I$(CUDA_ROOT)/include

all: gpp
EXEC=gpp

#LDFLAGS += -L/opt/cray/pe/papi/5.6.0.6/lib -lpapi

gpp: $(gpp_obj)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(gpp_obj) -o $(EXEC).ex

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

gpp.cpp: ../ComplexClass/CustomComplex.h

clean:
	rm -f *.o *.lst $(EXEC).ex *.ptx *.cub
