include ../Makefile.inc


ff_obj = ff.o

CXXFLAGS += -I$(CUDA_ROOT)/include

std_complex=y
ifeq ($(std_complex),y)
	CXXFLAGS += -Dstd_complex
endif

all: ff
EXEC=ff

ff: $(ff_obj)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(ff_obj) -o $(EXEC).ex

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

ff.cpp: ../ComplexClass/CustomComplex.h

clean:
	rm -f *.o *.lst $(EXEC).ex *.ptx *.cub
