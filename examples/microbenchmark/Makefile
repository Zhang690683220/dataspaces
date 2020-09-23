DATASPACES_DIR=/home/bz186/dataspaces/build_rdma/bin
DATASPACES_INC=-I/home/bz186/dataspaces/build_rdma/include
DATASPACES_LIB=-L/home/bz186/dataspaces/build_rdma/lib -ldspaces -ldscommon -ldart   -lm  -lrt   -libverbs  -lrdmacm
CXXFLAGS= --std=c++11
CXX=mpicxx


all: stdio dataspaces

stdio: stdio.cpp Timer.hpp
	$(CXX) -o $@ $< $(CXXFLAGS)

dataspaces: dataspaces.cpp Timer.hpp
	$(CXX) -o $@ $^ $(DATASPACES_INC) $(DATASPACES_LIB) $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f *.o stdio dataspaces

