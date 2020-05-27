DATASPACES_DIR=/home/bz186/myrepo/dataspaces/build/bin
DATASPACES_INC=-I/home/bz186/myrepo/dataspaces/include
DATASPACES_LIB=-L/home/bz186/myrepo/dataspaces/build/lib -ldspaces -ldscommon -ldart   -lm  -lrt   -libverbs  -lrdmacm  -L/home/bz186/zfp/lib/ -lzfp
CC=mpicc

DEPS=put.c get.c
OBJ=put.o get.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(DATASPACES_INC) $(DATASPACES_LIB)


all: put get

put: put.o
	$(CC) -o $@ $^ $(DATASPACES_INC) $(DATASPACES_LIB)


get: get.o
	$(CC) -o $@ $^ $(DATASPACES_INC) $(DATASPACES_LIB)

.PHONY: clean

clean:
	rm -f *.o put get
