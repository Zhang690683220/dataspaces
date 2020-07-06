DATASPACES_DIR=/home/bz186/myrepo/dataspaces/build/bin
DATASPACES_INC=-I/home/bz186/myrepo/dataspaces/build/include
DATASPACES_LIB=-L/home/bz186/myrepo/dataspaces/build/lib -ldspaces -ldscommon -ldart   -lm  -lrt   -libverbs  -lrdmacm  -L/home/bz186/zfp/lib/ -lzfp
ZFP_INC=-I/home/bz186/myrepo/dataspaces/include -I/home/bz186/zfp/include -I/home/bz186/myrepo/dataspaces

CC=mpicc

DEPS=put.c get.c test.c put_server.c get_server.c
OBJ=put.o get.o test.o put_server.o get_server.o

%.o: %.c $(DEPS)
	$(CC) -g -c -o $@ $< $(ZFP_INC) $(DATASPACES_INC) $(DATASPACES_LIB)


all: put get test put_server get_server

put: put.o
	$(CC) -o $@ $^ $(ZFP_INC) $(DATASPACES_INC) $(DATASPACES_LIB)


get: get.o
	$(CC) -o $@ $^ $(ZFP_INC) $(DATASPACES_INC) $(DATASPACES_LIB)

put_server: put_server.o
	$(CC) -o $@ $^ $(ZFP_INC) $(DATASPACES_INC) $(DATASPACES_LIB)


get_server: get_server.o
	$(CC) -o $@ $^ $(ZFP_INC) $(DATASPACES_INC) $(DATASPACES_LIB)

test: test.o
	$(CC) -g -o $@ $^ $(ZFP_INC) $(DATASPACES_INC) $(DATASPACES_LIB)

.PHONY: clean

clean:
	rm -f *.o put get test put_server get_server
