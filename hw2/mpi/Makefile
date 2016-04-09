CC=mpicc
CFLAGS=-L. -L$CILKHOME/lib64 -Wl,-rpath=. -lmul

all: libmul.so
	$(CC) $(CFLAGS) -o csort distributedsort.cpp
libmul.so:
	icc -m64 -fPIC -shared -o libmul.so parallelmerge.cpp
clean:
	rm -rf multiply-mpi-cilk core*
