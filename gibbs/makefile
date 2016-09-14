JAVA_HOME=/home/yulele/software/jdk1.7.0_79

IDIR =./include
CC=g++
CFLAGS=-I$(IDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -std=c++11 -O3 -fPIC -march=native -g

ODIR=obj
LDIR=lib
BDIR=bin
SDIR=src

LIBS=-fopenmp

_DEPS = utils.h csr.h xorshift.h warplda.h cgslda.h traverse_map.h lda_jni_Utils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_WARPLDA_OBJ = utils.o csr.o warplda.o 
WARPLDA_OBJ = $(patsubst %,$(ODIR)/%,$(_WARPLDA_OBJ))

_CGSLDA_OBJ = utils.o csr.o cgslda.o 
CGSLDA_OBJ = $(patsubst %,$(ODIR)/%,$(_CGSLDA_OBJ))

_FTREELDA_OBJ = utils.o csr.o ftreelda.o 
FTREELDA_OBJ = $(patsubst %,$(ODIR)/%,$(_FTREELDA_OBJ))

_OPENMP_OBJ = openmp.o
OPENMP_OBJ = $(patsubst %,$(ODIR)/%,$(_OPENMP_OBJ))

_WARPJNI_OBJ = utils.o csr.o warp_jni.o
WARPJNI_OBJ = $(patsubst %,$(ODIR)/%,$(_WARPJNI_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)


all: warplda cgslda ftreelda

test: openmp

jni: warp_jni

warplda: $(WARPLDA_OBJ)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS) $(LIBS)

cgslda: $(CGSLDA_OBJ)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS) $(LIBS)

ftreelda: $(FTREELDA_OBJ)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS) $(LIBS)

openmp: $(OPENMP_OBJ)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS) $(LIBS)

warp_jni: $(WARPJNI_OBJ)
	$(CC) -shared -o lib$@.so $^ $(CFLAGS) $(LIBS)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS) $(LIBS)



.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(BDIR)/*
