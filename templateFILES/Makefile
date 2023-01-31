#Makefile
UNAME:=$(shell uname)

SRC=src
OBJ=obj
LIBS=libs
BIN=bin

TOOLCHAIN=gnu
MPICXX=
CCXXFLAGS=-std=c++17
CCDEFINEFLAGS=-DHDF5_FLOAT	-DHDF5_SINGLE_FILE	-DCHECKS_ON
HDF5_INSTALL=
ARMADILLO_INSTALL=

INCLUDE_HDF5=-I$(HDF5_INSTALL)/include
INCLUDES=-I$(ARMADILLO_INSTALL)/include	-I$(ARMADILLO_INSTALL)/lib	-I$(ARMADILLO_INSTALL)/lib64
EXTLIB=-L$(HDF5_INSTALL)/lib
LIBSHDF=$(EXTLIB) $(HDF5_INSTALL)/lib/libhdf5.a
LIB_FLAGS=-L$(ARMADILLO_INSTALL)/lib	-L$(ARMADILLO_INSTALL)/lib64
LIB=-lz	-ldl	-lm	-lhdf5_cpp	-larmadillo #-lgfortran

# Set the relevant flag based on the compiler TOOLCHAIN=gnu or ibm
ifeq ($(TOOLCHAIN),gnu)
	OMP_LIB_FLAGS=-fopenmp
else
	OMP_LIB_FLAGS=-qopenmp
endif
CCFLAGS=$(CCXXFLAGS) $(OMP_LIB_FLAGS)	-O3

BINARIES = $(addprefix $(OBJ)/,\
main.o	            \
types.o   			    \
initialize.o		    \
units.o						  \
initDistribution.o  \
outputHDF5.o 			  \
PIC.o 						  \
particleBC.o			  \
collisionOperator.o \
fieldSolve.o				\
rfOperator.o				\
mpi_main.o)

all: $(BIN)/PICOS++

$(OBJ)/%.o : $(SRC)/%.cpp
	$(MPICXX)	-c	$(CCFLAGS)	$(CCDEFINEFLAGS)	$<	$(INCLUDES)	$(INCLUDE_HDF5)	$(LIB_FLAGS) -o $@

$(BIN)/PICOS++: $(BINARIES)
	$(MPICXX)	$(BINARIES)	-o	$@	$(LIB_FLAGS)	$(LIBSHDF)	$(LIB)	$(OMP_LIB_FLAGS)

info:
	@echo	"COMPILING IN: "$(UNAME)

clean:
	@echo	"CLEANING PREVIOUS COMPILATION"
	rm -f $(BIN)/PICOS++	$(OBJ)/*.o	$(SRC)/*~
