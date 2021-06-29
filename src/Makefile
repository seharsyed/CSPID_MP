#
# Makefile.in
#
 
#-----------------------------------------------------------------------------
# Shell
#-----------------------------------------------------------------------------

SHELL = /bin/sh

#-----------------------------------------------------------------------------
# Platform
#-----------------------------------------------------------------------------

PLAT = LINUX

#-----------------------------------------------------------------------------
# Compilers:
#-----------------------------------------------------------------------------

CC = gcc
FC = gfortran
LOADER = $(FC)


#-----------------------------------------------------------------------------
# Libraries
#-----------------------------------------------------------------------------

#I have added the path to the static libaries, will resolve the dynamic library 
#linking issue later. 
BLAS = -lblas
CBLAS = -lcblas
LIBES = -lm
LAPACK = -llapack
LAPACKE = -llapacke
#SPBLAS = /home/academic/snaveed/CSPID/spblas_0_8/libsparseblas.a 
#----------------------------------------------------------------------------
#Include files (Header files)
#----------------------------------------------------------------------------

INCLUDES = -I../include

#-----------------------------------------------------------------------------
# Flags for Compilers
#-----------------------------------------------------------------------------

CFLAGS = -Wall
FFLAGS = -O3  

#-----------------------------------------------------------------------------
# Archive programs and flags
#-----------------------------------------------------------------------------

ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlibi

dlvl = ./.

#--------------------------------------------------------------------------
#Object Files
#--------------------------------------------------------------------------

TARGET = spid   #Name of the files that should provide the output
RM = rm -f


#---------------------------------------------------------------------------
#CSPID compilation
#-------------------------------------------------------------------------

#LIB = $(LIBES) $(LAPACKE) $(LAPACK) $(BLAS)
LIB =  $(LIBES) $(LAPACKE) $(CBLAS) $(SPBLAS) $(LAPACK) $(BLAS)
COMMON_DEPS = mmio.c 

all: $(TARGET)

dbgmres: dbgmres.c $(COMMON_DEPS) 
	$(CC) $(FLAGS) $(INCLUDES) -o dbgmres dbgmres.c mmio.c $(LIB)

bgmres: bgmres.c $(COMMON_DEPS) 
	$(CC) $(FLAGS) $(INCLUDES) -o bgmres bgmres.c mmio.c $(LIB)

spid: spid.c $(COMMON_DEPS)
	 $(CC) $(FLAGS) $(INCLUDES) -o spid spid.c mmio.c $(LIB)
	
clean:
	$(RM) $(TARGET)	

cleanobj:
		rm -f *.o



