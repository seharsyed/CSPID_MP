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
 
BLLIB        = /usr/local/lib/libblas.a
CBLIB        = /usr/local/lib/libcblas.a
LAPACKLIB    = /usr/local/lib/liblapack.a
TMGLIB       = /usr/local/lib/libtmglib.a
LAPACKELIB   = /usr/local/lib/liblapacke.a


#----------------------------------------------------------------------------
#Include files (Header files)
#----------------------------------------------------------------------------

INCLUDES = ../include

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

TARGET = cspid   #Name of the files that should provide the output
RM = rm -f


#---------------------------------------------------------------------------
#CSPID compilation
#-------------------------------------------------------------------------

LIBRARIES =	 $(LAPACKELIB) $(LAPACKLIB) $(CBLIB) $(BLLIB);

COMMON_DEPS = mmio.c 

all: $(TARGET)

cspid: cspid.c $(COMMON_DEPS) $(INCLUDES)
	$(CC) $(FLAGS) -o cspid cspid.c mmio.c
	
clean:
	$(RM) $(TARGET)	

cleanobj:
		rm -f *.o



