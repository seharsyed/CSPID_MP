# CSPID
CSPID is the library for implementing the Block Generalized Minimal Residual Method with Spectral two-grid Preconditioning techniques. 

The main folder is divided into these parts: 

1. data: contains all the data files or test matrix instances

2. include: contains all the header files

3. lib: contains the static or dynamic library (cblas, arpack, lapack)

4. src: all the source files 
 
Dependencies: ARPACK -- sudo apt-get install libsuitesparse-dev libarpack2-dev 
                        sudo apt-get install liblapack3
                        sudo apt-get install liblapack-dev
                        sudo apt-get install libopenblas-base
                        sudo apt-get install libopenblas-dev
                        sudo apt-get install liblapacke-dev libatlas-base-dev 
     
TODO:                   add dependencies of arpack, libtool, mkl
                        Write bash file for installion                   
                       
The code is compiled by accessing the src folder. 
once inside the src folder 
   
    type: 
  
  1) make clean (It removes all the previous object files)  
  
  2)  a) For Block Gmres without deflation type 
        make bgmres
      b) For Block GMRES with deflation
        make dgmres
      c) For Block GMRES with Spectral Preconditioning and Initial Deflation
         make spid
 
  Output: 
     
TODO: add .mtx file in bash installation 

 1) For BGMRES output type 
    ./bgmres filename.mtx 
 2) For DBGMRES output
    ./dbgmres filename.mtx
 3) For Spectral Preconditioning and Initial Deflation type 
    ./spid filename.mtx
 
          (Here the matrix instance used for examples is matrixA.m)


----------------------------------------------------------------------
//ARPACK -NG Instructions for linking

Libraries have been installed in:
   /usr/local/lib

If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the '-LLIBDIR'
flag during linking and do at least one of the following:
   - add LIBDIR to the 'LD_LIBRARY_PATH' environment variable
     during execution
   - add LIBDIR to the 'LD_RUN_PATH' environment variable
     during linking
   - use the '-Wl,-rpath -Wl,LIBDIR' linker flag
   - have your system administrator add LIBDIR to '/etc/ld.so.conf'

See any operating system documentation about shared libraries for
more information, such as the ld(1) and ld.so(8) manual pages.

 There are five documents within the DOCUMENT subdirectory.
  In summary,
  
    ex-nonsym.doc, ex-sym.doc  and ex-complex.doc
     -------------  ----------      --------------
    Example Templates of how to invoke the different computational
    modes offered by [D,S]NAUPD, [D,S]SAUPD and [C,Z]NAUPD.
  
    stat.doc
    --------
    File that gets timing statistics for the different parts
    of the Arnoldi update iteration codes within ARPACK.
  
    debug.doc
    ---------
    File that explains the different printing options of the
    Arnoldi update iteration codes within ARPACK.

