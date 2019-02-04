# CSPID
CSPID is the library for implementing the Block Generalized Minimal Residual Method with Spectral two-grid Preconditioning techniques. 

The main folder is divided into these parts: 

1. data: contains all the data files or test matrix instances

2. include: contains all the header files

3. lib: contains the static or dynamic library (cblas, arpack, lapack)

4. src: all the source files 

Dependencies: ARPACK -- apt-get install libsuitesparse-dev libarpack2-dev 
             #WIll add dependencies for lapacke and cblas

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
     

 1) For BGMRES output type 
    ./bgmres filename.mtx 
 2) For DBGMRES output
    ./dbgmres filename.mtx
 3) For Spectral Preconditioning and Initial Deflation type 
    ./spid filename.mtx
 
          (Here the matrix instance used for examples is matrixA.mtx)
