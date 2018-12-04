# CSPID
CSPID is the library for implementing the Block Generalized Minimal Residual Method with Spectral two-grid Preconditioning techniques. 

The main folder is divided into these parts: 

1. data: contains all the data files or test matrix instances

2. include: contains all the header files

3. lib: contains the static or dynamic library (cblas, arpack, lapack)

4. src: all the source files 

The code is compiled by accessing the src folder. 
once inside the src folder 
   
    type: 
  
  1) make clean (It removes all the previous object files)  
  
  2) make all  (compiles the C implemenation) 
  
  3) ./cspid matrixA.mtx
          (Here the matrix instance is matrixA.mtx)
  4) Enter the no of restart
            'type' 10
 5) Enter the no of restart 
          'type' 10
