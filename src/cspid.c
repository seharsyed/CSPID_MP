/*************************************
  This is going to be the main file
 *************************************
 Author:  Sehar Naveed
 Institute:  Unibz
 Date: September 14, 2018
***************************************/


#include <stdio.h>

//Testing the first routine "Hello CSPID"//

int main ()
{


FILE *fpointer;

fpointer = fopen("output.txt", "w"); 
//This is writing to an output file//

fprintf(fpointer, "Hello CSPID\n");

fclose(fpointer);


/*********************************
 * TODO
*********************************/

/*I need to read data from a file
 
  The format should be 
  1) open input file
  2) get the size of the matrix
  3) Allocate the memory to the array
  4) Read data */
  
return 0;

}
