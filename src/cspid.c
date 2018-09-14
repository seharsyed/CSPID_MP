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

fprintf(fpointer, "Hello CSPID\n");

fclose(fpointer);

return 0;

}
