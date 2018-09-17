/*************************************
  This is going to be the main file
 *************************************
 Author:  Sehar Naveed
 Institute:  Unibz
 Date: September 14, 2018
***************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
//Testing the first routine "Hello CSPID"//

int main (int arg, char** args)
{

    //printf("%s :\n\n\n",args[1]);
//FILE *fpointer;

//fpointer = fopen("output.txt", "w");
//This is writing to an output file//

//fprintf(fpointer, "Hello CSPID\n");

//fclose(fpointer);


/*********************************
 * TODO
*********************************/

/*I need to read data from a file
 
  The format should be 
  1) open input file
  2) get the size of the matrix
  3) Allocate the memory to the array
  4) Read data */
    
    
    /**
     File Reading
     **/
    
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    if(arg <2){
        printf("Missing File name\n");
        return 0;
    }
    fp = fopen(args[1], "r");
    if (fp == NULL)
        exit(0);
    
//    int r = 3, c = 3, i=0, j, count;
//
//    float *arr[r];
//    for (i=0; i<r; i++)
//        arr[i] = (float *)malloc(c * sizeof(float));
//
//
//    for(i = 0; i < r; i++)
//        arr[i] = (*arr + c * i);
    //char* split = NULL;
    int i=0, j=0;
    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s\n", line);
        char* split = strsep(&line, " ") ;
        printf("%s\t", split);
        split = strsep(&line, " ") ;
        printf("%s\t", split);
        split = strsep(&line, " ") ;
        printf("%s\t\n", split);
        
    }
    
//
//    for (i = 0; i <  r; i++){
//        for (j = 0; j < c; j++){
//            printf("%f \t ",arr[i][j]);
//        }
//        printf("\n");
//    }
//    free(arr);
    fclose(fp);
    if (line)
        free(line);
    exit(0);
  
return 0;

}
