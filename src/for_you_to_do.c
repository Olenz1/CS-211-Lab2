#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;

}

int mydgetrf(double *A, int *ipiv, int n) 
{
    return 0;
}

void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
//     int i1, j1, k1;
// 		for (i = 0; i < n; i += b)
// 			for (j = 0; j < n; j += b)
// 				for (k = 0; k < n; k += b)
// 					for (i1 = i; i1 < i + b; i1++)
// 						for (j1 = j; j1 < j + b; j1++) {
// 							register double r = C[i1 * n + j1];
// 							for (k1 = k; k1 < k + b; k1++)
// 								r += A[i1 * n + k1] * B[k1 * n + j1];
// 							C[i1 * n + j1] = r;
// 						}
    return;
}

int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
   
    return 0;
}

