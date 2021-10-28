#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, t, j, k;
    for (i = 0; i < n - 1; i ++)
    {
        //pivoting
        int maxind = i;
        int max = abs(A[i * n + i]);
        for (t = i; t < n; t ++)
            if (abs(A[t * n + i]) > max)
            {
                maxind = t;
                max = abs(A[t * n + i]);
            }
        if (max == 0)   return -1;
        else if (maxind != i)
        {
            //save pivoting infomation
            int temps = ipiv[i];
            ipiv[i] = ipiv[maxind];
            ipiv[maxind] = temps;
            //swap rows
            int j;
            for (j = 0; j < n; j ++)
            {
                double tempv = A[n * i + j];
                A[i * n + j] = A[maxind * n + j];
                A[maxind * n + j] = tempv;
            }
        }
        //factorization
        for (j = i; j < n; j ++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            for (k = i; k < n; k ++ )
                A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
        }
    }

    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    int y[n], x[n];
  //forward substitution
    int i, j;
    y[0] = B[ipiv[0]];
    for (i = 1; i < n; i ++)
    {
        int sum = 0;
        for (j = 0; j < i - 1; j ++) {
            sum += y[j] * A[i * n + j];
        }
        y[i] = B[ipiv[i]] - sum;

    }
    if (UPLO == 'L')  //backward substitution
    {
        x[n - 1] = y[n - 1] / A[n * n - n + n - 1];
        for (i = n - 1 - 1; i >= 0; i--)
        {
            int sum = 0;
            for (j = i; j < n; j++) {
                sum += x[j] * A[i * n + j];
            }
            x[i] = (y[i] - sum) / A[i * n + i];
        }
    }
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
    int i1, j1, k1;
		for (i = 0; i < n; i += b)
			for (j = 0; j < n; j += b)
				for (k = 0; k < n; k += b)
					for (i1 = i; i1 < i + b; i1++)
						for (j1 = j; j1 < j + b; j1++) {
							register double r = C[i1 * n + j1];
							for (k1 = k; k1 < k + b; k1++)
								r += A[i1 * n + k1] * B[k1 * n + j1];
							C[i1 * n + j1] = r;
						}
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, t, j, k, j1, k1;
    for (ib = 0; ib < n - 1; ib += b)
    {
        int end = ib + b - 1;
        int maxind = ib;
        int max = abs(A[ib * n + ib]);
        for (t = ib; t < n; t += b)
            if (abs(A[t * n + ib]) > max)
            {
                maxind = t;
                max = abs(A[t * n + ib]);
            }
        if (max == 0)   return -1;
        else {
            if (maxind != ib)
            {
                //save pivoting infomation
                int temps = ipiv[ib];
                ipiv[ib] = ipiv[maxind];
                ipiv[maxind] = temps;
                //swap rows
                int j;
                for (j = 0; j < n; j += b)
                {
                    double tempv = A[n * ib + j];
                    A[ib * n + j] = A[maxind * n + j];
                    A[maxind * n + j] = tempv;
                }
            }
        }
        //factorization
        // for (j = ib; j < n; j ++)
        //     A[j * n + ib] /= A[ib * n + ib];
        // for (j = ib; j < n; j ++)
        //     for (k = ib; k < n; k ++)
        //         A[j * n + k] = A[j * n + k] - A[j * n + ib] * A[ib * n + k]
        
        for (j = ib - 1; ib < end; ib += b)
        {
            for (k = end; k < n; k += b)
                for (j1 = ib - 1; j1 < end; j1 += b)
                    A[j * n + k] = A[j * n + k] / A[j * n + j1];
            for (k = end; k < n; k += b)  
                for (k1 = end; k1 < n; k1 += b)
                    A[k * n + k1] -= A[k * n + j] * A[j * n + k];
        }
    }
    return 0;
}

// for i = 1 to n-1
//      A(i+1:n,i) = A(i+1:n,i) / A(i,i)         … BLAS 1 (scale a vector)
//      A(i+1:n,i+1:n) = A(i+1:n , i+1:n )  … BLAS 2 (rank-1 update)
//               - A(i+1:n , i) * A(i , i+1:n)

//     i = i
//     j = i+1:n
//     k = i+1:n

//         for (j = i; j < n; j ++)
//             A[j * n + i] /= A[i * n + i];
//         for (j = i; j < n; j ++)
//             for (k = i; k < n; k ++)
//                 A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k]

// for ib = 1:n-1   b
//      A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n)
//      A(end+1:n , end+1:n ) = A(end+1:n , end+1:n )
//                   - A(end+1:n , ib:end) * A(ib:end , end+1:n) 

//     i = ib
//     j = ib:end
//     k = end+1:n

// for (j = ib - 1; ib < end; ib += b)
//     for (k = end; k < n; k += b)
//         A[j * n + k] = LL-1 * A[j * n + k];
//     for (k = end; k < n; k += b)  
//         for (k1 = end; k1 < n; k1 += b)
//             A[k * n + k1] -= A[k * n + j] * A[j * n + k];
