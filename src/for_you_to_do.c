    // int ib, t, j, k, j1, k1;
    // for (ib = 0; ib < n - 1; ib += b)
    // {
    //     int end = ib + b - 1;
    //     int maxind = ib;
    //     double max = fabs(A[ib * n + ib]);
    //     for (t = ib + 1; t < n; t++)
    //         if (fabs(A[t * n + ib]) > max)
    //         {
    //             maxind = t;
    //             max = fabs(A[t * n + ib]);
    //         }
    //     if (max == 0)   return -1;
    //     else {
    //         if (maxind != ib)
    //         {
    //             //save pivoting infomation
    //             int temps = ipiv[ib];
    //             ipiv[ib] = ipiv[maxind];
    //             ipiv[maxind] = temps;
    //             //swap rows
    //             int j;
    //             for (j = 0; j < n; j++)
    //             {
    //                 double tempv = A[n * ib + j];
    //                 A[ib * n + j] = A[maxind * n + j];
    //                 A[maxind * n + j] = tempv;
    //             }
    //         }
    //     }
    //     //factorization

    //     for (j = ib + 1; ib < end; ib ++)
    //     {
    //         for (k = end + 1; k < n; k ++)
    //             for (j1 = ib; ib < end; ib ++)
    //                 A[j * n + k] = A[j * n + k] / A[j * n + j1];
    //         for (k = end + 1; k < n; k ++)
    //             for (k1 = end + 1; k1 < n; k1 ++)
    //                 A[k * n + k1] -= A[k * n + j] * A[j * n + k1];
    //     }
    // }
