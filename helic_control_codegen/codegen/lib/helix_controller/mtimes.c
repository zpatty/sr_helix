/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mtimes.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 23-Oct-2024 10:17:18
 */

/* Include Files */
#include "mtimes.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : const double A[30]
 *                const double B[30]
 *                double C[100]
 * Return Type  : void
 */
void mtimes(const double A[30], const double B[30], double C[100])
{
  int i;
  int j;
  for (j = 0; j < 10; j++) {
    double d;
    double d1;
    double d2;
    int coffset;
    coffset = j * 10;
    d = B[j];
    d1 = B[j + 10];
    d2 = B[j + 20];
    for (i = 0; i < 10; i++) {
      int aoffset;
      aoffset = i * 3;
      C[coffset + i] =
          (A[aoffset] * d + A[aoffset + 1] * d1) + A[aoffset + 2] * d2;
    }
  }
}

/*
 * File trailer for mtimes.c
 *
 * [EOF]
 */
