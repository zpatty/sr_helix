/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mldivide.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 20:27:52
 */

/* Include Files */
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "xgetrf.h"
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const double A[100]
 *                double B[10]
 * Return Type  : void
 */
void mldivide(const double A[100], double B[10])
{
  double b_A[100];
  double temp;
  int ipiv[10];
  int i;
  int info;
  int k;
  int kAcol;
  memcpy(&b_A[0], &A[0], 100U * sizeof(double));
  xgetrf(b_A, ipiv);
  for (i = 0; i < 9; i++) {
    info = ipiv[i];
    if (info != i + 1) {
      temp = B[i];
      B[i] = B[info - 1];
      B[info - 1] = temp;
    }
  }
  for (k = 0; k < 10; k++) {
    kAcol = 10 * k;
    if (B[k] != 0.0) {
      info = k + 2;
      for (i = info; i < 11; i++) {
        B[i - 1] -= B[k] * b_A[(i + kAcol) - 1];
      }
    }
  }
  for (k = 9; k >= 0; k--) {
    kAcol = 10 * k;
    temp = B[k];
    if (temp != 0.0) {
      temp /= b_A[k + kAcol];
      B[k] = temp;
      for (i = 0; i < k; i++) {
        B[i] -= B[k] * b_A[i + kAcol];
      }
    }
  }
}

/*
 * File trailer for mldivide.c
 *
 * [EOF]
 */
