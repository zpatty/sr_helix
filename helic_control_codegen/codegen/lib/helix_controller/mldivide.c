/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mldivide.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 23-Oct-2024 10:17:18
 */

/* Include Files */
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "xgetrf.h"
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const double A[100]
 *                double B[30]
 * Return Type  : void
 */
void b_mldivide(const double A[100], double B[30])
{
  double b_A[100];
  double temp;
  int ipiv[10];
  int b_i;
  int i;
  int info;
  int j;
  int jBcol;
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
      temp = B[i + 10];
      B[i + 10] = B[info + 9];
      B[info + 9] = temp;
      temp = B[i + 20];
      B[i + 20] = B[info + 19];
      B[info + 19] = temp;
    }
  }
  for (j = 0; j < 3; j++) {
    jBcol = 10 * j;
    for (k = 0; k < 10; k++) {
      kAcol = 10 * k;
      info = k + jBcol;
      if (B[info] != 0.0) {
        b_i = k + 2;
        for (i = b_i; i < 11; i++) {
          int i1;
          i1 = (i + jBcol) - 1;
          B[i1] -= B[info] * b_A[(i + kAcol) - 1];
        }
      }
    }
  }
  for (j = 0; j < 3; j++) {
    jBcol = 10 * j;
    for (k = 9; k >= 0; k--) {
      kAcol = 10 * k;
      info = k + jBcol;
      temp = B[info];
      if (temp != 0.0) {
        B[info] = temp / b_A[k + kAcol];
        for (i = 0; i < k; i++) {
          b_i = i + jBcol;
          B[b_i] -= B[info] * b_A[i + kAcol];
        }
      }
    }
  }
}

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
