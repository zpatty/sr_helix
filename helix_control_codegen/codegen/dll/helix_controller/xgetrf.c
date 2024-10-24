/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgetrf.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 24-Oct-2024 11:42:25
 */

/* Include Files */
#include "xgetrf.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : double A[100]
 *                int ipiv[10]
 * Return Type  : int
 */
int xgetrf(double A[100], int ipiv[10])
{
  int i;
  int info;
  int j;
  int jA;
  int jp1j;
  int k;
  for (i = 0; i < 10; i++) {
    ipiv[i] = i + 1;
  }
  info = 0;
  for (j = 0; j < 9; j++) {
    double smax;
    int a;
    int b_tmp;
    int mmj_tmp;
    mmj_tmp = 8 - j;
    b_tmp = j * 11;
    jp1j = b_tmp + 2;
    jA = 10 - j;
    a = 0;
    smax = fabs(A[b_tmp]);
    for (k = 2; k <= jA; k++) {
      double s;
      s = fabs(A[(b_tmp + k) - 1]);
      if (s > smax) {
        a = k - 1;
        smax = s;
      }
    }
    if (A[b_tmp + a] != 0.0) {
      if (a != 0) {
        jA = j + a;
        ipiv[j] = jA + 1;
        for (k = 0; k < 10; k++) {
          a = j + k * 10;
          smax = A[a];
          i = jA + k * 10;
          A[a] = A[i];
          A[i] = smax;
        }
      }
      i = (b_tmp - j) + 10;
      for (jA = jp1j; jA <= i; jA++) {
        A[jA - 1] /= A[b_tmp];
      }
    } else {
      info = j + 1;
    }
    jA = b_tmp;
    for (jp1j = 0; jp1j <= mmj_tmp; jp1j++) {
      smax = A[(b_tmp + jp1j * 10) + 10];
      if (smax != 0.0) {
        i = jA + 12;
        a = (jA - j) + 20;
        for (k = i; k <= a; k++) {
          A[k - 1] += A[((b_tmp + k) - jA) - 11] * -smax;
        }
      }
      jA += 10;
    }
  }
  if ((info == 0) && (!(A[99] != 0.0))) {
    info = 10;
  }
  return info;
}

/*
 * File trailer for xgetrf.c
 *
 * [EOF]
 */
