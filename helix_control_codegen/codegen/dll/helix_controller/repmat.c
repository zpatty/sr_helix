/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: repmat.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 24-Oct-2024 11:42:25
 */

/* Include Files */
#include "repmat.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : const double a[3]
 *                double b[9]
 * Return Type  : void
 */
void repmat(const double a[3], double b[9])
{
  double d;
  double d1;
  double d2;
  int itilerow;
  d = a[0];
  d1 = a[1];
  d2 = a[2];
  for (itilerow = 0; itilerow < 3; itilerow++) {
    int ibcol;
    ibcol = itilerow * 3;
    b[ibcol] = d;
    b[ibcol + 1] = d1;
    b[ibcol + 2] = d2;
  }
}

/*
 * File trailer for repmat.c
 *
 * [EOF]
 */
