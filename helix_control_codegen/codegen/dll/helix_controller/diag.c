/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: diag.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 24-Oct-2024 11:42:25
 */

/* Include Files */
#include "diag.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const double v[10]
 *                double d[100]
 * Return Type  : void
 */
void diag(const double v[10], double d[100])
{
  int j;
  memset(&d[0], 0, 100U * sizeof(double));
  for (j = 0; j < 10; j++) {
    d[j + 10 * j] = v[j];
  }
}

/*
 * File trailer for diag.c
 *
 * [EOF]
 */
