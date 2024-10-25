/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: blkdiag.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 23-Oct-2024 10:17:18
 */

/* Include Files */
#include "blkdiag.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const double varargin_2[9]
 *                const double varargin_3[9]
 *                const double varargin_4[9]
 *                double y[100]
 * Return Type  : void
 */
void blkdiag(const double varargin_2[9], const double varargin_3[9],
             const double varargin_4[9], double y[100])
{
  int i;
  memset(&y[0], 0, 100U * sizeof(double));
  y[0] = 1.0;
  for (i = 0; i < 3; i++) {
    int b_y_tmp;
    int c_y_tmp;
    int d_y_tmp;
    int y_tmp;
    y_tmp = 10 * (i + 1);
    y[y_tmp + 1] = varargin_2[3 * i];
    b_y_tmp = 10 * (i + 4);
    y[b_y_tmp + 4] = varargin_3[3 * i];
    c_y_tmp = 10 * (i + 7);
    y[c_y_tmp + 7] = varargin_4[3 * i];
    d_y_tmp = 3 * i + 1;
    y[y_tmp + 2] = varargin_2[d_y_tmp];
    y[b_y_tmp + 5] = varargin_3[d_y_tmp];
    y[c_y_tmp + 8] = varargin_4[d_y_tmp];
    d_y_tmp = 3 * i + 2;
    y[y_tmp + 3] = varargin_2[d_y_tmp];
    y[b_y_tmp + 6] = varargin_3[d_y_tmp];
    y[c_y_tmp + 9] = varargin_4[d_y_tmp];
  }
}

/*
 * File trailer for blkdiag.c
 *
 * [EOF]
 */
