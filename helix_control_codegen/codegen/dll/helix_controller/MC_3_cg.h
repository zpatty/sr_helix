/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: MC_3_cg.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 23:34:49
 */

#ifndef MC_3_CG_H
#define MC_3_CG_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void MC_3_cg(const double q[10], const double qd[10], double m, double r,
             double L0, double M[100], double C[10]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for MC_3_cg.h
 *
 * [EOF]
 */
