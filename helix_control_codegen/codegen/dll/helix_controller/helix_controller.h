/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: helix_controller.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 20:27:52
 */

#ifndef HELIX_CONTROLLER_H
#define HELIX_CONTROLLER_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void helix_controller(
    const double q[10], const double dq[10], const double qd[10],
    const double dqd[10], const double ddqd[10], double d, double m, double r,
    double kb, double ks, double bb, double bs, double bm, double L0,
    const double Kp[100], const double KD[100], double Kpx, double KDx,
    const double xd[3], const double dxd[3], const double dxr[3],
    double conv_pcc, double conv_motor, double tau[10], double tau_r[10],
    double x[3], double M[100], double C[10], double A[100], double cq[3]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for helix_controller.h
 *
 * [EOF]
 */
