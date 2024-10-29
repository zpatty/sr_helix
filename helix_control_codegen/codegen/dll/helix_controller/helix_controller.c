/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: helix_controller.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 20:27:52
 */

/* Include Files */
#include "helix_controller.h"
#include "J_r.h"
#include "MC_3_cg.h"
#include "helix_controller_rtwutil.h"
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "xgetrf.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const double q[10]
 *                const double dq[10]
 *                const double qd[10]
 *                const double dqd[10]
 *                const double ddqd[10]
 *                double d
 *                double m
 *                double r
 *                double kb
 *                double ks
 *                double bb
 *                double bs
 *                double bm
 *                double L0
 *                const double Kp[100]
 *                const double KD[100]
 *                double Kpx
 *                double KDx
 *                const double xd[3]
 *                const double dxd[3]
 *                const double dxr[3]
 *                double conv_pcc
 *                double conv_motor
 *                double tau[10]
 *                double tau_r[10]
 *                double x[3]
 *                double M[100]
 *                double C[10]
 *                double A[100]
 *                double cq[3]
 * Return Type  : void
 */
void helix_controller(const double q[10], const double dq[10],
                      const double qd[10], const double dqd[10],
                      const double ddqd[10], double d, double m, double r,
                      double kb, double ks, double bb, double bs, double bm,
                      double L0, const double Kp[100], const double KD[100],
                      double Kpx, double KDx, const double xd[3],
                      const double dxd[3], const double dxr[3], double conv_pcc,
                      double conv_motor, double tau[10], double tau_r[10],
                      double x[3], double M[100], double C[10], double A[100],
                      double cq[3])
{
  static const signed char b[9] = {-1, 0, 0, 0, -1, 0, 0, 0, -1};
  static const signed char iv[9] = {0, 1, 0, -1, 0, 0, 0, 0, 1};
  double D[100];
  double K[100];
  double KDr[100];
  double Kpr[100];
  double b_A[100];
  double conversion[100];
  double y[100];
  double J[30];
  double Jbar[30];
  double b_tmp[30];
  double b_tmp_tmp[30];
  double Ai[27];
  double b_C[10];
  double b_K[10];
  double b_qd[10];
  double c_C[10];
  double Aq[9];
  double qp[9];
  double a[3];
  double Dq_tmp;
  double a21;
  double b_d;
  double b_del_tmp;
  double d1;
  double del;
  double del_tmp;
  double theta;
  int ipiv[10];
  int b_i;
  int i;
  int ibcol;
  int kAcol;
  int r1;
  int r2;
  int r3;
  int rtemp;
  signed char b_I[100];
  MC_3_cg(q, dq, m, r, L0, M, C);
  for (i = 0; i < 3; i++) {
    b_i = 3 * i;
    qp[3 * i] = q[b_i + 1];
    Aq[b_i] = kb;
    qp[3 * i + 1] = q[b_i + 2];
    Aq[b_i + 1] = kb;
    qp[3 * i + 2] = q[b_i + 3];
    Aq[b_i + 2] = ks;
  }
  b_C[0] = 0.0;
  memcpy(&b_C[1], &Aq[0], 9U * sizeof(double));
  memset(&K[0], 0, 100U * sizeof(double));
  for (r3 = 0; r3 < 10; r3++) {
    K[r3 + 10 * r3] = b_C[r3];
  }
  for (rtemp = 0; rtemp < 3; rtemp++) {
    ibcol = rtemp * 3;
    Aq[ibcol] = bb;
    Aq[ibcol + 1] = bb;
    Aq[ibcol + 2] = bs;
  }
  b_C[0] = bm;
  memcpy(&b_C[1], &Aq[0], 9U * sizeof(double));
  memset(&D[0], 0, 100U * sizeof(double));
  for (r3 = 0; r3 < 10; r3++) {
    D[r3 + 10 * r3] = b_C[r3];
  }
  /*  Fci = cell(1,4); */
  /*  Ai = cell(1,4); */
  for (i = 0; i < 3; i++) {
    double b_Aq[9];
    double Dq;
    double L_tmp;
    b_d = qp[3 * i];
    d1 = qp[3 * i + 1];
    del_tmp = b_d * b_d;
    b_del_tmp = d1 * d1;
    del = sqrt(del_tmp + b_del_tmp);
    theta = del / d;
    Dq_tmp = sin(del);
    Dq = del - Dq_tmp;
    L_tmp = qp[3 * i + 2] + L0;
    if ((b_d < 1.0E-5) && (d1 < 1.0E-5)) {
      cq[i] = L_tmp / 3.0;
      for (b_i = 0; b_i < 9; b_i++) {
        Aq[b_i] = iv[b_i];
      }
    } else {
      cq[i] = 2.0 * (L_tmp / theta - d) * sin(theta / 6.0);
      a21 = rt_powd_snf(del, 3.0);
      Aq[0] = b_d * d1 * Dq / a21;
      Aq[3] = (-del_tmp * del - b_del_tmp * Dq_tmp) / a21;
      Aq[6] = b_d * Dq * L_tmp / a21;
      Aq[1] = (b_del_tmp * del + del_tmp * Dq_tmp) / a21;
      Aq[4] = -b_d * d1 * Dq / a21;
      Aq[7] = d1 * Dq * L_tmp / a21;
      Aq[2] = 0.0;
      Aq[5] = 0.0;
      Aq[8] = Dq_tmp / del;
    }
    for (b_i = 0; b_i < 3; b_i++) {
      b_d = Aq[b_i];
      d1 = Aq[b_i + 3];
      theta = Aq[b_i + 6];
      for (r2 = 0; r2 < 3; r2++) {
        b_Aq[b_i + 3 * r2] =
            (b_d * (double)b[3 * r2] + d1 * (double)b[3 * r2 + 1]) +
            theta * (double)b[3 * r2 + 2];
      }
    }
    Aq[0] = d * 0.49999999999999994;
    Aq[3] = d * 0.49999999999999994;
    Aq[6] = -d;
    Aq[1] = -d * 0.86602540378443871;
    Aq[4] = d * 0.86602540378443871;
    Aq[7] = 0.0;
    Aq[2] = 1.0;
    Aq[5] = 1.0;
    Aq[8] = 1.0;
    for (b_i = 0; b_i < 3; b_i++) {
      b_d = b_Aq[b_i];
      d1 = b_Aq[b_i + 3];
      theta = b_Aq[b_i + 6];
      for (r2 = 0; r2 < 3; r2++) {
        Ai[(b_i + 3 * r2) + 9 * i] =
            (b_d * Aq[3 * r2] + d1 * Aq[3 * r2 + 1]) + theta * Aq[3 * r2 + 2];
      }
    }
  }
  /*  A = blkdiag(1,Ai(:,:,1),Ai(:,:,2),Ai(:,:,3)); */
  memset(&A[0], 0, 100U * sizeof(double));
  A[0] = 1.0;
  a[0] = conv_pcc;
  a[1] = conv_pcc;
  a[2] = conv_pcc;
  for (i = 0; i < 3; i++) {
    b_i = 3 * i;
    for (ibcol = 0; ibcol < 3; ibcol++) {
      r1 = 3 * ibcol + 9 * i;
      rtemp = ibcol + b_i;
      kAcol = b_i + 10 * (rtemp + 1);
      A[kAcol + 1] = Ai[r1];
      A[kAcol + 2] = Ai[r1 + 1];
      A[kAcol + 3] = Ai[r1 + 2];
      Aq[rtemp] = a[ibcol];
    }
  }
  b_C[0] = conv_motor;
  memcpy(&b_C[1], &Aq[0], 9U * sizeof(double));
  memset(&conversion[0], 0, 100U * sizeof(double));
  for (r3 = 0; r3 < 10; r3++) {
    conversion[r3 + 10 * r3] = b_C[r3];
    b_d = 0.0;
    b_qd[r3] = qd[r3] - q[r3];
    d1 = 0.0;
    theta = 0.0;
    for (b_i = 0; b_i < 10; b_i++) {
      r2 = r3 + 10 * b_i;
      b_d += M[r2] * ddqd[b_i];
      theta += K[r2] * qd[b_i];
      d1 += D[r2] * dqd[b_i];
    }
    b_d += C[r3];
    b_C[r3] = b_d;
    c_C[r3] = (b_d + theta) + d1;
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = 0.0;
    for (r2 = 0; r2 < 10; r2++) {
      b_d += Kp[b_i + 10 * r2] * b_qd[r2];
    }
    b_K[b_i] = b_d;
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_qd[b_i] = dqd[b_i] - dq[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = 0.0;
    for (r2 = 0; r2 < 10; r2++) {
      b_d += KD[b_i + 10 * r2] * b_qd[r2];
    }
    b_C[b_i] = (c_C[b_i] + b_K[b_i]) + b_d;
  }
  mldivide(A, b_C);
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = 0.0;
    for (r2 = 0; r2 < 10; r2++) {
      b_d += conversion[b_i + 10 * r2] * b_C[r2];
    }
    tau[b_i] = b_d;
  }
  memcpy(&Kpr[0], &Kp[0], 100U * sizeof(double));
  Kpr[0] = 0.0;
  memcpy(&KDr[0], &KD[0], 100U * sizeof(double));
  KDr[0] = 0.0;
  memcpy(&c_C[0], &q[0], 10U * sizeof(double));
  J_r(c_C, L0, d, J, x);
  for (b_i = 0; b_i < 3; b_i++) {
    for (r2 = 0; r2 < 10; r2++) {
      b_tmp_tmp[r2 + 10 * b_i] = J[b_i + 3 * r2];
    }
  }
  memcpy(&b_tmp[0], &b_tmp_tmp[0], 30U * sizeof(double));
  memcpy(&b_A[0], &M[0], 100U * sizeof(double));
  xgetrf(b_A, ipiv);
  for (i = 0; i < 9; i++) {
    b_i = ipiv[i];
    if (b_i != i + 1) {
      theta = b_tmp[i];
      b_tmp[i] = b_tmp[b_i - 1];
      b_tmp[b_i - 1] = theta;
      theta = b_tmp[i + 10];
      b_tmp[i + 10] = b_tmp[b_i + 9];
      b_tmp[b_i + 9] = theta;
      theta = b_tmp[i + 20];
      b_tmp[i + 20] = b_tmp[b_i + 19];
      b_tmp[b_i + 19] = theta;
    }
  }
  for (r3 = 0; r3 < 3; r3++) {
    ibcol = 10 * r3;
    for (rtemp = 0; rtemp < 10; rtemp++) {
      kAcol = 10 * rtemp;
      b_i = rtemp + ibcol;
      if (b_tmp[b_i] != 0.0) {
        r2 = rtemp + 2;
        for (i = r2; i < 11; i++) {
          r1 = (i + ibcol) - 1;
          b_tmp[r1] -= b_tmp[b_i] * b_A[(i + kAcol) - 1];
        }
      }
    }
  }
  for (r3 = 0; r3 < 3; r3++) {
    ibcol = 10 * r3;
    for (rtemp = 9; rtemp >= 0; rtemp--) {
      kAcol = 10 * rtemp;
      b_i = rtemp + ibcol;
      b_d = b_tmp[b_i];
      if (b_d != 0.0) {
        b_tmp[b_i] = b_d / b_A[rtemp + kAcol];
        for (i = 0; i < rtemp; i++) {
          r1 = i + ibcol;
          b_tmp[r1] -= b_tmp[b_i] * b_A[i + kAcol];
        }
      }
    }
  }
  for (b_i = 0; b_i < 3; b_i++) {
    for (r2 = 0; r2 < 3; r2++) {
      b_d = 0.0;
      for (rtemp = 0; rtemp < 10; rtemp++) {
        b_d += J[b_i + 3 * rtemp] * b_tmp[rtemp + 10 * r2];
      }
      Aq[b_i + 3 * r2] = b_d;
    }
  }
  r1 = 0;
  r2 = 1;
  r3 = 2;
  theta = fabs(Aq[0]);
  a21 = fabs(Aq[1]);
  if (a21 > theta) {
    theta = a21;
    r1 = 1;
    r2 = 0;
  }
  if (fabs(Aq[2]) > theta) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }
  Aq[r2] /= Aq[r1];
  Aq[r3] /= Aq[r1];
  Aq[r2 + 3] -= Aq[r2] * Aq[r1 + 3];
  Aq[r3 + 3] -= Aq[r3] * Aq[r1 + 3];
  Aq[r2 + 6] -= Aq[r2] * Aq[r1 + 6];
  Aq[r3 + 6] -= Aq[r3] * Aq[r1 + 6];
  if (fabs(Aq[r3 + 3]) > fabs(Aq[r2 + 3])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }
  Aq[r3 + 3] /= Aq[r2 + 3];
  Aq[r3 + 6] -= Aq[r3 + 3] * Aq[r2 + 6];
  qp[3 * r1] = 1.0 / Aq[r1];
  theta = Aq[r1 + 3];
  qp[3 * r2] = 0.0 - qp[3 * r1] * theta;
  a21 = Aq[r1 + 6];
  qp[3 * r3] = 0.0 - qp[3 * r1] * a21;
  del_tmp = Aq[r2 + 3];
  qp[3 * r2] /= del_tmp;
  b_del_tmp = Aq[r2 + 6];
  qp[3 * r3] -= qp[3 * r2] * b_del_tmp;
  del = Aq[r3 + 6];
  qp[3 * r3] /= del;
  Dq_tmp = Aq[r3 + 3];
  qp[3 * r2] -= qp[3 * r3] * Dq_tmp;
  qp[3 * r1] -= qp[3 * r3] * Aq[r3];
  qp[3 * r1] -= qp[3 * r2] * Aq[r2];
  ibcol = 3 * r1 + 1;
  qp[ibcol] = 0.0 / Aq[r1];
  rtemp = 3 * r2 + 1;
  qp[rtemp] = 1.0 - qp[ibcol] * theta;
  kAcol = 3 * r3 + 1;
  qp[kAcol] = 0.0 - qp[ibcol] * a21;
  qp[rtemp] /= del_tmp;
  qp[kAcol] -= qp[rtemp] * b_del_tmp;
  qp[kAcol] /= del;
  qp[rtemp] -= qp[kAcol] * Dq_tmp;
  qp[ibcol] -= qp[kAcol] * Aq[r3];
  qp[ibcol] -= qp[rtemp] * Aq[r2];
  ibcol = 3 * r1 + 2;
  qp[ibcol] = 0.0 / Aq[r1];
  rtemp = 3 * r2 + 2;
  qp[rtemp] = 0.0 - qp[ibcol] * theta;
  kAcol = 3 * r3 + 2;
  qp[kAcol] = 1.0 - qp[ibcol] * a21;
  qp[rtemp] /= del_tmp;
  qp[kAcol] -= qp[rtemp] * b_del_tmp;
  qp[kAcol] /= del;
  qp[rtemp] -= qp[kAcol] * Dq_tmp;
  qp[ibcol] -= qp[kAcol] * Aq[r3];
  qp[ibcol] -= qp[rtemp] * Aq[r2];
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = b_tmp[b_i];
    d1 = b_tmp[b_i + 10];
    theta = b_tmp[b_i + 20];
    for (r2 = 0; r2 < 3; r2++) {
      Jbar[b_i + 10 * r2] =
          (b_d * qp[3 * r2] + d1 * qp[3 * r2 + 1]) + theta * qp[3 * r2 + 2];
    }
  }
  for (r3 = 0; r3 < 10; r3++) {
    rtemp = r3 * 10;
    b_d = Jbar[r3];
    d1 = Jbar[r3 + 10];
    theta = Jbar[r3 + 20];
    for (i = 0; i < 10; i++) {
      ibcol = i * 3;
      a21 = (J[ibcol] * b_d + J[ibcol + 1] * d1) + J[ibcol + 2] * theta;
      r1 = rtemp + i;
      b_A[r1] = a21;
      y[r1] = a21;
    }
  }
  memset(&b_I[0], 0, 100U * sizeof(signed char));
  for (rtemp = 0; rtemp < 10; rtemp++) {
    b_I[rtemp + 10 * rtemp] = 1;
    b_d = 0.0;
    d1 = 0.0;
    for (b_i = 0; b_i < 10; b_i++) {
      r2 = rtemp + 10 * b_i;
      b_d += K[r2] * q[b_i];
      d1 += D[r2] * dq[b_i];
    }
    b_K[rtemp] = b_d + d1;
    b_d = b_tmp_tmp[rtemp];
    d1 = b_tmp_tmp[rtemp + 10];
    theta = b_tmp_tmp[rtemp + 20];
    for (b_i = 0; b_i < 3; b_i++) {
      b_tmp[rtemp + 10 * b_i] =
          (b_d * qp[3 * b_i] + d1 * qp[3 * b_i + 1]) + theta * qp[3 * b_i + 2];
    }
  }
  a[0] = Kpx * (xd[0] - x[0]) + KDx * (dxd[0] - dxr[0]);
  a[1] = Kpx * (xd[1] - x[1]) + KDx * (dxd[1] - dxr[1]);
  a[2] = Kpx * (xd[2] - x[2]) + KDx * (dxd[2] - dxr[2]);
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = 0.0;
    for (r2 = 0; r2 < 10; r2++) {
      b_d += b_A[b_i + 10 * r2] * b_K[r2];
    }
    b_C[b_i] = (C[b_i] + b_d) + ((b_tmp[b_i] * a[0] + b_tmp[b_i + 10] * a[1]) +
                                 b_tmp[b_i + 20] * a[2]);
  }
  mldivide(A, b_C);
  for (b_i = 0; b_i < 100; b_i++) {
    Kpr[b_i] = -Kpr[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = 0.0;
    d1 = 0.0;
    for (r2 = 0; r2 < 10; r2++) {
      rtemp = b_i + 10 * r2;
      b_d += Kpr[rtemp] * q[r2];
      d1 += KDr[rtemp] * dq[r2];
    }
    b_qd[b_i] = d1;
    b_K[b_i] = b_d;
  }
  for (b_i = 0; b_i < 100; b_i++) {
    y[b_i] = (double)b_I[b_i] - y[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_K[b_i] -= b_qd[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = 0.0;
    for (r2 = 0; r2 < 10; r2++) {
      b_d += y[b_i + 10 * r2] * b_K[r2];
    }
    b_qd[b_i] = b_C[b_i] + b_d;
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_d = 0.0;
    for (r2 = 0; r2 < 10; r2++) {
      b_d += conversion[b_i + 10 * r2] * b_qd[r2];
    }
    tau_r[b_i] = b_d;
  }
}

/*
 * File trailer for helix_controller.c
 *
 * [EOF]
 */
