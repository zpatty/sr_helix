/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: helix_controller.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 17:40:15
 */

/* Include Files */
#include "helix_controller.h"
#include "MC_3_cg.h"
#include "cosd.h"
#include "helix_controller_emxutil.h"
#include "helix_controller_rtwutil.h"
#include "helix_controller_types.h"
#include "rt_nonfinite.h"
#include "xgeqp3.h"
#include <math.h>

/* Function Declarations */
static void binary_expand_op(double in1[10], const emxArray_real_T *in2,
                             const emxArray_real_T *in3,
                             const emxArray_real_T *in4,
                             const emxArray_real_T *in5, const double in6[100],
                             const double in7[10], const double in8[10],
                             const double in9[100], const double in10[10],
                             const double in11[10]);

/* Function Definitions */
/*
 * Arguments    : double in1[10]
 *                const emxArray_real_T *in2
 *                const emxArray_real_T *in3
 *                const emxArray_real_T *in4
 *                const emxArray_real_T *in5
 *                const double in6[100]
 *                const double in7[10]
 *                const double in8[10]
 *                const double in9[100]
 *                const double in10[10]
 *                const double in11[10]
 * Return Type  : void
 */
static void binary_expand_op(double in1[10], const emxArray_real_T *in2,
                             const emxArray_real_T *in3,
                             const emxArray_real_T *in4,
                             const emxArray_real_T *in5, const double in6[100],
                             const double in7[10], const double in8[10],
                             const double in9[100], const double in10[10],
                             const double in11[10])
{
  double b_in7[10];
  const double *in2_data;
  const double *in3_data;
  const double *in4_data;
  const double *in5_data;
  double d;
  int i;
  int i1;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  int stride_3_0;
  in5_data = in5->data;
  in4_data = in4->data;
  in3_data = in3->data;
  in2_data = in2->data;
  for (i = 0; i < 10; i++) {
    b_in7[i] = in7[i] - in8[i];
  }
  for (i = 0; i < 10; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      d += in6[i + 10 * i1] * b_in7[i1];
    }
    in1[i] = d;
  }
  for (i = 0; i < 10; i++) {
    b_in7[i] = in10[i] - in11[i];
  }
  stride_0_0 = (in2->size[0] != 1);
  stride_1_0 = (in3->size[0] != 1);
  stride_2_0 = (in4->size[0] != 1);
  stride_3_0 = (in5->size[0] != 1);
  for (i = 0; i < 10; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      d += in9[i + 10 * i1] * b_in7[i1];
    }
    in1[i] = ((((in2_data[i * stride_0_0] + in3_data[i * stride_1_0]) +
                in4_data[i * stride_2_0]) +
               in5_data[i * stride_3_0]) +
              in1[i]) +
             d;
  }
}

/*
 * Arguments    : const double q[10]
 *                const double dq[10]
 *                const double qd[10]
 *                const double dqd[10]
 *                const double ddqd[10]
 *                double N
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
 *                emxArray_real_T *tau
 *                double *tau_r
 *                double x[36]
 *                emxArray_real_T *M
 *                emxArray_real_T *C
 *                emxArray_real_T *A
 *                emxArray_real_T *cq
 * Return Type  : void
 */
void helix_controller(const double q[10], const double dq[10],
                      const double qd[10], const double dqd[10],
                      const double ddqd[10], double N, double d, double m,
                      double r, double kb, double ks, double bb, double bs,
                      double bm, double L0, const double Kp[100],
                      const double KD[100], double Kpx, double KDx,
                      const double xd[3], const double dxd[3],
                      const double dxr[3], double conv_pcc, double conv_motor,
                      emxArray_real_T *tau, double *tau_r, double x[36],
                      emxArray_real_T *M, emxArray_real_T *C,
                      emxArray_real_T *A, emxArray_real_T *cq)
{
  static const signed char c_b[9] = {-1, 0, 0, 0, -1, 0, 0, 0, -1};
  static const signed char iv[9] = {0, 1, 0, -1, 0, 0, 0, 0, 1};
  emxArray_int32_T *jpvt;
  emxArray_real_T *Ai;
  emxArray_real_T *D;
  emxArray_real_T *J;
  emxArray_real_T *K;
  emxArray_real_T *b;
  emxArray_real_T *b_C;
  emxArray_real_T *b_b;
  emxArray_real_T *b_v;
  emxArray_real_T *conversion;
  emxArray_real_T *qp;
  emxArray_real_T *v;
  double c_C[10];
  double b_del_tmp;
  double dx;
  double theta;
  double *Ai_data;
  double *C_data;
  double *D_data;
  double *K_data;
  double *M_data;
  double *b_data;
  double *conversion_data;
  double *cq_data;
  double *qp_data;
  int b_i;
  int i;
  int ibcol;
  int ibtile;
  int itilerow;
  int j;
  int jp1j;
  int k;
  int loop_ub_tmp;
  int mc_tmp;
  int *jpvt_data;
  (void)Kpx;
  (void)KDx;
  (void)xd;
  (void)dxd;
  (void)dxr;
  emxInit_real_T(&qp, 2);
  i = qp->size[0] * qp->size[1];
  qp->size[0] = 3;
  jp1j = (int)(N - 1.0);
  qp->size[1] = (int)(N - 1.0);
  emxEnsureCapacity_real_T(qp, i);
  qp_data = qp->data;
  ibtile = 3 * (int)(N - 1.0);
  for (i = 0; i < ibtile; i++) {
    qp_data[i] = 0.0;
  }
  for (b_i = 0; b_i < jp1j; b_i++) {
    theta = 3.0 * (((double)b_i + 1.0) - 1.0) + 2.0;
    qp_data[3 * b_i] = q[(int)theta - 1];
    qp_data[3 * b_i + 1] = q[(int)(theta + 1.0) - 1];
    qp_data[3 * b_i + 2] = q[(int)(theta + 2.0) - 1];
  }
  emxInit_real_T(&J, 2);
  MC_3_cg(q, dq, m, r, L0, N, M, C, J, x);
  C_data = C->data;
  M_data = M->data;
  emxFree_real_T(&J);
  emxInit_real_T(&b, 1);
  i = b->size[0];
  b->size[0] = ibtile;
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  for (itilerow = 0; itilerow < jp1j; itilerow++) {
    ibcol = itilerow * 3;
    b_data[ibcol] = kb;
    b_data[ibcol + 1] = kb;
    b_data[ibcol + 2] = ks;
  }
  emxInit_real_T(&v, 1);
  i = v->size[0];
  v->size[0] = b->size[0] + 1;
  emxEnsureCapacity_real_T(v, i);
  Ai_data = v->data;
  Ai_data[0] = 0.0;
  ibcol = b->size[0];
  for (i = 0; i < ibcol; i++) {
    Ai_data[i + 1] = b_data[i];
  }
  ibcol = v->size[0];
  emxInit_real_T(&K, 2);
  i = K->size[0] * K->size[1];
  K->size[0] = v->size[0];
  K->size[1] = v->size[0];
  emxEnsureCapacity_real_T(K, i);
  K_data = K->data;
  loop_ub_tmp = v->size[0] * v->size[0];
  for (i = 0; i < loop_ub_tmp; i++) {
    K_data[i] = 0.0;
  }
  for (j = 0; j < ibcol; j++) {
    K_data[j + K->size[0] * j] = Ai_data[j];
  }
  i = b->size[0];
  b->size[0] = ibtile;
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  for (itilerow = 0; itilerow < jp1j; itilerow++) {
    ibcol = itilerow * 3;
    b_data[ibcol] = bb;
    b_data[ibcol + 1] = bb;
    b_data[ibcol + 2] = bs;
  }
  i = v->size[0];
  v->size[0] = b->size[0] + 1;
  emxEnsureCapacity_real_T(v, i);
  Ai_data = v->data;
  Ai_data[0] = bm;
  ibcol = b->size[0];
  for (i = 0; i < ibcol; i++) {
    Ai_data[i + 1] = b_data[i];
  }
  ibcol = v->size[0];
  emxInit_real_T(&D, 2);
  i = D->size[0] * D->size[1];
  D->size[0] = v->size[0];
  D->size[1] = v->size[0];
  emxEnsureCapacity_real_T(D, i);
  D_data = D->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    D_data[i] = 0.0;
  }
  for (j = 0; j < ibcol; j++) {
    D_data[j + D->size[0] * j] = Ai_data[j];
  }
  /*  Fci = cell(1,4); */
  /*  Ai = cell(1,4); */
  emxInit_real_T(&Ai, 3);
  i = Ai->size[0] * Ai->size[1] * Ai->size[2];
  Ai->size[0] = 3;
  Ai->size[1] = 3;
  Ai->size[2] = (int)(N - 1.0);
  emxEnsureCapacity_real_T(Ai, i);
  Ai_data = Ai->data;
  ibcol = 9 * (int)(N - 1.0);
  for (i = 0; i < ibcol; i++) {
    Ai_data[i] = 0.0;
  }
  i = cq->size[0];
  cq->size[0] = (int)(N - 1.0);
  emxEnsureCapacity_real_T(cq, i);
  cq_data = cq->data;
  for (i = 0; i < jp1j; i++) {
    cq_data[i] = 0.0;
  }
  for (b_i = 0; b_i < jp1j; b_i++) {
    double Aq[9];
    double b_Aq[9];
    double Dq;
    double Dq_tmp;
    double L_tmp;
    double del;
    double del_tmp;
    double dy;
    dx = qp_data[3 * b_i];
    dy = qp_data[3 * b_i + 1];
    del_tmp = dx * dx;
    b_del_tmp = dy * dy;
    del = sqrt(del_tmp + b_del_tmp);
    theta = del / d;
    Dq_tmp = sin(del);
    Dq = del - Dq_tmp;
    L_tmp = qp_data[3 * b_i + 2] + L0;
    if ((dx < 1.0E-5) && (dy < 1.0E-5)) {
      cq_data[b_i] = L_tmp / 3.0;
      for (i = 0; i < 9; i++) {
        Aq[i] = iv[i];
      }
    } else {
      cq_data[b_i] = 2.0 * (L_tmp / theta - d) * sin(theta / 6.0);
      theta = rt_powd_snf(del, 3.0);
      Aq[0] = dx * dy * Dq / theta;
      Aq[3] = (-del_tmp * del - b_del_tmp * Dq_tmp) / theta;
      Aq[6] = dx * Dq * L_tmp / theta;
      Aq[1] = (b_del_tmp * del + del_tmp * Dq_tmp) / theta;
      Aq[4] = -dx * dy * Dq / theta;
      Aq[7] = dy * Dq * L_tmp / theta;
      Aq[2] = 0.0;
      Aq[5] = 0.0;
      Aq[8] = Dq_tmp / del;
    }
    dx = 60.0;
    b_cosd(&dx);
    b_del_tmp = 30.0;
    b_cosd(&b_del_tmp);
    for (i = 0; i < 3; i++) {
      theta = Aq[i];
      dy = Aq[i + 3];
      del_tmp = Aq[i + 6];
      for (k = 0; k < 3; k++) {
        b_Aq[i + 3 * k] =
            (theta * (double)c_b[3 * k] + dy * (double)c_b[3 * k + 1]) +
            del_tmp * (double)c_b[3 * k + 2];
      }
    }
    theta = d * dx;
    Aq[0] = theta;
    Aq[3] = theta;
    Aq[6] = -d;
    Aq[1] = -d * b_del_tmp;
    Aq[4] = d * b_del_tmp;
    Aq[7] = 0.0;
    Aq[2] = 1.0;
    Aq[5] = 1.0;
    Aq[8] = 1.0;
    for (i = 0; i < 3; i++) {
      theta = b_Aq[i];
      dy = b_Aq[i + 3];
      del_tmp = b_Aq[i + 6];
      for (k = 0; k < 3; k++) {
        Ai_data[(i + 3 * k) + 9 * b_i] =
            (theta * Aq[3 * k] + dy * Aq[3 * k + 1]) + del_tmp * Aq[3 * k + 2];
      }
    }
  }
  emxFree_real_T(&qp);
  /*  A = blkdiag(1,Ai(:,:,1),Ai(:,:,2),Ai(:,:,3)); */
  i = (int)((N - 1.0) * 3.0 + 1.0);
  k = A->size[0] * A->size[1];
  A->size[0] = i;
  A->size[1] = i;
  emxEnsureCapacity_real_T(A, k);
  qp_data = A->data;
  itilerow = i * i;
  for (i = 0; i < itilerow; i++) {
    qp_data[i] = 0.0;
  }
  qp_data[0] = 1.0;
  for (b_i = 0; b_i < jp1j; b_i++) {
    double a[3];
    theta = 3.0 * (((double)b_i + 1.0) - 1.0) + 2.0;
    a[0] = theta;
    a[1] = theta + 1.0;
    a[2] = theta + 2.0;
    for (i = 0; i < 3; i++) {
      ibcol = (int)a[i] - 1;
      k = 3 * i + 9 * b_i;
      qp_data[((int)theta + A->size[0] * ibcol) - 1] = Ai_data[k];
      qp_data[((int)(theta + 1.0) + A->size[0] * ibcol) - 1] = Ai_data[k + 1];
      qp_data[((int)(theta + 2.0) + A->size[0] * ibcol) - 1] = Ai_data[k + 2];
    }
  }
  emxFree_real_T(&Ai);
  emxInit_real_T(&b_b, 2);
  i = b_b->size[0] * b_b->size[1];
  b_b->size[0] = 1;
  b_b->size[1] = ibtile;
  emxEnsureCapacity_real_T(b_b, i);
  b_data = b_b->data;
  for (ibcol = 0; ibcol < jp1j; ibcol++) {
    ibtile = ibcol * 3;
    b_data[ibtile] = conv_pcc;
    b_data[ibtile + 1] = conv_pcc;
    b_data[ibtile + 2] = conv_pcc;
  }
  emxInit_real_T(&b_v, 2);
  i = b_v->size[0] * b_v->size[1];
  b_v->size[0] = 1;
  b_v->size[1] = b_b->size[1] + 1;
  emxEnsureCapacity_real_T(b_v, i);
  Ai_data = b_v->data;
  Ai_data[0] = conv_motor;
  ibcol = b_b->size[1];
  for (i = 0; i < ibcol; i++) {
    Ai_data[i + 1] = b_data[i];
  }
  emxFree_real_T(&b_b);
  ibcol = b_v->size[1];
  emxInit_real_T(&conversion, 2);
  i = conversion->size[0] * conversion->size[1];
  conversion->size[0] = b_v->size[1];
  conversion->size[1] = b_v->size[1];
  emxEnsureCapacity_real_T(conversion, i);
  conversion_data = conversion->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    conversion_data[i] = 0.0;
  }
  for (j = 0; j < ibcol; j++) {
    conversion_data[j + conversion->size[0] * j] = Ai_data[j];
  }
  emxFree_real_T(&b_v);
  ibcol = M->size[0] - 1;
  jp1j = M->size[1];
  i = b->size[0];
  b->size[0] = M->size[0];
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  for (b_i = 0; b_i <= ibcol; b_i++) {
    b_data[b_i] = 0.0;
  }
  for (k = 0; k < jp1j; k++) {
    ibtile = k * M->size[0];
    for (b_i = 0; b_i <= ibcol; b_i++) {
      b_data[b_i] += M_data[ibtile + b_i] * ddqd[k];
    }
  }
  mc_tmp = K->size[0] - 1;
  jp1j = K->size[1];
  i = v->size[0];
  v->size[0] = K->size[0];
  emxEnsureCapacity_real_T(v, i);
  Ai_data = v->data;
  for (b_i = 0; b_i <= mc_tmp; b_i++) {
    Ai_data[b_i] = 0.0;
  }
  for (k = 0; k < jp1j; k++) {
    ibtile = k * K->size[0];
    for (b_i = 0; b_i <= mc_tmp; b_i++) {
      Ai_data[b_i] += K_data[ibtile + b_i] * qd[k];
    }
  }
  jp1j = D->size[1];
  emxInit_real_T(&b_C, 1);
  i = b_C->size[0];
  b_C->size[0] = D->size[0];
  emxEnsureCapacity_real_T(b_C, i);
  cq_data = b_C->data;
  for (b_i = 0; b_i <= mc_tmp; b_i++) {
    cq_data[b_i] = 0.0;
  }
  for (k = 0; k < jp1j; k++) {
    ibtile = k * D->size[0];
    for (b_i = 0; b_i <= mc_tmp; b_i++) {
      cq_data[b_i] += D_data[ibtile + b_i] * dqd[k];
    }
  }
  emxFree_real_T(&D);
  if (C->size[0] == 1) {
    i = b->size[0];
  } else {
    i = C->size[0];
  }
  if (i == 1) {
    k = v->size[0];
  } else {
    k = i;
  }
  if (k == 1) {
    loop_ub_tmp = b_C->size[0];
  } else {
    loop_ub_tmp = k;
  }
  if ((C->size[0] == b->size[0]) && (i == v->size[0]) && (k == b_C->size[0]) &&
      (loop_ub_tmp == 10)) {
    double b_dqd[10];
    double b_qd[10];
    for (i = 0; i < 10; i++) {
      b_qd[i] = qd[i] - q[i];
      b_dqd[i] = dqd[i] - dq[i];
    }
    for (i = 0; i < 10; i++) {
      theta = 0.0;
      for (k = 0; k < 10; k++) {
        theta += Kp[i + 10 * k] * b_qd[k];
      }
      c_C[i] = (((C_data[i] + b_data[i]) + Ai_data[i]) + cq_data[i]) + theta;
    }
    for (i = 0; i < 10; i++) {
      theta = 0.0;
      for (k = 0; k < 10; k++) {
        theta += KD[i + 10 * k] * b_dqd[k];
      }
      c_C[i] += theta;
    }
  } else {
    binary_expand_op(c_C, C, b, v, b_C, Kp, qd, q, KD, dqd, dq);
  }
  emxFree_real_T(&b_C);
  emxFree_real_T(&v);
  if (A->size[1] == 0) {
    b->size[0] = 0;
  } else if (A->size[1] == 10) {
    signed char ipiv[10];
    i = K->size[0] * K->size[1];
    K->size[0] = A->size[0];
    K->size[1] = 10;
    emxEnsureCapacity_real_T(K, i);
    K_data = K->data;
    for (i = 0; i < 100; i++) {
      K_data[i] = qp_data[i];
    }
    for (i = 0; i < 10; i++) {
      ipiv[i] = (signed char)(i + 1);
    }
    for (j = 0; j < 9; j++) {
      signed char i1;
      loop_ub_tmp = 8 - j;
      itilerow = j * 11;
      jp1j = itilerow + 2;
      ibcol = 10 - j;
      ibtile = 0;
      theta = fabs(K_data[itilerow]);
      for (k = 2; k <= ibcol; k++) {
        dx = fabs(K_data[(itilerow + k) - 1]);
        if (dx > theta) {
          ibtile = k - 1;
          theta = dx;
        }
      }
      if (K_data[itilerow + ibtile] != 0.0) {
        if (ibtile != 0) {
          ibcol = j + ibtile;
          ipiv[j] = (signed char)(ibcol + 1);
          for (k = 0; k < 10; k++) {
            ibtile = j + k * 10;
            theta = K_data[ibtile];
            i = ibcol + k * 10;
            K_data[ibtile] = K_data[i];
            K_data[i] = theta;
          }
        }
        i = (itilerow - j) + 10;
        for (b_i = jp1j; b_i <= i; b_i++) {
          K_data[b_i - 1] /= K_data[itilerow];
        }
      }
      ibcol = itilerow;
      for (ibtile = 0; ibtile <= loop_ub_tmp; ibtile++) {
        theta = K_data[(itilerow + ibtile * 10) + 10];
        if (theta != 0.0) {
          i = ibcol + 12;
          k = (ibcol - j) + 20;
          for (jp1j = i; jp1j <= k; jp1j++) {
            K_data[jp1j - 1] +=
                K_data[((itilerow + jp1j) - ibcol) - 11] * -theta;
          }
        }
        ibcol += 10;
      }
      i1 = ipiv[j];
      if (i1 != j + 1) {
        theta = c_C[j];
        c_C[j] = c_C[i1 - 1];
        c_C[i1 - 1] = theta;
      }
    }
    for (k = 0; k < 10; k++) {
      ibcol = 10 * k;
      if (c_C[k] != 0.0) {
        i = k + 2;
        for (b_i = i; b_i < 11; b_i++) {
          c_C[b_i - 1] -= c_C[k] * K_data[(b_i + ibcol) - 1];
        }
      }
    }
    for (k = 9; k >= 0; k--) {
      ibcol = 10 * k;
      theta = c_C[k];
      if (theta != 0.0) {
        theta /= K_data[k + ibcol];
        c_C[k] = theta;
        for (b_i = 0; b_i < k; b_i++) {
          c_C[b_i] -= c_C[k] * K_data[b_i + ibcol];
        }
      }
    }
    i = b->size[0];
    b->size[0] = 10;
    emxEnsureCapacity_real_T(b, i);
    b_data = b->data;
    for (i = 0; i < 10; i++) {
      b_data[i] = c_C[i];
    }
  } else {
    i = K->size[0] * K->size[1];
    K->size[0] = A->size[0];
    K->size[1] = A->size[1];
    emxEnsureCapacity_real_T(K, i);
    K_data = K->data;
    for (i = 0; i < itilerow; i++) {
      K_data[i] = qp_data[i];
    }
    emxInit_int32_T(&jpvt, 2);
    xgeqp3(K, tau, jpvt);
    jpvt_data = jpvt->data;
    cq_data = tau->data;
    K_data = K->data;
    jp1j = 0;
    if (K->size[0] < K->size[1]) {
      ibcol = K->size[0];
      ibtile = K->size[1];
    } else {
      ibcol = K->size[1];
      ibtile = K->size[0];
    }
    if (ibcol > 0) {
      theta =
          fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)ibtile) *
          fabs(K_data[0]);
      while ((jp1j < ibcol) &&
             (!(fabs(K_data[jp1j + K->size[0] * jp1j]) <= theta))) {
        jp1j++;
      }
    }
    i = b->size[0];
    b->size[0] = K->size[1];
    emxEnsureCapacity_real_T(b, i);
    b_data = b->data;
    ibcol = K->size[1];
    for (i = 0; i < ibcol; i++) {
      b_data[i] = 0.0;
    }
    if (K->size[1] > 10) {
      i = 10;
    } else {
      i = K->size[1];
    }
    for (j = 0; j < i; j++) {
      if (cq_data[j] != 0.0) {
        theta = c_C[j];
        k = j + 2;
        for (b_i = k; b_i < 11; b_i++) {
          theta += K_data[(b_i + K->size[0] * j) - 1] * c_C[b_i - 1];
        }
        theta *= cq_data[j];
        if (theta != 0.0) {
          c_C[j] -= theta;
          for (b_i = k; b_i < 11; b_i++) {
            c_C[b_i - 1] -= K_data[(b_i + K->size[0] * j) - 1] * theta;
          }
        }
      }
    }
    i = (unsigned char)jp1j;
    for (b_i = 0; b_i < i; b_i++) {
      b_data[jpvt_data[b_i] - 1] = c_C[b_i];
    }
    for (j = jp1j; j >= 1; j--) {
      i = jpvt_data[j - 1];
      b_data[i - 1] /= K_data[(j + K->size[0] * (j - 1)) - 1];
      k = (unsigned char)(j - 1);
      for (b_i = 0; b_i < k; b_i++) {
        b_data[jpvt_data[b_i] - 1] -=
            b_data[i - 1] * K_data[b_i + K->size[0] * (j - 1)];
      }
    }
    emxFree_int32_T(&jpvt);
  }
  emxFree_real_T(&K);
  jp1j = conversion->size[1];
  i = tau->size[0];
  tau->size[0] = conversion->size[0];
  emxEnsureCapacity_real_T(tau, i);
  cq_data = tau->data;
  for (b_i = 0; b_i <= mc_tmp; b_i++) {
    cq_data[b_i] = 0.0;
  }
  for (k = 0; k < jp1j; k++) {
    ibtile = k * conversion->size[0];
    for (b_i = 0; b_i <= mc_tmp; b_i++) {
      cq_data[b_i] += conversion_data[ibtile + b_i] * b_data[k];
    }
  }
  emxFree_real_T(&b);
  emxFree_real_T(&conversion);
  /*  [J, x] = J_r(q,L0,d); */
  /*  tau_r = conversion*(A\(C + J'*Jbar'*(K*q + D*dq) + J'*Lam*(Kpx*(xd-x) +
   * KDx*(dxd - dxr))) + (eye(10) - J'*Jbar')*(-Kpr * q - KDr * dq)); */
  *tau_r = 0.0;
}

/*
 * File trailer for helix_controller.c
 *
 * [EOF]
 */
