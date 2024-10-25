/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: helix_controller.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 23-Oct-2024 10:17:18
 */

/* Include Files */
#include "helix_controller.h"
#include "J_r.h"
#include "blkdiag.h"
#include "diag.h"
#include "helix_controller_rtwutil.h"
#include "mldivide.h"
#include "mrdivide_helper.h"
#include "mtimes.h"
#include "repmat.h"
#include "rt_nonfinite.h"
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
  static const double b_A[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  static const double b_b[6] = {-0.0, -0.0, -9.81, -0.0, -0.0, -0.0};
  static const signed char b[36] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                    0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  static const signed char c_b[9] = {-1, 0, 0, 0, -1, 0, 0, 0, -1};
  static const signed char iv1[9] = {0, 1, 0, -1, 0, 0, 0, 0, 1};
  static const signed char b_a[6] = {0, 0, 0, 0, 0, 1};
  static const signed char d_a[6] = {0, 0, 0, 0, 0, 2};
  static const signed char iv[6] = {0, 0, 0, 0, 0, 1};
  double Xup[144];
  double b_I[144];
  double D[100];
  double K[100];
  double KDr[100];
  double Kpr[100];
  double conversion[100];
  double y_tmp[100];
  double b_Smod[54];
  double a_tmp[36];
  double c_I[36];
  double J[30];
  double b_b_tmp[30];
  double b_tmp[30];
  double b_tmp_tmp[30];
  double Ai[27];
  double a[24];
  double f[24];
  double v[24];
  double c_Xup[18];
  double fh3[18];
  double g[16];
  double b_C[10];
  double b_K[10];
  double b_bm[10];
  double b_v[10];
  double Smod[9];
  double b_skw[9];
  double c_skw[9];
  double qp[9];
  double skw[9];
  double d_I[6];
  double vJ[6];
  double b_dq[3];
  double Dq;
  double Smod_tmp_tmp;
  double b_del_tmp;
  double c_a;
  double del;
  double del_tmp;
  double f_fh3_tmp;
  double fh3_tmp;
  double g_tmp;
  double theta;
  double vJ_tmp;
  int I_tmp;
  int b_i;
  int i;
  int i1;
  int i2;
  int ibtile;
  int j;
  signed char e_I[100];
  signed char Xtree[36];
  signed char qi_data[9];
  signed char tmp_data[9];
  qp[0] = q[1];
  qp[3] = q[4];
  qp[6] = q[7];
  qp[1] = q[2];
  qp[4] = q[5];
  qp[7] = q[8];
  qp[2] = q[3];
  qp[5] = q[6];
  qp[8] = q[9];
  /*  q = state */
  /*  qd = velocity */
  /*  qdd = acceleration */
  /*  m = mass of a module; */
  /*  mm = mass of a motor; */
  /*  hm = height of a motor; */
  /*  rm = "radius" a motor; */
  /*  r = radius of the hex plate; */
  /*  L0 = Length of a module; */
  /*  d = distance to cable; */
  /*  N = number of links (number of modules plus number of motors) */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  This function calculates the Mass and Coriolis + Gravity Matrices */
  /*  given the current state and geometric parameters of the pushpuppet robot
   */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  qi = {1, 2:4, 5, 6:8, 9, 10:12, 13, 14:16}; */
  memset(&b_I[0], 0, 144U * sizeof(double));
  vJ[0] = m;
  vJ[1] = m;
  vJ[2] = m;
  theta = r * r;
  vJ_tmp = 0.25 * m * theta;
  vJ[3] = vJ_tmp;
  vJ[4] = vJ_tmp;
  vJ[5] = 0.5 * m * theta;
  memset(&b_I[0], 0, 36U * sizeof(double));
  for (j = 0; j < 6; j++) {
    b_I[j + 6 * j] = vJ[j];
  }
  for (i = 0; i < 36; i++) {
    fh3_tmp = b_I[i];
    b_I[i + 36] = fh3_tmp;
    b_I[i + 72] = fh3_tmp;
    b_I[i + 108] = fh3_tmp;
  }
  memset(&Xup[0], 0, 144U * sizeof(double));
  memset(&v[0], 0, 24U * sizeof(double));
  memset(&a[0], 0, 24U * sizeof(double));
  memset(&f[0], 0, 24U * sizeof(double));
  memset(&C[0], 0, 10U * sizeof(double));
  for (i = 0; i < 36; i++) {
    Xtree[i] = b[i];
  }
  /* [diag(ones(3,1)) [0;0;0]; 0 0 0 1]; */
  Dq = sin(q[0]);
  g_tmp = cos(q[0]);
  g[0] = g_tmp;
  g[4] = -Dq;
  g[8] = 0.0;
  g[12] = 0.0;
  g[1] = Dq;
  g[5] = g_tmp;
  g[9] = 0.0;
  g[13] = 0.0;
  g[2] = 0.0;
  g[3] = 0.0;
  g[6] = 0.0;
  g[7] = 0.0;
  g[10] = 1.0;
  g[11] = 0.0;
  g[14] = 0.0;
  g[15] = 1.0;
  /*  convert transform to rotation and translation */
  /*  Adjoint calculator */
  /*  This function calculates the adjoint given simply the transform between */
  /*  joints */
  /*  get skew symmetric matrix of translation */
  for (i = 0; i < 3; i++) {
    skw[3 * i] = g[i];
    skw[3 * i + 1] = g[i + 4];
    skw[3 * i + 2] = g[i + 8];
  }
  for (i = 0; i < 9; i++) {
    b_skw[i] = -skw[i];
  }
  Smod[0] = 0.0;
  Smod[3] = -0.0;
  Smod[6] = 0.0;
  Smod[1] = 0.0;
  Smod[4] = 0.0;
  Smod[7] = -0.0;
  Smod[2] = -0.0;
  Smod[5] = 0.0;
  Smod[8] = 0.0;
  for (i = 0; i < 3; i++) {
    fh3_tmp = b_skw[i];
    Smod_tmp_tmp = b_skw[i + 3];
    theta = b_skw[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      c_skw[i + 3 * i1] =
          (fh3_tmp * Smod[3 * i1] + Smod_tmp_tmp * Smod[3 * i1 + 1]) +
          theta * Smod[3 * i1 + 2];
      c_I[i1 + 6 * i] = skw[i1 + 3 * i];
    }
  }
  for (i = 0; i < 3; i++) {
    I_tmp = 6 * (i + 3);
    c_I[I_tmp] = c_skw[3 * i];
    c_I[6 * i + 3] = 0.0;
    c_I[I_tmp + 3] = skw[3 * i];
    ibtile = 3 * i + 1;
    c_I[I_tmp + 1] = c_skw[ibtile];
    c_I[6 * i + 4] = 0.0;
    c_I[I_tmp + 4] = skw[ibtile];
    ibtile = 3 * i + 2;
    c_I[I_tmp + 2] = c_skw[ibtile];
    c_I[6 * i + 5] = 0.0;
    c_I[I_tmp + 5] = skw[ibtile];
  }
  for (b_i = 0; b_i < 6; b_i++) {
    for (i = 0; i < 6; i++) {
      fh3_tmp = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        fh3_tmp += c_I[b_i + 6 * i1] * (double)b[i1 + 6 * i];
      }
      Xup[b_i + 6 * i] = fh3_tmp;
    }
    fh3_tmp = (double)b_a[b_i] * dq[0];
    vJ[b_i] = fh3_tmp;
    v[b_i] = fh3_tmp;
  }
  memset(&c_I[0], 0, 36U * sizeof(double));
  skw[0] = 0.0;
  skw[3] = -v[5];
  skw[6] = v[4];
  skw[1] = v[5];
  skw[4] = 0.0;
  skw[7] = -v[3];
  skw[2] = -v[4];
  skw[5] = v[3];
  skw[8] = 0.0;
  for (i = 0; i < 3; i++) {
    theta = skw[3 * i];
    c_I[6 * i] = theta;
    I_tmp = 6 * (i + 3);
    c_I[I_tmp + 3] = theta;
    theta = skw[3 * i + 1];
    c_I[6 * i + 1] = theta;
    c_I[I_tmp + 4] = theta;
    theta = skw[3 * i + 2];
    c_I[6 * i + 2] = theta;
    c_I[I_tmp + 5] = theta;
  }
  c_I[18] = 0.0;
  c_I[24] = -v[2];
  c_I[30] = v[1];
  c_I[19] = v[2];
  c_I[25] = 0.0;
  c_I[31] = -v[0];
  c_I[20] = -v[1];
  c_I[26] = v[0];
  c_I[32] = 0.0;
  for (i = 0; i < 6; i++) {
    fh3_tmp = 0.0;
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      i2 = i + 6 * i1;
      fh3_tmp += Xup[i2] * b_b[i1];
      theta = c_I[i2];
      Smod_tmp_tmp += theta * vJ[i1];
      a_tmp[i1 + 6 * i] = -theta;
    }
    a[i] = fh3_tmp + Smod_tmp_tmp;
  }
  for (i = 0; i < 6; i++) {
    d_I[i] = 0.0;
    vJ[i] = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      fh3_tmp = 0.0;
      for (i2 = 0; i2 < 6; i2++) {
        fh3_tmp += a_tmp[i + 6 * i2] * b_I[i2 + 6 * i1];
      }
      d_I[i] += b_I[i + 6 * i1] * a[i1];
      vJ[i] += fh3_tmp * v[i1];
    }
    f[i] = d_I[i] + vJ[i];
  }
  /*  Recursive Newton Euler to Calculate C+G */
  Smod[0] = 0.0;
  Smod[4] = 0.0;
  Smod[8] = 0.0;
  for (b_i = 0; b_i < 3; b_i++) {
    double b_Xup[6];
    double c_del_tmp;
    double d_del_tmp;
    double e_del_tmp;
    double f_del_tmp;
    I_tmp = b_i * 3 + 2;
    ibtile = b_i * 3 + 4;
    if (ibtile < I_tmp) {
    } else {
      ibtile -= I_tmp;
      for (i = 0; i <= ibtile; i++) {
        qi_data[i] = (signed char)(I_tmp + i);
      }
    }
    del_tmp = q[qi_data[0] - 1];
    b_del_tmp = q[qi_data[1] - 1];
    c_del_tmp = del_tmp * del_tmp;
    d_del_tmp = b_del_tmp * b_del_tmp;
    e_del_tmp = c_del_tmp + d_del_tmp;
    f_del_tmp = sqrt(e_del_tmp);
    if (f_del_tmp < 1.0E-6) {
      memset(&skw[0], 0, 9U * sizeof(double));
      skw[0] = 1.0;
      skw[4] = 1.0;
      skw[8] = 1.0;
      for (i = 0; i < 3; i++) {
        I_tmp = i << 2;
        g[I_tmp] = skw[3 * i];
        g[I_tmp + 1] = skw[3 * i + 1];
        g[I_tmp + 2] = skw[3 * i + 2];
      }
      double Smod_tmp;
      double b_fh3_tmp;
      double c_fh3_tmp;
      double d_fh3_tmp;
      double fh3_tmp_tmp;
      g[12] = 0.0;
      g[13] = 0.0;
      del = q[qi_data[2] - 1];
      Dq = L0 + del;
      g[14] = Dq;
      g[3] = 0.0;
      g[7] = 0.0;
      g[11] = 0.0;
      g[15] = 1.0;
      Smod_tmp = Dq / (2.0 * d);
      b_Smod[18 * b_i] = Smod_tmp;
      b_Smod[18 * b_i + 6] = 0.0;
      b_Smod[18 * b_i + 12] = 0.0;
      b_Smod[18 * b_i + 1] = 0.0;
      b_Smod[18 * b_i + 7] = Smod_tmp;
      b_Smod[18 * b_i + 13] = 0.0;
      b_Smod[18 * b_i + 3] = 0.0;
      b_Smod[18 * b_i + 9] = -1.0 / d;
      b_Smod[18 * b_i + 15] = 0.0;
      b_Smod[18 * b_i + 4] = 1.0 / d;
      b_Smod[18 * b_i + 10] = 0.0;
      b_Smod[18 * b_i + 16] = 0.0;
      b_fh3_tmp = dq[qi_data[2] - 1] / (2.0 * d);
      fh3[0] = b_fh3_tmp;
      fh3[6] = 0.0;
      c_fh3_tmp = dq[qi_data[0] - 1];
      fh3[12] = -c_fh3_tmp / (2.0 * d);
      fh3[1] = 0.0;
      fh3[7] = b_fh3_tmp;
      b_fh3_tmp = dq[qi_data[1] - 1];
      fh3[13] = -b_fh3_tmp / (2.0 * d);
      fh3_tmp_tmp = d * d;
      d_fh3_tmp = 6.0 * fh3_tmp_tmp;
      fh3[2] = c_fh3_tmp * Dq / d_fh3_tmp;
      fh3[8] = (L0 * b_fh3_tmp + b_fh3_tmp * del) / d_fh3_tmp;
      fh3[14] = 0.0;
      b_Smod[18 * b_i + 2] = 0.0;
      b_Smod[18 * b_i + 5] = 0.0;
      fh3[3] = 0.0;
      fh3[4] = 0.0;
      b_Smod[18 * b_i + 8] = 0.0;
      b_Smod[18 * b_i + 11] = 0.0;
      fh3[9] = 0.0;
      fh3[10] = 0.0;
      b_Smod[18 * b_i + 14] = 1.0;
      b_Smod[18 * b_i + 17] = 0.0;
      fh3[15] = 0.0;
      fh3[16] = 0.0;
      d_fh3_tmp = 2.0 * fh3_tmp_tmp;
      fh3[5] = b_fh3_tmp / d_fh3_tmp;
      fh3[11] = -c_fh3_tmp / d_fh3_tmp;
      fh3[17] = 0.0;
    } else {
      double Smod_tmp;
      double Smod_tmp_tmp_tmp;
      double b_Smod_tmp;
      double b_Smod_tmp_tmp;
      double b_a_tmp;
      double b_fh3_tmp;
      double b_fh3_tmp_tmp;
      double b_fh3_tmp_tmp_tmp;
      double b_g_tmp_tmp;
      double c_Smod_tmp;
      double c_Smod_tmp_tmp;
      double c_fh3_tmp;
      double c_fh3_tmp_tmp;
      double d_Smod_tmp;
      double d_Smod_tmp_tmp;
      double d_fh3_tmp;
      double d_fh3_tmp_tmp;
      double e_Smod_tmp;
      double e_Smod_tmp_tmp;
      double e_fh3_tmp;
      double e_fh3_tmp_tmp;
      double f_Smod_tmp;
      double f_Smod_tmp_tmp;
      double fh3_tmp_tmp;
      double fh3_tmp_tmp_tmp;
      double g_Smod_tmp;
      double g_Smod_tmp_tmp;
      double g_fh3_tmp;
      double g_tmp_tmp;
      double h_Smod_tmp;
      double h_Smod_tmp_tmp;
      double h_fh3_tmp;
      double i_Smod_tmp;
      double i_fh3_tmp;
      double j_Smod_tmp;
      double j_fh3_tmp;
      double k_Smod_tmp;
      double k_fh3_tmp;
      double l_Smod_tmp;
      double l_fh3_tmp;
      double m_Smod_tmp;
      double m_fh3_tmp;
      double n_Smod_tmp;
      double n_fh3_tmp;
      double o_Smod_tmp;
      double o_fh3_tmp;
      double p_Smod_tmp;
      double p_fh3_tmp;
      double q_Smod_tmp;
      double r_Smod_tmp;
      vJ_tmp = f_del_tmp * f_del_tmp;
      b_a_tmp = L0 + q[qi_data[2] - 1];
      theta = d * b_a_tmp;
      c_a = theta / vJ_tmp;
      del = f_del_tmp / d;
      g_tmp_tmp = cos(del);
      g[0] = c_del_tmp / vJ_tmp * (g_tmp_tmp - 1.0) + 1.0;
      b_g_tmp_tmp = del_tmp * b_del_tmp;
      Dq = b_g_tmp_tmp / vJ_tmp * (g_tmp_tmp - 1.0);
      g[1] = Dq;
      g_tmp = sin(del);
      g[2] = -del_tmp / f_del_tmp * g_tmp;
      g[4] = Dq;
      g[5] = d_del_tmp / vJ_tmp * (g_tmp_tmp - 1.0) + 1.0;
      g[6] = -b_del_tmp / f_del_tmp * g_tmp;
      g[8] = del_tmp / f_del_tmp * g_tmp;
      g[9] = b_del_tmp / f_del_tmp * g_tmp;
      g[10] = g_tmp_tmp;
      g[12] = c_a * (del_tmp * (1.0 - g_tmp_tmp));
      g[13] = c_a * (b_del_tmp * (1.0 - g_tmp_tmp));
      Dq = f_del_tmp * g_tmp;
      g[14] = c_a * Dq;
      g[3] = 0.0;
      g[7] = 0.0;
      g[11] = 0.0;
      g[15] = 1.0;
      b_Smod_tmp_tmp = c_del_tmp * g_tmp;
      Smod_tmp = b_Smod_tmp_tmp * b_a_tmp;
      b_Smod_tmp = e_del_tmp * e_del_tmp;
      c_Smod_tmp = rt_powd_snf(e_del_tmp, 1.5);
      vJ_tmp = theta * (g_tmp_tmp - 1.0);
      d_Smod_tmp = vJ_tmp / e_del_tmp;
      c_Smod_tmp_tmp = d * g_tmp;
      e_Smod_tmp = g_tmp_tmp * f_del_tmp - c_Smod_tmp_tmp;
      d_Smod_tmp_tmp = rt_powd_snf(e_del_tmp, 3.0);
      Smod_tmp_tmp_tmp = c_del_tmp * d_del_tmp;
      theta = Smod_tmp_tmp_tmp * b_a_tmp * (g_tmp_tmp - 1.0);
      Smod_tmp_tmp = (Dq - 2.0 * d) + 2.0 * d * g_tmp_tmp;
      f_Smod_tmp = theta * Smod_tmp_tmp / d_Smod_tmp_tmp;
      e_Smod_tmp_tmp = c_del_tmp * g_tmp_tmp;
      g_Smod_tmp = d_del_tmp + e_Smod_tmp_tmp;
      h_Smod_tmp =
          (Smod_tmp / c_Smod_tmp - d_Smod_tmp) +
          2.0 * d * c_del_tmp * b_a_tmp * (g_tmp_tmp - 1.0) / b_Smod_tmp;
      b_Smod[18 * b_i] = (g_Smod_tmp * h_Smod_tmp / e_del_tmp -
                          Smod_tmp * e_Smod_tmp / b_Smod_tmp) +
                         f_Smod_tmp;
      b_Smod[18 * b_i + 6] = 0.0;
      Smod_tmp = d * del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * b_i + 12] = Smod_tmp / e_del_tmp;
      b_Smod[18 * b_i + 1] = 0.0;
      f_Smod_tmp_tmp = d_del_tmp * g_tmp;
      i_Smod_tmp = f_Smod_tmp_tmp * b_a_tmp;
      g_Smod_tmp_tmp = d_del_tmp * g_tmp_tmp;
      j_Smod_tmp = c_del_tmp + g_Smod_tmp_tmp;
      d_Smod_tmp =
          (i_Smod_tmp / c_Smod_tmp - d_Smod_tmp) +
          2.0 * d * d_del_tmp * b_a_tmp * (g_tmp_tmp - 1.0) / b_Smod_tmp;
      b_Smod[18 * b_i + 7] = (j_Smod_tmp * d_Smod_tmp / e_del_tmp -
                              i_Smod_tmp * e_Smod_tmp / b_Smod_tmp) +
                             f_Smod_tmp;
      f_Smod_tmp = d * b_del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * b_i + 13] = f_Smod_tmp / e_del_tmp;
      i_Smod_tmp = c_Smod_tmp - c_Smod_tmp_tmp * e_del_tmp;
      k_Smod_tmp = rt_powd_snf(e_del_tmp, 2.5);
      l_Smod_tmp = del_tmp * b_a_tmp;
      b_Smod[18 * b_i + 2] = l_Smod_tmp * i_Smod_tmp / k_Smod_tmp;
      m_Smod_tmp = b_del_tmp * b_a_tmp;
      b_Smod[18 * b_i + 8] = m_Smod_tmp * i_Smod_tmp / k_Smod_tmp;
      b_Smod[18 * b_i + 14] = c_Smod_tmp_tmp / f_del_tmp;
      n_Smod_tmp = d * c_Smod_tmp;
      h_Smod_tmp_tmp = f_del_tmp - c_Smod_tmp_tmp;
      o_Smod_tmp = b_g_tmp_tmp * h_Smod_tmp_tmp;
      b_Smod[18 * b_i + 3] = -o_Smod_tmp / n_Smod_tmp;
      p_Smod_tmp = d_del_tmp * f_del_tmp + d * c_del_tmp * g_tmp;
      b_Smod[18 * b_i + 9] = -p_Smod_tmp / n_Smod_tmp;
      b_Smod[18 * b_i + 15] = 0.0;
      q_Smod_tmp = c_del_tmp * f_del_tmp + d * d_del_tmp * g_tmp;
      b_Smod[18 * b_i + 4] = q_Smod_tmp / n_Smod_tmp;
      b_Smod[18 * b_i + 10] = o_Smod_tmp / n_Smod_tmp;
      b_Smod[18 * b_i + 16] = 0.0;
      o_Smod_tmp = b_del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * b_i + 5] = -o_Smod_tmp / e_del_tmp;
      r_Smod_tmp = del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * b_i + 11] = r_Smod_tmp / e_del_tmp;
      b_Smod[18 * b_i + 17] = 0.0;
      /*      J = [ -(d*(L0 + dL)*(cos(del/d)*dx^2 + dy^2)*(cos(del/d) -
       * 1))/del^4,                                                -(d*dx*dy*(L0
       * + dL)*(cos(del/d) - 1)^2)/del^4, - (d*dx*sin(del/d)^2)/del^2 -
       * (d*dx*dy^2*(cos(del/d) - 1)^2)/del^4 - (d*dx*(cos(del/d)*dx^2 +
       * dy^2)*(cos(del/d) - 1))/del^4; */
      /*                                                                                  -(d*dx*dy*(L0
       * + dL)*(cos(del/d) - 1)^2)/del^4,                               -(d*(L0
       * + dL)*(dx^2 + cos(del/d)*dy^2)*(cos(del/d) - 1))/del^4,
       * (d*dy*sin(del/d)^2)/del^2 - (d*dx^2*dy*(cos(del/d) - 1)^2)/del^4 -
       * (d*dy*(dx^2 + cos(del/d)*dy^2)*(cos(del/d) - 1))/del^4; */
      /*                                                                            -(d*dx*sin(del/d)*(L0
       * + dL)*(cos(del/d) - 1))/del^3, (d*dy*sin(del/d)*(L0 + dL)*(cos(del/d) -
       * 1))/del^3, (d*sin(del/d)*(del^2*cos(del/d) - dx^2*cos(del/d) +
       * dy^2*cos(del/d) + dx^2 - dy^2))/del^3; */
      /*                                                                                      (dx*dy*sin(del/d)*(cos(del/d)
       * - 1))/del^3, -(sin(del/d)*(del^2*cos(del/d) - dx^2*cos(del/d) +
       * 2*dy^2*cos(del/d) + dx^2 - 2*dy^2))/del^3, 0; */
      /*                                                                                    (sin(del/d)*(cos(del/d)*dx^2
       * + dy^2))/del^3, (dx*dy*sin(del/d)*(cos(del/d) - 1))/del^3, 0; */
      /*          (dy*(dx^2 + cos(del/d)*dy^2)*(cos(del/d) - 1))/del^4 -
       * (dy*sin(del/d)^2)/del^2 + (2*dx^2*dy*(cos(del/d) - 1)^2)/del^4,
       * (dx*(dx^2 + cos(del/d)*dy^2)*(cos(del/d) - 1))/del^4, 0]; */
      b_fh3_tmp = dq[qi_data[2] - 1];
      fh3_tmp_tmp_tmp = dq[qi_data[0] - 1];
      fh3_tmp_tmp = 2.0 * del_tmp * fh3_tmp_tmp_tmp;
      b_fh3_tmp_tmp_tmp = dq[qi_data[1] - 1];
      b_fh3_tmp_tmp = 2.0 * b_del_tmp * b_fh3_tmp_tmp_tmp;
      c_fh3_tmp = fh3_tmp_tmp + b_fh3_tmp_tmp;
      c_fh3_tmp_tmp = b_fh3_tmp * c_del_tmp;
      d_fh3_tmp = c_fh3_tmp_tmp * g_tmp;
      e_fh3_tmp = 2.0 * d * b_Smod_tmp;
      g_fh3_tmp = fh3_tmp_tmp * g_tmp * b_a_tmp;
      h_fh3_tmp = e_Smod_tmp_tmp * b_a_tmp;
      i_fh3_tmp = 2.0 * c_del_tmp;
      j_fh3_tmp = 2.0 * c_Smod_tmp;
      k_fh3_tmp = g_tmp * b_a_tmp * c_fh3_tmp / j_fh3_tmp -
                  d * b_fh3_tmp * (g_tmp_tmp - 1.0) / e_del_tmp;
      l_fh3_tmp = 2.0 * k_Smod_tmp;
      m_fh3_tmp = 2.0 * d * b_fh3_tmp;
      c_a = vJ_tmp * c_fh3_tmp / b_Smod_tmp;
      n_fh3_tmp = 2.0 * d * f_del_tmp;
      o_fh3_tmp = g_tmp * g_tmp;
      d_fh3_tmp_tmp = g_tmp_tmp * c_fh3_tmp;
      e_fh3_tmp_tmp = 2.0 * f_del_tmp;
      theta = theta *
              (g_tmp * c_fh3_tmp / e_fh3_tmp_tmp - d_fh3_tmp_tmp / (2.0 * d)) /
              d_Smod_tmp_tmp;
      del = c_fh3_tmp_tmp * d_del_tmp * (g_tmp_tmp - 1.0) * Smod_tmp_tmp /
            d_Smod_tmp_tmp;
      p_fh3_tmp = 2.0 * d * k_Smod_tmp;
      Dq = 3.0 * c_del_tmp * d_del_tmp * b_a_tmp * (g_tmp_tmp - 1.0) *
           c_fh3_tmp * Smod_tmp_tmp / rt_powd_snf(e_del_tmp, 4.0);
      f_fh3_tmp = fh3_tmp_tmp * d_del_tmp * b_a_tmp * (g_tmp_tmp - 1.0) *
                  Smod_tmp_tmp / d_Smod_tmp_tmp;
      fh3_tmp = i_fh3_tmp * b_del_tmp * b_fh3_tmp_tmp_tmp * b_a_tmp *
                (g_tmp_tmp - 1.0) * Smod_tmp_tmp / d_Smod_tmp_tmp;
      c_fh3_tmp_tmp = rt_powd_snf(e_del_tmp, 3.5);
      vJ_tmp = Smod_tmp_tmp_tmp * g_tmp * b_a_tmp * c_fh3_tmp * Smod_tmp_tmp /
               (2.0 * d * c_fh3_tmp_tmp);
      fh3[0] =
          ((((((((((((g_Smod_tmp *
                          ((((((((k_fh3_tmp + d_fh3_tmp / c_Smod_tmp) -
                                 5.0 * c_del_tmp * g_tmp * b_a_tmp * c_fh3_tmp /
                                     l_fh3_tmp) +
                                g_fh3_tmp / c_Smod_tmp) +
                               m_fh3_tmp * c_del_tmp * (g_tmp_tmp - 1.0) /
                                   b_Smod_tmp) +
                              c_a) +
                             4.0 * d * del_tmp * fh3_tmp_tmp_tmp * b_a_tmp *
                                 (g_tmp_tmp - 1.0) / b_Smod_tmp) -
                            4.0 * d * c_del_tmp * b_a_tmp * (g_tmp_tmp - 1.0) *
                                c_fh3_tmp / d_Smod_tmp_tmp) +
                           h_fh3_tmp * c_fh3_tmp / e_fh3_tmp) /
                          e_del_tmp +
                      h_Smod_tmp *
                          ((b_fh3_tmp_tmp + fh3_tmp_tmp * g_tmp_tmp) -
                           b_Smod_tmp_tmp * c_fh3_tmp / n_fh3_tmp) /
                          e_del_tmp) -
                     g_Smod_tmp * c_fh3_tmp * h_Smod_tmp / b_Smod_tmp) -
                    d_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
                   c_del_tmp * o_fh3_tmp * b_a_tmp * c_fh3_tmp / e_fh3_tmp) +
                  i_fh3_tmp * g_tmp * b_a_tmp * e_Smod_tmp * c_fh3_tmp /
                      d_Smod_tmp_tmp) -
                 theta) -
                g_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
               del) -
              h_fh3_tmp * e_Smod_tmp * c_fh3_tmp / p_fh3_tmp) -
             Dq) +
            f_fh3_tmp) +
           fh3_tmp) -
          vJ_tmp;
      fh3[6] = 0.0;
      d_fh3_tmp = del_tmp * g_tmp * c_fh3_tmp;
      fh3[12] = (d * fh3_tmp_tmp_tmp * (g_tmp_tmp - 1.0) / e_del_tmp -
                 d_fh3_tmp / j_fh3_tmp) -
                Smod_tmp * c_fh3_tmp / b_Smod_tmp;
      fh3[1] = 0.0;
      g_fh3_tmp = b_fh3_tmp * d_del_tmp * g_tmp;
      h_fh3_tmp = b_fh3_tmp_tmp * g_tmp * b_a_tmp;
      i_fh3_tmp = g_Smod_tmp_tmp * b_a_tmp;
      fh3[7] =
          ((((((((((((j_Smod_tmp *
                          ((((((((k_fh3_tmp + g_fh3_tmp / c_Smod_tmp) -
                                 5.0 * d_del_tmp * g_tmp * b_a_tmp * c_fh3_tmp /
                                     l_fh3_tmp) +
                                h_fh3_tmp / c_Smod_tmp) +
                               m_fh3_tmp * d_del_tmp * (g_tmp_tmp - 1.0) /
                                   b_Smod_tmp) +
                              c_a) +
                             4.0 * d * b_del_tmp * b_fh3_tmp_tmp_tmp * b_a_tmp *
                                 (g_tmp_tmp - 1.0) / b_Smod_tmp) -
                            4.0 * d * d_del_tmp * b_a_tmp * (g_tmp_tmp - 1.0) *
                                c_fh3_tmp / d_Smod_tmp_tmp) +
                           i_fh3_tmp * c_fh3_tmp / e_fh3_tmp) /
                          e_del_tmp +
                      d_Smod_tmp *
                          ((fh3_tmp_tmp + b_fh3_tmp_tmp * g_tmp_tmp) -
                           f_Smod_tmp_tmp * c_fh3_tmp / n_fh3_tmp) /
                          e_del_tmp) -
                     j_Smod_tmp * c_fh3_tmp * d_Smod_tmp / b_Smod_tmp) -
                    g_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
                   d_del_tmp * o_fh3_tmp * b_a_tmp * c_fh3_tmp / e_fh3_tmp) +
                  2.0 * d_del_tmp * g_tmp * b_a_tmp * e_Smod_tmp * c_fh3_tmp /
                      d_Smod_tmp_tmp) -
                 theta) -
                h_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
               del) -
              i_fh3_tmp * e_Smod_tmp * c_fh3_tmp / p_fh3_tmp) -
             Dq) +
            f_fh3_tmp) +
           fh3_tmp) -
          vJ_tmp;
      e_fh3_tmp = b_del_tmp * g_tmp * c_fh3_tmp;
      fh3[13] = (d * b_fh3_tmp_tmp_tmp * (g_tmp_tmp - 1.0) / e_del_tmp -
                 e_fh3_tmp / j_fh3_tmp) -
                f_Smod_tmp * c_fh3_tmp / b_Smod_tmp;
      vJ_tmp = c_Smod_tmp_tmp * c_fh3_tmp;
      g_fh3_tmp = (d_fh3_tmp_tmp * f_del_tmp / 2.0 -
                   3.0 * c_fh3_tmp * f_del_tmp / 2.0) +
                  vJ_tmp;
      h_fh3_tmp = 2.0 * c_fh3_tmp_tmp;
      fh3[2] = ((fh3_tmp_tmp_tmp * b_a_tmp * i_Smod_tmp / k_Smod_tmp -
                 l_Smod_tmp * g_fh3_tmp / k_Smod_tmp) +
                del_tmp * b_fh3_tmp * i_Smod_tmp / k_Smod_tmp) -
               5.0 * del_tmp * b_a_tmp * i_Smod_tmp * c_fh3_tmp / h_fh3_tmp;
      fh3[8] = ((b_fh3_tmp_tmp_tmp * b_a_tmp * i_Smod_tmp / k_Smod_tmp -
                 m_Smod_tmp * g_fh3_tmp / k_Smod_tmp) +
                b_del_tmp * b_fh3_tmp * i_Smod_tmp / k_Smod_tmp) -
               5.0 * b_del_tmp * b_a_tmp * i_Smod_tmp * c_fh3_tmp / h_fh3_tmp;
      fh3[14] = d_fh3_tmp_tmp / (2.0 * e_del_tmp) - vJ_tmp / j_fh3_tmp;
      b_fh3_tmp = b_g_tmp_tmp *
                  (c_fh3_tmp / e_fh3_tmp_tmp - d_fh3_tmp_tmp / e_fh3_tmp_tmp) /
                  n_Smod_tmp;
      g_fh3_tmp = del_tmp * b_fh3_tmp_tmp_tmp * h_Smod_tmp_tmp / n_Smod_tmp;
      h_fh3_tmp = fh3_tmp_tmp_tmp * b_del_tmp * h_Smod_tmp_tmp / n_Smod_tmp;
      i_fh3_tmp =
          3.0 * del_tmp * b_del_tmp * h_Smod_tmp_tmp * c_fh3_tmp / p_fh3_tmp;
      fh3[3] = ((i_fh3_tmp - g_fh3_tmp) - h_fh3_tmp) - b_fh3_tmp;
      fh3[9] = 3.0 * p_Smod_tmp * c_fh3_tmp / p_fh3_tmp -
               (((d_del_tmp * c_fh3_tmp / e_fh3_tmp_tmp +
                  b_fh3_tmp_tmp * f_del_tmp) +
                 2.0 * d * del_tmp * fh3_tmp_tmp_tmp * g_tmp) +
                e_Smod_tmp_tmp * c_fh3_tmp / e_fh3_tmp_tmp) /
                   n_Smod_tmp;
      fh3[15] = 0.0;
      fh3[4] =
          (((c_del_tmp * c_fh3_tmp / e_fh3_tmp_tmp + fh3_tmp_tmp * f_del_tmp) +
            2.0 * d * b_del_tmp * b_fh3_tmp_tmp_tmp * g_tmp) +
           g_Smod_tmp_tmp * c_fh3_tmp / e_fh3_tmp_tmp) /
              n_Smod_tmp -
          3.0 * q_Smod_tmp * c_fh3_tmp / p_fh3_tmp;
      fh3[10] = ((b_fh3_tmp + g_fh3_tmp) + h_fh3_tmp) - i_fh3_tmp;
      fh3[16] = 0.0;
      b_fh3_tmp = 2.0 * d * c_Smod_tmp;
      fh3[5] = (o_Smod_tmp * c_fh3_tmp / b_Smod_tmp -
                b_fh3_tmp_tmp_tmp * (g_tmp_tmp - 1.0) / e_del_tmp) +
               e_fh3_tmp / b_fh3_tmp;
      fh3[11] = (fh3_tmp_tmp_tmp * (g_tmp_tmp - 1.0) / e_del_tmp -
                 r_Smod_tmp * c_fh3_tmp / b_Smod_tmp) -
                d_fh3_tmp / b_fh3_tmp;
      fh3[17] = 0.0;
    }
    /*  convert transform to rotation and translation */
    /*  Adjoint calculator */
    /*  This function calculates the adjoint given simply the transform between
     */
    /*  joints */
    /*  get skew symmetric matrix of translation */
    for (i = 0; i < 3; i++) {
      skw[3 * i] = g[i];
      skw[3 * i + 1] = g[i + 4];
      skw[3 * i + 2] = g[i + 8];
    }
    for (i = 0; i < 9; i++) {
      b_skw[i] = -skw[i];
    }
    Smod[3] = -g[14];
    Smod[6] = g[13];
    Smod[1] = g[14];
    Smod[7] = -g[12];
    Smod[2] = -g[13];
    Smod[5] = g[12];
    for (i = 0; i < 3; i++) {
      fh3_tmp = b_skw[i];
      Smod_tmp_tmp = b_skw[i + 3];
      theta = b_skw[i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        c_skw[i + 3 * i1] =
            (fh3_tmp * Smod[3 * i1] + Smod_tmp_tmp * Smod[3 * i1 + 1]) +
            theta * Smod[3 * i1 + 2];
        c_I[i1 + 6 * i] = skw[i1 + 3 * i];
      }
    }
    for (i = 0; i < 3; i++) {
      I_tmp = 6 * (i + 3);
      c_I[I_tmp] = c_skw[3 * i];
      c_I[6 * i + 3] = 0.0;
      c_I[I_tmp + 3] = skw[3 * i];
      ibtile = 3 * i + 1;
      c_I[I_tmp + 1] = c_skw[ibtile];
      c_I[6 * i + 4] = 0.0;
      c_I[I_tmp + 4] = skw[ibtile];
      ibtile = 3 * i + 2;
      c_I[I_tmp + 2] = c_skw[ibtile];
      c_I[6 * i + 5] = 0.0;
      c_I[I_tmp + 5] = skw[ibtile];
    }
    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        fh3_tmp = 0.0;
        for (i2 = 0; i2 < 6; i2++) {
          fh3_tmp += c_I[i + 6 * i2] * (double)Xtree[i2 + 6 * i1];
        }
        Xup[(i + 6 * i1) + 36 * (b_i + 1)] = fh3_tmp;
      }
    }
    for (i = 0; i < 36; i++) {
      Xtree[i] = 0;
    }
    vJ_tmp = dq[qi_data[0] - 1];
    del = dq[qi_data[1] - 1];
    c_a = dq[qi_data[2] - 1];
    for (ibtile = 0; ibtile < 6; ibtile++) {
      Xtree[ibtile + 6 * ibtile] = 1;
      i = ibtile + 18 * b_i;
      fh3_tmp =
          (b_Smod[i] * vJ_tmp + b_Smod[i + 6] * del) + b_Smod[i + 12] * c_a;
      vJ[ibtile] = fh3_tmp;
      Smod_tmp_tmp = 0.0;
      for (i = 0; i < 6; i++) {
        Smod_tmp_tmp += Xup[(ibtile + 6 * i) + 36 * (b_i + 1)] * v[i + 6 * b_i];
      }
      b_Xup[ibtile] = Smod_tmp_tmp + fh3_tmp;
    }
    for (i = 0; i < 6; i++) {
      v[i + 6 * (b_i + 1)] = b_Xup[i];
    }
    memset(&c_I[0], 0, 36U * sizeof(double));
    skw[0] = 0.0;
    i = 6 * (b_i + 1);
    fh3_tmp = v[i + 5];
    skw[3] = -fh3_tmp;
    Smod_tmp_tmp = v[i + 4];
    skw[6] = Smod_tmp_tmp;
    skw[1] = fh3_tmp;
    skw[4] = 0.0;
    fh3_tmp = v[i + 3];
    skw[7] = -fh3_tmp;
    skw[2] = -Smod_tmp_tmp;
    skw[5] = fh3_tmp;
    skw[8] = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      theta = skw[3 * i1];
      c_I[6 * i1] = theta;
      I_tmp = 6 * (i1 + 3);
      c_I[I_tmp + 3] = theta;
      theta = skw[3 * i1 + 1];
      c_I[6 * i1 + 1] = theta;
      c_I[I_tmp + 4] = theta;
      theta = skw[3 * i1 + 2];
      c_I[6 * i1 + 2] = theta;
      c_I[I_tmp + 5] = theta;
    }
    c_I[18] = 0.0;
    fh3_tmp = v[i + 2];
    c_I[24] = -fh3_tmp;
    Smod_tmp_tmp = v[i + 1];
    c_I[30] = Smod_tmp_tmp;
    c_I[19] = fh3_tmp;
    c_I[25] = 0.0;
    fh3_tmp = v[i];
    c_I[31] = -fh3_tmp;
    c_I[20] = -Smod_tmp_tmp;
    c_I[26] = fh3_tmp;
    c_I[32] = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      fh3_tmp = 0.0;
      for (i2 = 0; i2 < 6; i2++) {
        fh3_tmp += Xup[(i1 + 6 * i2) + 36 * (b_i + 1)] * a[i2 + 6 * b_i];
      }
      b_Xup[i1] = fh3_tmp;
      d_I[i1] = (fh3[i1] * vJ_tmp + fh3[i1 + 6] * del) + fh3[i1 + 12] * c_a;
    }
    for (i1 = 0; i1 < 6; i1++) {
      fh3_tmp = 0.0;
      for (i2 = 0; i2 < 6; i2++) {
        Smod_tmp_tmp = c_I[i1 + 6 * i2];
        fh3_tmp += Smod_tmp_tmp * vJ[i2];
        a_tmp[i2 + 6 * i1] = -Smod_tmp_tmp;
      }
      a[i1 + i] = (b_Xup[i1] + d_I[i1]) + fh3_tmp;
    }
    for (i1 = 0; i1 < 6; i1++) {
      d_I[i1] = 0.0;
      vJ[i1] = 0.0;
      for (i2 = 0; i2 < 6; i2++) {
        fh3_tmp = 0.0;
        for (ibtile = 0; ibtile < 6; ibtile++) {
          fh3_tmp +=
              a_tmp[i1 + 6 * ibtile] * b_I[(ibtile + 6 * i2) + 36 * (b_i + 1)];
        }
        I_tmp = i2 + i;
        d_I[i1] += b_I[(i1 + 6 * i2) + 36 * (b_i + 1)] * a[I_tmp];
        vJ[i1] += fh3_tmp * v[I_tmp];
      }
      f[i1 + i] = d_I[i1] + vJ[i1];
    }
  }
  /*  Composite Rigid Body Algorithm to calculate M */
  /*  composite inertia calculation */
  for (b_i = 0; b_i < 4; b_i++) {
    if (4 - b_i == 1) {
      c_a = 0.0;
      for (i = 0; i < 6; i++) {
        c_a += (double)d_a[i] * f[i];
      }
      C[0] = c_a;
    } else {
      ibtile = (2 - b_i) * 3;
      for (i = 0; i < 3; i++) {
        qi_data[i] = (signed char)((ibtile + i) + 2);
      }
      for (i = 0; i < 3; i++) {
        tmp_data[i] = (signed char)(qi_data[i] - 1);
      }
      for (i = 0; i < 3; i++) {
        fh3_tmp = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          fh3_tmp += 2.0 * b_Smod[(i1 + 6 * i) + 18 * (2 - b_i)] *
                     f[i1 + 6 * (3 - b_i)];
        }
        b_dq[i] = fh3_tmp;
      }
      for (i = 0; i < 3; i++) {
        C[tmp_data[i]] = b_dq[i];
      }
      for (i = 0; i < 6; i++) {
        fh3_tmp = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          fh3_tmp += Xup[(i1 + 6 * i) + 36 * (3 - b_i)] * f[i1 + 6 * (3 - b_i)];
        }
        d_I[i] = f[i + 6 * (2 - b_i)] + fh3_tmp;
      }
      for (i = 0; i < 6; i++) {
        f[i + 6 * (2 - b_i)] = d_I[i];
      }
    }
    if (4 - b_i != 1) {
      i = 36 * (3 - b_i);
      for (i1 = 0; i1 < 6; i1++) {
        for (i2 = 0; i2 < 6; i2++) {
          fh3_tmp = 0.0;
          for (ibtile = 0; ibtile < 6; ibtile++) {
            fh3_tmp += Xup[(ibtile + 6 * i1) + i] * b_I[(ibtile + 6 * i2) + i];
          }
          a_tmp[i1 + 6 * i2] = fh3_tmp;
        }
        for (i2 = 0; i2 < 6; i2++) {
          fh3_tmp = 0.0;
          for (ibtile = 0; ibtile < 6; ibtile++) {
            fh3_tmp += a_tmp[i1 + 6 * ibtile] *
                       Xup[(ibtile + 6 * i2) + 36 * (3 - b_i)];
          }
          I_tmp = i1 + 6 * i2;
          c_I[I_tmp] = b_I[I_tmp + 36 * (2 - b_i)] + fh3_tmp;
        }
      }
      memcpy(&b_I[b_i * -36 + 72], &c_I[0], 36U * sizeof(double));
    }
  }
  memset(&M[0], 0, 100U * sizeof(double));
  /*  fh3 = zeros(6,3); */
  /*  fh1 = zeros(6,1); */
  for (b_i = 0; b_i < 4; b_i++) {
    if (b_i + 1 == 1) {
      fh3_tmp = 0.0;
      for (i = 0; i < 6; i++) {
        Smod_tmp_tmp = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          Smod_tmp_tmp += b_I[i + 6 * i1] * (double)b_a[i1];
        }
        fh3_tmp += (double)iv[i] * Smod_tmp_tmp;
      }
      M[0] = fh3_tmp;
    } else {
      signed char b_tmp_data[9];
      ibtile = (b_i - 1) * 3;
      for (i = 0; i < 3; i++) {
        qi_data[i] = (signed char)((ibtile + i) + 2);
      }
      for (i = 0; i < 6; i++) {
        for (i1 = 0; i1 < 3; i1++) {
          fh3_tmp = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            fh3_tmp += b_I[(i + 6 * i2) + 36 * b_i] *
                       b_Smod[(i2 + 6 * i1) + 18 * (b_i - 1)];
          }
          fh3[i + 6 * i1] = fh3_tmp;
        }
      }
      for (i = 0; i < 3; i++) {
        i1 = qi_data[i];
        tmp_data[i] = (signed char)(i1 - 1);
        b_tmp_data[i] = (signed char)(i1 - 1);
      }
      for (i = 0; i < 3; i++) {
        for (i1 = 0; i1 < 3; i1++) {
          fh3_tmp = 0.0;
          for (i2 = 0; i2 < 6; i2++) {
            fh3_tmp += b_Smod[(i2 + 6 * i) + 18 * (b_i - 1)] * fh3[i2 + 6 * i1];
          }
          Smod[i + 3 * i1] = fh3_tmp;
        }
      }
      for (i = 0; i < 3; i++) {
        for (i1 = 0; i1 < 3; i1++) {
          M[tmp_data[i1] + 10 * b_tmp_data[i]] = Smod[i1 + 3 * i];
        }
      }
      j = b_i - 1;
      while (j + 2 > 1) {
        for (i = 0; i < 6; i++) {
          for (i1 = 0; i1 < 3; i1++) {
            fh3_tmp = 0.0;
            for (i2 = 0; i2 < 6; i2++) {
              fh3_tmp += Xup[(i2 + 6 * i) + 36 * (j + 1)] * fh3[i2 + 6 * i1];
            }
            c_Xup[i + 6 * i1] = fh3_tmp;
          }
        }
        memcpy(&fh3[0], &c_Xup[0], 18U * sizeof(double));
        j--;
        if (j + 2 == 1) {
          for (i = 0; i < 3; i++) {
            fh3_tmp = 0.0;
            for (i1 = 0; i1 < 6; i1++) {
              fh3_tmp += fh3[i1 + 6 * i] * (double)b_a[i1];
            }
            b_dq[i] = fh3_tmp;
          }
          for (i = 0; i < 3; i++) {
            M[b_tmp_data[i]] = b_dq[i];
          }
          /* (S{j}' * fh).'; */
          for (i = 0; i < 3; i++) {
            fh3_tmp = 0.0;
            for (i1 = 0; i1 < 6; i1++) {
              fh3_tmp += (double)iv[i1] * fh3[i1 + 6 * i];
            }
            b_dq[i] = fh3_tmp;
          }
          for (i = 0; i < 3; i++) {
            M[10 * b_tmp_data[i]] = b_dq[i];
          }
        } else {
          int qi_size_idx_1;
          signed char c_tmp_data[6];
          signed char d_tmp_data[6];
          I_tmp = j * 3 + 2;
          ibtile = j * 3 + 4;
          if (ibtile < I_tmp) {
            qi_size_idx_1 = 0;
          } else {
            ibtile = (signed char)ibtile - (signed char)I_tmp;
            qi_size_idx_1 = ibtile + 1;
            for (i = 0; i <= ibtile; i++) {
              qi_data[i] = (signed char)((signed char)I_tmp + (signed char)i);
            }
          }
          for (i = 0; i < qi_size_idx_1; i++) {
            i1 = qi_data[i];
            c_tmp_data[i] = (signed char)i1;
            d_tmp_data[i] = (signed char)(i1 - 1);
          }
          for (i = 0; i < 3; i++) {
            for (i1 = 0; i1 < 3; i1++) {
              fh3_tmp = 0.0;
              for (i2 = 0; i2 < 6; i2++) {
                fh3_tmp += fh3[i2 + 6 * i] * b_Smod[(i2 + 6 * i1) + 18 * j];
              }
              Smod[i + 3 * i1] = fh3_tmp;
            }
          }
          /* (S{j}' * fh).'; */
          for (i = 0; i < qi_size_idx_1; i++) {
            for (i1 = 0; i1 < 3; i1++) {
              M[b_tmp_data[i1] + 10 * d_tmp_data[i]] = Smod[i1 + 3 * i];
            }
            d_tmp_data[i] = (signed char)(c_tmp_data[i] - 1);
          }
          for (i = 0; i < 3; i++) {
            for (i1 = 0; i1 < 3; i1++) {
              fh3_tmp = 0.0;
              for (i2 = 0; i2 < 6; i2++) {
                fh3_tmp += b_Smod[(i2 + 6 * i) + 18 * j] * fh3[i2 + 6 * i1];
              }
              Smod[i + 3 * i1] = fh3_tmp;
            }
          }
          for (i = 0; i < 3; i++) {
            for (i1 = 0; i1 < qi_size_idx_1; i1++) {
              M[d_tmp_data[i1] + 10 * b_tmp_data[i]] =
                  Smod[i1 + qi_size_idx_1 * i];
            }
          }
        }
      }
    }
  }
  /*   */
  /*  for i = 1:N */
  /*      if mod(i,2) == 0 */
  /*          qi = 2 + (2*i-4):4 + (2*i-4); */
  /*          fh3 = IC(:,:,i) * Smod(:,:,i/2); */
  /*          H(qi,qi) = Smod(:,:,i/2)' * fh3; */
  /*          j = i; */
  /*          while j < N+1 && j > 1 */
  /*              fh3 = Xup(:,:,j)' * fh3; */
  /*              j = j - 1; */
  /*              if mod(j,2) == 0 */
  /*                  qj = 2 + (2*j-4):4 + (2*j-4); */
  /*                  H(qi,qj) = fh3' * Smod(:,:,j/2); %(S{j}' * fh).'; */
  /*                  H(qj,qi) = Smod(:,:,j/2)' * fh3; */
  /*              elseif j == 1 */
  /*                  qj = 1 + (2*j-2); */
  /*                  H(qi,qj) = fh3' * Smotz; %(S{j}' * fh).'; */
  /*                  H(qj,qi) = Smotz' * fh3; */
  /*              else */
  /*                  qj = 1 + (2*j-2); */
  /*                  H(qi,qj) = fh3' * Smoty; %(S{j}' * fh).'; */
  /*                  H(qj,qi) = Smoty' * fh3; */
  /*              end */
  /*          end */
  /*      else */
  /*          qi = 1 + (2*i-2); */
  /*          if i == 1 */
  /*              fh1 = IC(:,:,i) * Smotz; */
  /*              H(qi,qi) = Smotz' * fh1; */
  /*          else */
  /*              fh1 = IC(:,:,i) * Smoty; */
  /*              H(qi,qi) = Smoty' * fh1; */
  /*          end */
  /*   */
  /*          j = i; */
  /*          while j < N+1 && j > 1 */
  /*              fh1 = Xup(:,:,j)' * fh1; */
  /*              j = j - 1; */
  /*              if mod(j,2) == 0 */
  /*                  qj = 2 + (2*j-4):4 + (2*j-4); */
  /*                  H(qi,qj) = fh1' * Smod(:,:,j/2); %(S{j}' * fh).'; */
  /*                  H(qj,qi) = Smod(:,:,j/2)' * fh1; */
  /*              elseif j == 1 */
  /*                  qj = 1 + (2*j-2); */
  /*                  H(qi,qj) = fh1' * Smotz; %(S{j}' * fh).'; */
  /*                  H(qj,qi) = Smotz' * fh1; */
  /*              else */
  /*                  qj = 1 + (2*j-2); */
  /*                  H(qi,qj) = fh1' * Smoty; %(S{j}' * fh).'; */
  /*                  H(qj,qi) = Smoty' * fh1; */
  /*              end */
  /*          end */
  /*      end */
  /*   */
  /*  end */
  b_dq[0] = kb;
  b_dq[1] = kb;
  b_dq[2] = ks;
  repmat(b_dq, skw);
  b_bm[0] = 0.0;
  memcpy(&b_bm[1], &skw[0], 9U * sizeof(double));
  diag(b_bm, K);
  b_dq[0] = bb;
  b_dq[1] = bb;
  b_dq[2] = bs;
  repmat(b_dq, skw);
  b_bm[0] = bm;
  memcpy(&b_bm[1], &skw[0], 9U * sizeof(double));
  diag(b_bm, D);
  /*  Fci = cell(1,4); */
  /*  Ai = cell(1,4); */
  for (b_i = 0; b_i < 3; b_i++) {
    fh3_tmp = qp[3 * b_i];
    Smod_tmp_tmp = qp[3 * b_i + 1];
    del_tmp = fh3_tmp * fh3_tmp;
    b_del_tmp = Smod_tmp_tmp * Smod_tmp_tmp;
    del = sqrt(del_tmp + b_del_tmp);
    theta = del / d;
    c_a = sin(del);
    Dq = del - c_a;
    f_fh3_tmp = qp[3 * b_i + 2] + L0;
    if ((fh3_tmp < 1.0E-5) && (Smod_tmp_tmp < 1.0E-5)) {
      cq[b_i] = f_fh3_tmp / 3.0;
      for (i = 0; i < 9; i++) {
        skw[i] = iv1[i];
      }
    } else {
      cq[b_i] = 2.0 * (f_fh3_tmp / theta - d) * sin(theta / 6.0);
      vJ_tmp = rt_powd_snf(del, 3.0);
      skw[0] = fh3_tmp * Smod_tmp_tmp * Dq / vJ_tmp;
      skw[3] = (-del_tmp * del - b_del_tmp * c_a) / vJ_tmp;
      skw[6] = fh3_tmp * Dq * f_fh3_tmp / vJ_tmp;
      skw[1] = (b_del_tmp * del + del_tmp * c_a) / vJ_tmp;
      skw[4] = -fh3_tmp * Smod_tmp_tmp * Dq / vJ_tmp;
      skw[7] = Smod_tmp_tmp * Dq * f_fh3_tmp / vJ_tmp;
      skw[2] = 0.0;
      skw[5] = 0.0;
      skw[8] = c_a / del;
    }
    Smod[0] = d * 0.49999999999999994;
    Smod[3] = d * 0.49999999999999994;
    Smod[6] = -d;
    Smod[1] = -d * 0.86602540378443871;
    Smod[4] = d * 0.86602540378443871;
    Smod[7] = 0.0;
    for (i = 0; i < 3; i++) {
      fh3_tmp = skw[i];
      Smod_tmp_tmp = skw[i + 3];
      theta = skw[i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        b_skw[i + 3 * i1] = (fh3_tmp * (double)c_b[3 * i1] +
                             Smod_tmp_tmp * (double)c_b[3 * i1 + 1]) +
                            theta * (double)c_b[3 * i1 + 2];
      }
      Smod[3 * i + 2] = 1.0;
    }
    for (i = 0; i < 3; i++) {
      fh3_tmp = b_skw[i];
      Smod_tmp_tmp = b_skw[i + 3];
      theta = b_skw[i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        Ai[(i + 3 * i1) + 9 * b_i] =
            (fh3_tmp * Smod[3 * i1] + Smod_tmp_tmp * Smod[3 * i1 + 1]) +
            theta * Smod[3 * i1 + 2];
      }
    }
  }
  blkdiag(&Ai[0], &Ai[9], &Ai[18], A);
  for (I_tmp = 0; I_tmp < 3; I_tmp++) {
    ibtile = I_tmp * 3;
    skw[ibtile] = conv_pcc;
    skw[ibtile + 1] = conv_pcc;
    skw[ibtile + 2] = conv_pcc;
  }
  b_v[0] = conv_motor;
  memcpy(&b_v[1], &skw[0], 9U * sizeof(double));
  memset(&conversion[0], 0, 100U * sizeof(double));
  for (j = 0; j < 10; j++) {
    conversion[j + 10 * j] = b_v[j];
    fh3_tmp = 0.0;
    b_bm[j] = qd[j] - q[j];
    Smod_tmp_tmp = 0.0;
    theta = 0.0;
    for (i = 0; i < 10; i++) {
      i1 = j + 10 * i;
      fh3_tmp += M[i1] * ddqd[i];
      theta += K[i1] * qd[i];
      Smod_tmp_tmp += D[i1] * dqd[i];
    }
    b_C[j] = ((C[j] + fh3_tmp) + theta) + Smod_tmp_tmp;
  }
  for (i = 0; i < 10; i++) {
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      fh3_tmp += Kp[i + 10 * i1] * b_bm[i1];
    }
    b_K[i] = fh3_tmp;
  }
  for (i = 0; i < 10; i++) {
    b_bm[i] = dqd[i] - dq[i];
  }
  for (i = 0; i < 10; i++) {
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      fh3_tmp += KD[i + 10 * i1] * b_bm[i1];
    }
    b_v[i] = (b_C[i] + b_K[i]) + fh3_tmp;
  }
  mldivide(A, b_v);
  for (i = 0; i < 10; i++) {
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      fh3_tmp += conversion[i + 10 * i1] * b_v[i1];
    }
    tau[i] = fh3_tmp;
  }
  memcpy(&Kpr[0], &Kp[0], 100U * sizeof(double));
  Kpr[0] = 0.0;
  memcpy(&KDr[0], &KD[0], 100U * sizeof(double));
  KDr[0] = 0.0;
  memcpy(&b_C[0], &q[0], 10U * sizeof(double));
  J_r(b_C, L0, d, J, x);
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 10; i1++) {
      b_tmp_tmp[i1 + 10 * i] = J[i + 3 * i1];
    }
  }
  memcpy(&b_tmp[0], &b_tmp_tmp[0], 30U * sizeof(double));
  b_mldivide(M, b_tmp);
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      fh3_tmp = 0.0;
      for (i2 = 0; i2 < 10; i2++) {
        fh3_tmp += J[i + 3 * i2] * b_tmp[i2 + 10 * i1];
      }
      Smod[i + 3 * i1] = fh3_tmp;
    }
  }
  mrdiv(b_A, Smod, skw);
  for (i = 0; i < 10; i++) {
    fh3_tmp = b_tmp[i];
    Smod_tmp_tmp = b_tmp[i + 10];
    theta = b_tmp[i + 20];
    for (i1 = 0; i1 < 3; i1++) {
      b_b_tmp[i + 10 * i1] =
          (fh3_tmp * skw[3 * i1] + Smod_tmp_tmp * skw[3 * i1 + 1]) +
          theta * skw[3 * i1 + 2];
    }
  }
  mtimes(J, b_b_tmp, y_tmp);
  memset(&e_I[0], 0, 100U * sizeof(signed char));
  for (ibtile = 0; ibtile < 10; ibtile++) {
    e_I[ibtile + 10 * ibtile] = 1;
    fh3_tmp = 0.0;
    Smod_tmp_tmp = 0.0;
    for (i = 0; i < 10; i++) {
      i1 = ibtile + 10 * i;
      fh3_tmp += K[i1] * q[i];
      Smod_tmp_tmp += D[i1] * dq[i];
    }
    b_K[ibtile] = fh3_tmp + Smod_tmp_tmp;
    fh3_tmp = b_tmp_tmp[ibtile];
    Smod_tmp_tmp = b_tmp_tmp[ibtile + 10];
    theta = b_tmp_tmp[ibtile + 20];
    for (i = 0; i < 3; i++) {
      b_tmp[ibtile + 10 * i] =
          (fh3_tmp * skw[3 * i] + Smod_tmp_tmp * skw[3 * i + 1]) +
          theta * skw[3 * i + 2];
    }
  }
  b_dq[0] = Kpx * (xd[0] - x[0]) + KDx * (dxd[0] - dxr[0]);
  b_dq[1] = Kpx * (xd[1] - x[1]) + KDx * (dxd[1] - dxr[1]);
  b_dq[2] = Kpx * (xd[2] - x[2]) + KDx * (dxd[2] - dxr[2]);
  for (i = 0; i < 10; i++) {
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      fh3_tmp += y_tmp[i + 10 * i1] * b_K[i1];
    }
    b_v[i] =
        (C[i] + fh3_tmp) + ((b_tmp[i] * b_dq[0] + b_tmp[i + 10] * b_dq[1]) +
                            b_tmp[i + 20] * b_dq[2]);
  }
  mldivide(A, b_v);
  for (i = 0; i < 100; i++) {
    Kpr[i] = -Kpr[i];
  }
  for (i = 0; i < 10; i++) {
    fh3_tmp = 0.0;
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      i2 = i + 10 * i1;
      fh3_tmp += Kpr[i2] * q[i1];
      Smod_tmp_tmp += KDr[i2] * dq[i1];
    }
    b_bm[i] = Smod_tmp_tmp;
    b_K[i] = fh3_tmp;
  }
  for (i = 0; i < 100; i++) {
    y_tmp[i] = (double)e_I[i] - y_tmp[i];
  }
  for (i = 0; i < 10; i++) {
    b_K[i] -= b_bm[i];
  }
  for (i = 0; i < 10; i++) {
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      fh3_tmp += y_tmp[i + 10 * i1] * b_K[i1];
    }
    b_bm[i] = b_v[i] + fh3_tmp;
  }
  for (i = 0; i < 10; i++) {
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      fh3_tmp += conversion[i + 10 * i1] * b_bm[i1];
    }
    tau_r[i] = fh3_tmp;
  }
}

/*
 * File trailer for helix_controller.c
 *
 * [EOF]
 */
