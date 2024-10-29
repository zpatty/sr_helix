/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: helix_controller.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 23:07:32
 */

/* Include Files */
#include "helix_controller.h"
#include "J_r.h"
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
  static const double b_b[6] = {-0.0, -0.0, -4.905, -0.0, -0.0, -0.0};
  static const double v[6] = {0.738,
                              0.738,
                              0.738,
                              0.0020730419999999998,
                              0.0020730419999999998,
                              0.0041460839999999995};
  static const signed char b[36] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                    0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  static const signed char c_b[9] = {-1, 0, 0, 0, -1, 0, 0, 0, -1};
  static const signed char d_a[9] = {0, 1, 0, -1, 0, 0, 0, 0, 1};
  static const signed char iv1[9] = {0, 1, 0, -1, 0, 0, 0, 0, 1};
  static const signed char b_a[6] = {0, 0, 0, 0, 0, 1};
  static const signed char c_a[6] = {0, 0, 0, 0, 0, 2};
  static const signed char iv[6] = {0, 0, 0, 0, 0, 1};
  double XJ[144];
  double Xup[144];
  double b_I[144];
  double D[100];
  double K[100];
  double KDr[100];
  double Kpr[100];
  double conversion[100];
  double y_tmp[100];
  double b_Smod[54];
  double I_tmp[36];
  double b_I_tmp[36];
  double J[30];
  double b_b_tmp[30];
  double b_tmp[30];
  double b_tmp_tmp[30];
  double Ai[27];
  double a[24];
  double b_v[24];
  double f[24];
  double c_Xup[18];
  double fh3[18];
  double g[16];
  double b_C[10];
  double b_K[10];
  double b_qd[10];
  double c_C[10];
  double Al[9];
  double Smod[9];
  double b_skw[9];
  double qp[9];
  double skw[9];
  double c_I[6];
  double vJ[6];
  double b_dq[3];
  double Dq;
  double Smod_tmp_tmp;
  double b_a_tmp;
  double b_del_tmp;
  double b_g_tmp;
  double del;
  double del_tmp;
  double fh3_tmp;
  double g_tmp;
  double theta;
  double vJ_tmp_tmp;
  int XJ_tmp;
  int b_XJ_tmp;
  int b_i;
  int c_XJ_tmp;
  int i;
  int i1;
  int j;
  signed char d_I[100];
  signed char Xtree[36];
  signed char qi_data[9];
  signed char tmp_data[9];
  for (i = 0; i < 3; i++) {
    b_i = 3 * i;
    qp[3 * i] = q[b_i + 1];
    qp[3 * i + 1] = q[b_i + 2];
    qp[3 * i + 2] = q[b_i + 3];
  }
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
  /*  for i=1:N-1 */
  /*      q(i*3+1) = -q(i*3+1); */
  /*      qd(i*3+1) = -qd(i*3+1); */
  /*  end */
  /*  qi = {1, 2:4, 5, 6:8, 9, 10:12, 13, 14:16}; */
  memset(&b_I[0], 0, 144U * sizeof(double));
  memset(&b_I[0], 0, 36U * sizeof(double));
  for (j = 0; j < 6; j++) {
    b_I[j + 6 * j] = v[j];
  }
  vJ[0] = m;
  vJ[1] = m;
  vJ[2] = m;
  vJ_tmp_tmp = r * r;
  theta = 0.25 * m * vJ_tmp_tmp;
  vJ[3] = theta;
  vJ[4] = theta;
  vJ[5] = 0.5 * m * vJ_tmp_tmp;
  memset(&I_tmp[0], 0, 36U * sizeof(double));
  for (j = 0; j < 6; j++) {
    I_tmp[j + 6 * j] = vJ[j];
  }
  memcpy(&b_I[36], &I_tmp[0], 36U * sizeof(double));
  for (i = 0; i < 2; i++) {
    memcpy(&b_I[i * 36 + 72], &I_tmp[0], 36U * sizeof(double));
  }
  memset(&XJ[0], 0, 144U * sizeof(double));
  memset(&Xup[0], 0, 144U * sizeof(double));
  memset(&b_v[0], 0, 24U * sizeof(double));
  memset(&a[0], 0, 24U * sizeof(double));
  memset(&f[0], 0, 24U * sizeof(double));
  memset(&C[0], 0, 10U * sizeof(double));
  for (b_i = 0; b_i < 36; b_i++) {
    Xtree[b_i] = b[b_i];
  }
  /* [diag(ones(3,1)) [0;0;0]; 0 0 0 1]; */
  g_tmp = sin(q[0]);
  b_g_tmp = cos(q[0]);
  g[0] = b_g_tmp;
  g[4] = -g_tmp;
  g[8] = 0.0;
  g[12] = 0.0;
  g[1] = g_tmp;
  g[5] = b_g_tmp;
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
  for (b_i = 0; b_i < 3; b_i++) {
    skw[3 * b_i] = g[b_i];
    skw[3 * b_i + 1] = g[b_i + 4];
    skw[3 * b_i + 2] = g[b_i + 8];
  }
  for (b_i = 0; b_i < 9; b_i++) {
    b_skw[b_i] = -skw[b_i];
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
  for (b_i = 0; b_i < 3; b_i++) {
    Smod_tmp_tmp = b_skw[b_i];
    fh3_tmp = b_skw[b_i + 3];
    del = b_skw[b_i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Al[b_i + 3 * i1] =
          (Smod_tmp_tmp * Smod[3 * i1] + fh3_tmp * Smod[3 * i1 + 1]) +
          del * Smod[3 * i1 + 2];
      XJ[i1 + 6 * b_i] = skw[i1 + 3 * b_i];
    }
  }
  for (b_i = 0; b_i < 3; b_i++) {
    XJ_tmp = 6 * (b_i + 3);
    XJ[XJ_tmp] = Al[3 * b_i];
    XJ[6 * b_i + 3] = 0.0;
    XJ[XJ_tmp + 3] = skw[3 * b_i];
    b_XJ_tmp = 3 * b_i + 1;
    XJ[XJ_tmp + 1] = Al[b_XJ_tmp];
    XJ[6 * b_i + 4] = 0.0;
    XJ[XJ_tmp + 4] = skw[b_XJ_tmp];
    b_XJ_tmp = 3 * b_i + 2;
    XJ[XJ_tmp + 2] = Al[b_XJ_tmp];
    XJ[6 * b_i + 5] = 0.0;
    XJ[XJ_tmp + 5] = skw[b_XJ_tmp];
  }
  for (i = 0; i < 6; i++) {
    for (b_i = 0; b_i < 6; b_i++) {
      Smod_tmp_tmp = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        Smod_tmp_tmp += XJ[i + 6 * i1] * (double)b[i1 + 6 * b_i];
      }
      Xup[i + 6 * b_i] = Smod_tmp_tmp;
    }
    Smod_tmp_tmp = (double)b_a[i] * dq[0];
    vJ[i] = Smod_tmp_tmp;
    b_v[i] = Smod_tmp_tmp;
  }
  memset(&I_tmp[0], 0, 36U * sizeof(double));
  skw[0] = 0.0;
  skw[3] = -b_v[5];
  skw[6] = b_v[4];
  skw[1] = b_v[5];
  skw[4] = 0.0;
  skw[7] = -b_v[3];
  skw[2] = -b_v[4];
  skw[5] = b_v[3];
  skw[8] = 0.0;
  for (b_i = 0; b_i < 3; b_i++) {
    theta = skw[3 * b_i];
    I_tmp[6 * b_i] = theta;
    XJ_tmp = 6 * (b_i + 3);
    I_tmp[XJ_tmp + 3] = theta;
    theta = skw[3 * b_i + 1];
    I_tmp[6 * b_i + 1] = theta;
    I_tmp[XJ_tmp + 4] = theta;
    theta = skw[3 * b_i + 2];
    I_tmp[6 * b_i + 2] = theta;
    I_tmp[XJ_tmp + 5] = theta;
  }
  I_tmp[18] = 0.0;
  I_tmp[24] = -b_v[2];
  I_tmp[30] = b_v[1];
  I_tmp[19] = b_v[2];
  I_tmp[25] = 0.0;
  I_tmp[31] = -b_v[0];
  I_tmp[20] = -b_v[1];
  I_tmp[26] = b_v[0];
  I_tmp[32] = 0.0;
  for (b_i = 0; b_i < 6; b_i++) {
    Smod_tmp_tmp = 0.0;
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      b_XJ_tmp = b_i + 6 * i1;
      Smod_tmp_tmp += Xup[b_XJ_tmp] * b_b[i1];
      del = I_tmp[b_XJ_tmp];
      fh3_tmp += del * vJ[i1];
      b_I_tmp[i1 + 6 * b_i] = -del;
    }
    a[b_i] = Smod_tmp_tmp + fh3_tmp;
  }
  for (b_i = 0; b_i < 6; b_i++) {
    c_I[b_i] = 0.0;
    vJ[b_i] = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      Smod_tmp_tmp = 0.0;
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        Smod_tmp_tmp += b_I_tmp[b_i + 6 * b_XJ_tmp] * b_I[b_XJ_tmp + 6 * i1];
      }
      c_I[b_i] += b_I[b_i + 6 * i1] * a[i1];
      vJ[b_i] += Smod_tmp_tmp * b_v[i1];
    }
    f[b_i] = c_I[b_i] + vJ[b_i];
  }
  /*  Recursive Newton Euler to Calculate C+G */
  Smod[0] = 0.0;
  Smod[4] = 0.0;
  Smod[8] = 0.0;
  for (i = 0; i < 3; i++) {
    double b_Xup[6];
    double c_del_tmp;
    double d_del_tmp;
    double e_del_tmp;
    double f_del_tmp;
    double n_fh3_tmp;
    double o_fh3_tmp;
    double p_fh3_tmp;
    b_XJ_tmp = i * 3 + 2;
    XJ_tmp = i * 3 + 4;
    if (XJ_tmp < b_XJ_tmp) {
    } else {
      XJ_tmp -= b_XJ_tmp;
      for (b_i = 0; b_i <= XJ_tmp; b_i++) {
        qi_data[b_i] = (signed char)(b_XJ_tmp + b_i);
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
      for (b_i = 0; b_i < 3; b_i++) {
        XJ_tmp = b_i << 2;
        g[XJ_tmp] = skw[3 * b_i];
        g[XJ_tmp + 1] = skw[3 * b_i + 1];
        g[XJ_tmp + 2] = skw[3 * b_i + 2];
      }
      double Smod_tmp;
      double b_fh3_tmp;
      double c_fh3_tmp;
      double d_fh3_tmp;
      g[12] = 0.0;
      g[13] = 0.0;
      Dq = q[qi_data[2] - 1];
      g_tmp = L0 + Dq;
      g[14] = g_tmp;
      g[3] = 0.0;
      g[7] = 0.0;
      g[11] = 0.0;
      g[15] = 1.0;
      Smod_tmp = g_tmp / (2.0 * r);
      b_Smod[18 * i] = Smod_tmp;
      b_Smod[18 * i + 6] = 0.0;
      b_Smod[18 * i + 12] = 0.0;
      b_Smod[18 * i + 1] = 0.0;
      b_Smod[18 * i + 7] = Smod_tmp;
      b_Smod[18 * i + 13] = 0.0;
      b_Smod[18 * i + 3] = 0.0;
      b_Smod[18 * i + 9] = -1.0 / r;
      b_Smod[18 * i + 15] = 0.0;
      b_Smod[18 * i + 4] = 1.0 / r;
      b_Smod[18 * i + 10] = 0.0;
      b_Smod[18 * i + 16] = 0.0;
      b_fh3_tmp = dq[qi_data[2] - 1] / (2.0 * r);
      fh3[0] = b_fh3_tmp;
      fh3[6] = 0.0;
      c_fh3_tmp = dq[qi_data[0] - 1];
      fh3[12] = -c_fh3_tmp / (2.0 * r);
      fh3[1] = 0.0;
      fh3[7] = b_fh3_tmp;
      b_fh3_tmp = dq[qi_data[1] - 1];
      fh3[13] = -b_fh3_tmp / (2.0 * r);
      d_fh3_tmp = 6.0 * vJ_tmp_tmp;
      fh3[2] = c_fh3_tmp * g_tmp / d_fh3_tmp;
      fh3[8] = (L0 * b_fh3_tmp + b_fh3_tmp * Dq) / d_fh3_tmp;
      fh3[14] = 0.0;
      b_Smod[18 * i + 2] = 0.0;
      b_Smod[18 * i + 5] = 0.0;
      fh3[3] = 0.0;
      fh3[4] = 0.0;
      b_Smod[18 * i + 8] = 0.0;
      b_Smod[18 * i + 11] = 0.0;
      fh3[9] = 0.0;
      fh3[10] = 0.0;
      b_Smod[18 * i + 14] = 1.0;
      b_Smod[18 * i + 17] = 0.0;
      fh3[15] = 0.0;
      fh3[16] = 0.0;
      d_fh3_tmp = 2.0 * vJ_tmp_tmp;
      fh3[5] = b_fh3_tmp / d_fh3_tmp;
      fh3[11] = -c_fh3_tmp / d_fh3_tmp;
      fh3[17] = 0.0;
    } else {
      double Smod_tmp;
      double Smod_tmp_tmp_tmp;
      double a_tmp;
      double b_Smod_tmp;
      double b_Smod_tmp_tmp;
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
      double f_fh3_tmp;
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
      double o_Smod_tmp;
      double p_Smod_tmp;
      double q_Smod_tmp;
      double r_Smod_tmp;
      theta = f_del_tmp * f_del_tmp;
      a_tmp = L0 + q[qi_data[2] - 1];
      b_a_tmp = r * a_tmp;
      del = b_a_tmp / theta;
      Dq = f_del_tmp / r;
      g_tmp_tmp = cos(Dq);
      g[0] = c_del_tmp / theta * (g_tmp_tmp - 1.0) + 1.0;
      b_g_tmp_tmp = del_tmp * b_del_tmp;
      g_tmp = b_g_tmp_tmp / theta * (g_tmp_tmp - 1.0);
      g[1] = g_tmp;
      b_g_tmp = sin(Dq);
      g[2] = -del_tmp / f_del_tmp * b_g_tmp;
      g[4] = g_tmp;
      g[5] = d_del_tmp / theta * (g_tmp_tmp - 1.0) + 1.0;
      g[6] = -b_del_tmp / f_del_tmp * b_g_tmp;
      g[8] = del_tmp / f_del_tmp * b_g_tmp;
      g[9] = b_del_tmp / f_del_tmp * b_g_tmp;
      g[10] = g_tmp_tmp;
      g[12] = del * (del_tmp * (1.0 - g_tmp_tmp));
      g[13] = del * (b_del_tmp * (1.0 - g_tmp_tmp));
      g_tmp = f_del_tmp * b_g_tmp;
      g[14] = del * g_tmp;
      g[3] = 0.0;
      g[7] = 0.0;
      g[11] = 0.0;
      g[15] = 1.0;
      b_Smod_tmp_tmp = c_del_tmp * b_g_tmp;
      Smod_tmp = b_Smod_tmp_tmp * a_tmp;
      b_Smod_tmp = e_del_tmp * e_del_tmp;
      c_Smod_tmp = rt_powd_snf(e_del_tmp, 1.5);
      theta = b_a_tmp * (g_tmp_tmp - 1.0);
      d_Smod_tmp = theta / e_del_tmp;
      c_Smod_tmp_tmp = r * b_g_tmp;
      e_Smod_tmp = g_tmp_tmp * f_del_tmp - c_Smod_tmp_tmp;
      d_Smod_tmp_tmp = rt_powd_snf(e_del_tmp, 3.0);
      Smod_tmp_tmp_tmp = c_del_tmp * d_del_tmp;
      del = Smod_tmp_tmp_tmp * a_tmp * (g_tmp_tmp - 1.0);
      Smod_tmp_tmp = (g_tmp - 2.0 * r) + 2.0 * r * g_tmp_tmp;
      f_Smod_tmp = del * Smod_tmp_tmp / d_Smod_tmp_tmp;
      e_Smod_tmp_tmp = c_del_tmp * g_tmp_tmp;
      g_Smod_tmp = d_del_tmp + e_Smod_tmp_tmp;
      h_Smod_tmp = (Smod_tmp / c_Smod_tmp - d_Smod_tmp) +
                   2.0 * r * c_del_tmp * a_tmp * (g_tmp_tmp - 1.0) / b_Smod_tmp;
      b_Smod[18 * i] = (g_Smod_tmp * h_Smod_tmp / e_del_tmp -
                        Smod_tmp * e_Smod_tmp / b_Smod_tmp) +
                       f_Smod_tmp;
      b_Smod[18 * i + 6] = 0.0;
      Smod_tmp = r * del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * i + 12] = Smod_tmp / e_del_tmp;
      b_Smod[18 * i + 1] = 0.0;
      f_Smod_tmp_tmp = d_del_tmp * b_g_tmp;
      i_Smod_tmp = f_Smod_tmp_tmp * a_tmp;
      g_Smod_tmp_tmp = d_del_tmp * g_tmp_tmp;
      j_Smod_tmp = c_del_tmp + g_Smod_tmp_tmp;
      d_Smod_tmp = (i_Smod_tmp / c_Smod_tmp - d_Smod_tmp) +
                   2.0 * r * d_del_tmp * a_tmp * (g_tmp_tmp - 1.0) / b_Smod_tmp;
      b_Smod[18 * i + 7] = (j_Smod_tmp * d_Smod_tmp / e_del_tmp -
                            i_Smod_tmp * e_Smod_tmp / b_Smod_tmp) +
                           f_Smod_tmp;
      f_Smod_tmp = r * b_del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * i + 13] = f_Smod_tmp / e_del_tmp;
      i_Smod_tmp = c_Smod_tmp - c_Smod_tmp_tmp * e_del_tmp;
      k_Smod_tmp = rt_powd_snf(e_del_tmp, 2.5);
      l_Smod_tmp = del_tmp * a_tmp;
      b_Smod[18 * i + 2] = l_Smod_tmp * i_Smod_tmp / k_Smod_tmp;
      m_Smod_tmp = b_del_tmp * a_tmp;
      b_Smod[18 * i + 8] = m_Smod_tmp * i_Smod_tmp / k_Smod_tmp;
      b_Smod[18 * i + 14] = c_Smod_tmp_tmp / f_del_tmp;
      n_Smod_tmp = r * c_Smod_tmp;
      h_Smod_tmp_tmp = f_del_tmp - c_Smod_tmp_tmp;
      o_Smod_tmp = b_g_tmp_tmp * h_Smod_tmp_tmp;
      b_Smod[18 * i + 3] = -o_Smod_tmp / n_Smod_tmp;
      p_Smod_tmp = d_del_tmp * f_del_tmp + r * c_del_tmp * b_g_tmp;
      b_Smod[18 * i + 9] = -p_Smod_tmp / n_Smod_tmp;
      b_Smod[18 * i + 15] = 0.0;
      q_Smod_tmp = c_del_tmp * f_del_tmp + r * d_del_tmp * b_g_tmp;
      b_Smod[18 * i + 4] = q_Smod_tmp / n_Smod_tmp;
      b_Smod[18 * i + 10] = o_Smod_tmp / n_Smod_tmp;
      b_Smod[18 * i + 16] = 0.0;
      o_Smod_tmp = b_del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * i + 5] = -o_Smod_tmp / e_del_tmp;
      r_Smod_tmp = del_tmp * (g_tmp_tmp - 1.0);
      b_Smod[18 * i + 11] = r_Smod_tmp / e_del_tmp;
      b_Smod[18 * i + 17] = 0.0;
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
      d_fh3_tmp = c_fh3_tmp_tmp * b_g_tmp;
      e_fh3_tmp = 2.0 * r * b_Smod_tmp;
      f_fh3_tmp = fh3_tmp_tmp * b_g_tmp * a_tmp;
      g_fh3_tmp = e_Smod_tmp_tmp * a_tmp;
      h_fh3_tmp = 2.0 * c_del_tmp;
      i_fh3_tmp = 2.0 * c_Smod_tmp;
      j_fh3_tmp = b_g_tmp * a_tmp * c_fh3_tmp / i_fh3_tmp -
                  r * b_fh3_tmp * (g_tmp_tmp - 1.0) / e_del_tmp;
      k_fh3_tmp = 2.0 * k_Smod_tmp;
      l_fh3_tmp = 2.0 * r * b_fh3_tmp;
      fh3_tmp = theta * c_fh3_tmp / b_Smod_tmp;
      b_a_tmp = 2.0 * r * f_del_tmp;
      g_tmp = b_g_tmp * b_g_tmp;
      d_fh3_tmp_tmp = g_tmp_tmp * c_fh3_tmp;
      e_fh3_tmp_tmp = 2.0 * f_del_tmp;
      del = del *
            (b_g_tmp * c_fh3_tmp / e_fh3_tmp_tmp - d_fh3_tmp_tmp / (2.0 * r)) /
            d_Smod_tmp_tmp;
      Dq = c_fh3_tmp_tmp * d_del_tmp * (g_tmp_tmp - 1.0) * Smod_tmp_tmp /
           d_Smod_tmp_tmp;
      m_fh3_tmp = 2.0 * r * k_Smod_tmp;
      n_fh3_tmp = 3.0 * c_del_tmp * d_del_tmp * a_tmp * (g_tmp_tmp - 1.0) *
                  c_fh3_tmp * Smod_tmp_tmp / rt_powd_snf(e_del_tmp, 4.0);
      o_fh3_tmp = fh3_tmp_tmp * d_del_tmp * a_tmp * (g_tmp_tmp - 1.0) *
                  Smod_tmp_tmp / d_Smod_tmp_tmp;
      p_fh3_tmp = h_fh3_tmp * b_del_tmp * b_fh3_tmp_tmp_tmp * a_tmp *
                  (g_tmp_tmp - 1.0) * Smod_tmp_tmp / d_Smod_tmp_tmp;
      c_fh3_tmp_tmp = rt_powd_snf(e_del_tmp, 3.5);
      theta = Smod_tmp_tmp_tmp * b_g_tmp * a_tmp * c_fh3_tmp * Smod_tmp_tmp /
              (2.0 * r * c_fh3_tmp_tmp);
      fh3[0] =
          ((((((((((((g_Smod_tmp *
                          ((((((((j_fh3_tmp + d_fh3_tmp / c_Smod_tmp) -
                                 5.0 * c_del_tmp * b_g_tmp * a_tmp * c_fh3_tmp /
                                     k_fh3_tmp) +
                                f_fh3_tmp / c_Smod_tmp) +
                               l_fh3_tmp * c_del_tmp * (g_tmp_tmp - 1.0) /
                                   b_Smod_tmp) +
                              fh3_tmp) +
                             4.0 * r * del_tmp * fh3_tmp_tmp_tmp * a_tmp *
                                 (g_tmp_tmp - 1.0) / b_Smod_tmp) -
                            4.0 * r * c_del_tmp * a_tmp * (g_tmp_tmp - 1.0) *
                                c_fh3_tmp / d_Smod_tmp_tmp) +
                           g_fh3_tmp * c_fh3_tmp / e_fh3_tmp) /
                          e_del_tmp +
                      h_Smod_tmp *
                          ((b_fh3_tmp_tmp + fh3_tmp_tmp * g_tmp_tmp) -
                           b_Smod_tmp_tmp * c_fh3_tmp / b_a_tmp) /
                          e_del_tmp) -
                     g_Smod_tmp * c_fh3_tmp * h_Smod_tmp / b_Smod_tmp) -
                    d_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
                   c_del_tmp * g_tmp * a_tmp * c_fh3_tmp / e_fh3_tmp) +
                  h_fh3_tmp * b_g_tmp * a_tmp * e_Smod_tmp * c_fh3_tmp /
                      d_Smod_tmp_tmp) -
                 del) -
                f_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
               Dq) -
              g_fh3_tmp * e_Smod_tmp * c_fh3_tmp / m_fh3_tmp) -
             n_fh3_tmp) +
            o_fh3_tmp) +
           p_fh3_tmp) -
          theta;
      fh3[6] = 0.0;
      d_fh3_tmp = del_tmp * b_g_tmp * c_fh3_tmp;
      fh3[12] = (r * fh3_tmp_tmp_tmp * (g_tmp_tmp - 1.0) / e_del_tmp -
                 d_fh3_tmp / i_fh3_tmp) -
                Smod_tmp * c_fh3_tmp / b_Smod_tmp;
      fh3[1] = 0.0;
      f_fh3_tmp = b_fh3_tmp * d_del_tmp * b_g_tmp;
      g_fh3_tmp = b_fh3_tmp_tmp * b_g_tmp * a_tmp;
      h_fh3_tmp = g_Smod_tmp_tmp * a_tmp;
      fh3[7] =
          ((((((((((((j_Smod_tmp *
                          ((((((((j_fh3_tmp + f_fh3_tmp / c_Smod_tmp) -
                                 5.0 * d_del_tmp * b_g_tmp * a_tmp * c_fh3_tmp /
                                     k_fh3_tmp) +
                                g_fh3_tmp / c_Smod_tmp) +
                               l_fh3_tmp * d_del_tmp * (g_tmp_tmp - 1.0) /
                                   b_Smod_tmp) +
                              fh3_tmp) +
                             4.0 * r * b_del_tmp * b_fh3_tmp_tmp_tmp * a_tmp *
                                 (g_tmp_tmp - 1.0) / b_Smod_tmp) -
                            4.0 * r * d_del_tmp * a_tmp * (g_tmp_tmp - 1.0) *
                                c_fh3_tmp / d_Smod_tmp_tmp) +
                           h_fh3_tmp * c_fh3_tmp / e_fh3_tmp) /
                          e_del_tmp +
                      d_Smod_tmp *
                          ((fh3_tmp_tmp + b_fh3_tmp_tmp * g_tmp_tmp) -
                           f_Smod_tmp_tmp * c_fh3_tmp / b_a_tmp) /
                          e_del_tmp) -
                     j_Smod_tmp * c_fh3_tmp * d_Smod_tmp / b_Smod_tmp) -
                    f_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
                   d_del_tmp * g_tmp * a_tmp * c_fh3_tmp / e_fh3_tmp) +
                  2.0 * d_del_tmp * b_g_tmp * a_tmp * e_Smod_tmp * c_fh3_tmp /
                      d_Smod_tmp_tmp) -
                 del) -
                g_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
               Dq) -
              h_fh3_tmp * e_Smod_tmp * c_fh3_tmp / m_fh3_tmp) -
             n_fh3_tmp) +
            o_fh3_tmp) +
           p_fh3_tmp) -
          theta;
      e_fh3_tmp = b_del_tmp * b_g_tmp * c_fh3_tmp;
      fh3[13] = (r * b_fh3_tmp_tmp_tmp * (g_tmp_tmp - 1.0) / e_del_tmp -
                 e_fh3_tmp / i_fh3_tmp) -
                f_Smod_tmp * c_fh3_tmp / b_Smod_tmp;
      theta = c_Smod_tmp_tmp * c_fh3_tmp;
      f_fh3_tmp = (d_fh3_tmp_tmp * f_del_tmp / 2.0 -
                   3.0 * c_fh3_tmp * f_del_tmp / 2.0) +
                  theta;
      g_fh3_tmp = 2.0 * c_fh3_tmp_tmp;
      fh3[2] = ((fh3_tmp_tmp_tmp * a_tmp * i_Smod_tmp / k_Smod_tmp -
                 l_Smod_tmp * f_fh3_tmp / k_Smod_tmp) +
                del_tmp * b_fh3_tmp * i_Smod_tmp / k_Smod_tmp) -
               5.0 * del_tmp * a_tmp * i_Smod_tmp * c_fh3_tmp / g_fh3_tmp;
      fh3[8] = ((b_fh3_tmp_tmp_tmp * a_tmp * i_Smod_tmp / k_Smod_tmp -
                 m_Smod_tmp * f_fh3_tmp / k_Smod_tmp) +
                b_del_tmp * b_fh3_tmp * i_Smod_tmp / k_Smod_tmp) -
               5.0 * b_del_tmp * a_tmp * i_Smod_tmp * c_fh3_tmp / g_fh3_tmp;
      fh3[14] = d_fh3_tmp_tmp / (2.0 * e_del_tmp) - theta / i_fh3_tmp;
      b_fh3_tmp = b_g_tmp_tmp *
                  (c_fh3_tmp / e_fh3_tmp_tmp - d_fh3_tmp_tmp / e_fh3_tmp_tmp) /
                  n_Smod_tmp;
      f_fh3_tmp = del_tmp * b_fh3_tmp_tmp_tmp * h_Smod_tmp_tmp / n_Smod_tmp;
      g_fh3_tmp = fh3_tmp_tmp_tmp * b_del_tmp * h_Smod_tmp_tmp / n_Smod_tmp;
      h_fh3_tmp =
          3.0 * del_tmp * b_del_tmp * h_Smod_tmp_tmp * c_fh3_tmp / m_fh3_tmp;
      fh3[3] = ((h_fh3_tmp - f_fh3_tmp) - g_fh3_tmp) - b_fh3_tmp;
      fh3[9] = 3.0 * p_Smod_tmp * c_fh3_tmp / m_fh3_tmp -
               (((d_del_tmp * c_fh3_tmp / e_fh3_tmp_tmp +
                  b_fh3_tmp_tmp * f_del_tmp) +
                 2.0 * r * del_tmp * fh3_tmp_tmp_tmp * b_g_tmp) +
                e_Smod_tmp_tmp * c_fh3_tmp / e_fh3_tmp_tmp) /
                   n_Smod_tmp;
      fh3[15] = 0.0;
      fh3[4] =
          (((c_del_tmp * c_fh3_tmp / e_fh3_tmp_tmp + fh3_tmp_tmp * f_del_tmp) +
            2.0 * r * b_del_tmp * b_fh3_tmp_tmp_tmp * b_g_tmp) +
           g_Smod_tmp_tmp * c_fh3_tmp / e_fh3_tmp_tmp) /
              n_Smod_tmp -
          3.0 * q_Smod_tmp * c_fh3_tmp / m_fh3_tmp;
      fh3[10] = ((b_fh3_tmp + f_fh3_tmp) + g_fh3_tmp) - h_fh3_tmp;
      fh3[16] = 0.0;
      b_fh3_tmp = 2.0 * r * c_Smod_tmp;
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
    for (b_i = 0; b_i < 3; b_i++) {
      skw[3 * b_i] = g[b_i];
      skw[3 * b_i + 1] = g[b_i + 4];
      skw[3 * b_i + 2] = g[b_i + 8];
    }
    for (b_i = 0; b_i < 9; b_i++) {
      b_skw[b_i] = -skw[b_i];
    }
    Smod[3] = -g[14];
    Smod[6] = g[13];
    Smod[1] = g[14];
    Smod[7] = -g[12];
    Smod[2] = -g[13];
    Smod[5] = g[12];
    for (b_i = 0; b_i < 3; b_i++) {
      Smod_tmp_tmp = b_skw[b_i];
      fh3_tmp = b_skw[b_i + 3];
      del = b_skw[b_i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        Al[b_i + 3 * i1] =
            (Smod_tmp_tmp * Smod[3 * i1] + fh3_tmp * Smod[3 * i1 + 1]) +
            del * Smod[3 * i1 + 2];
        XJ[(i1 + 6 * b_i) + 36 * (i + 1)] = skw[i1 + 3 * b_i];
      }
    }
    XJ_tmp = 36 * (i + 1);
    for (b_i = 0; b_i < 3; b_i++) {
      b_XJ_tmp = 6 * (b_i + 3) + XJ_tmp;
      XJ[b_XJ_tmp] = Al[3 * b_i];
      c_XJ_tmp = 6 * b_i + XJ_tmp;
      XJ[c_XJ_tmp + 3] = 0.0;
      XJ[b_XJ_tmp + 3] = skw[3 * b_i];
      j = 3 * b_i + 1;
      XJ[b_XJ_tmp + 1] = Al[j];
      XJ[c_XJ_tmp + 4] = 0.0;
      XJ[b_XJ_tmp + 4] = skw[j];
      j = 3 * b_i + 2;
      XJ[b_XJ_tmp + 2] = Al[j];
      XJ[c_XJ_tmp + 5] = 0.0;
      XJ[b_XJ_tmp + 5] = skw[j];
    }
    for (b_i = 0; b_i < 6; b_i++) {
      for (i1 = 0; i1 < 6; i1++) {
        Smod_tmp_tmp = 0.0;
        for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
          Smod_tmp_tmp += XJ[(b_i + 6 * b_XJ_tmp) + 36 * (i + 1)] *
                          (double)Xtree[b_XJ_tmp + 6 * i1];
        }
        Xup[(b_i + 6 * i1) + 36 * (i + 1)] = Smod_tmp_tmp;
      }
    }
    for (b_i = 0; b_i < 36; b_i++) {
      Xtree[b_i] = 0;
    }
    b_a_tmp = dq[qi_data[0] - 1];
    Dq = dq[qi_data[1] - 1];
    g_tmp = dq[qi_data[2] - 1];
    for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
      Xtree[XJ_tmp + 6 * XJ_tmp] = 1;
      b_i = XJ_tmp + 18 * i;
      Smod_tmp_tmp = (b_Smod[b_i] * b_a_tmp + b_Smod[b_i + 6] * Dq) +
                     b_Smod[b_i + 12] * g_tmp;
      vJ[XJ_tmp] = Smod_tmp_tmp;
      fh3_tmp = 0.0;
      for (b_i = 0; b_i < 6; b_i++) {
        fh3_tmp += Xup[(XJ_tmp + 6 * b_i) + 36 * (i + 1)] * b_v[b_i + 6 * i];
      }
      b_Xup[XJ_tmp] = fh3_tmp + Smod_tmp_tmp;
    }
    for (b_i = 0; b_i < 6; b_i++) {
      b_v[b_i + 6 * (i + 1)] = b_Xup[b_i];
    }
    memset(&I_tmp[0], 0, 36U * sizeof(double));
    skw[0] = 0.0;
    b_i = 6 * (i + 1);
    Smod_tmp_tmp = b_v[b_i + 5];
    skw[3] = -Smod_tmp_tmp;
    fh3_tmp = b_v[b_i + 4];
    skw[6] = fh3_tmp;
    skw[1] = Smod_tmp_tmp;
    skw[4] = 0.0;
    del = b_v[b_i + 3];
    skw[7] = -del;
    skw[2] = -fh3_tmp;
    skw[5] = del;
    skw[8] = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      theta = skw[3 * i1];
      I_tmp[6 * i1] = theta;
      XJ_tmp = 6 * (i1 + 3);
      I_tmp[XJ_tmp + 3] = theta;
      theta = skw[3 * i1 + 1];
      I_tmp[6 * i1 + 1] = theta;
      I_tmp[XJ_tmp + 4] = theta;
      theta = skw[3 * i1 + 2];
      I_tmp[6 * i1 + 2] = theta;
      I_tmp[XJ_tmp + 5] = theta;
    }
    I_tmp[18] = 0.0;
    n_fh3_tmp = b_v[b_i + 2];
    I_tmp[24] = -n_fh3_tmp;
    o_fh3_tmp = b_v[b_i + 1];
    I_tmp[30] = o_fh3_tmp;
    I_tmp[19] = n_fh3_tmp;
    I_tmp[25] = 0.0;
    p_fh3_tmp = b_v[b_i];
    I_tmp[31] = -p_fh3_tmp;
    I_tmp[20] = -o_fh3_tmp;
    I_tmp[26] = p_fh3_tmp;
    I_tmp[32] = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      theta = 0.0;
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        theta += Xup[(i1 + 6 * b_XJ_tmp) + 36 * (i + 1)] * a[b_XJ_tmp + 6 * i];
      }
      b_Xup[i1] = theta;
      c_I[i1] = (fh3[i1] * b_a_tmp + fh3[i1 + 6] * Dq) + fh3[i1 + 12] * g_tmp;
    }
    for (i1 = 0; i1 < 6; i1++) {
      theta = 0.0;
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        theta += I_tmp[i1 + 6 * b_XJ_tmp] * vJ[b_XJ_tmp];
      }
      a[i1 + b_i] = (b_Xup[i1] + c_I[i1]) + theta;
    }
    memset(&I_tmp[0], 0, 36U * sizeof(double));
    skw[0] = 0.0;
    skw[3] = -Smod_tmp_tmp;
    skw[6] = fh3_tmp;
    skw[1] = Smod_tmp_tmp;
    skw[4] = 0.0;
    skw[7] = -del;
    skw[2] = -fh3_tmp;
    skw[5] = del;
    skw[8] = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      theta = skw[3 * i1];
      I_tmp[6 * i1] = theta;
      XJ_tmp = 6 * (i1 + 3);
      I_tmp[XJ_tmp + 3] = theta;
      theta = skw[3 * i1 + 1];
      I_tmp[6 * i1 + 1] = theta;
      I_tmp[XJ_tmp + 4] = theta;
      theta = skw[3 * i1 + 2];
      I_tmp[6 * i1 + 2] = theta;
      I_tmp[XJ_tmp + 5] = theta;
    }
    I_tmp[18] = 0.0;
    I_tmp[24] = -n_fh3_tmp;
    I_tmp[30] = o_fh3_tmp;
    I_tmp[19] = n_fh3_tmp;
    I_tmp[25] = 0.0;
    I_tmp[31] = -p_fh3_tmp;
    I_tmp[20] = -o_fh3_tmp;
    I_tmp[26] = p_fh3_tmp;
    I_tmp[32] = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        b_I_tmp[b_XJ_tmp + 6 * i1] = -I_tmp[i1 + 6 * b_XJ_tmp];
      }
    }
    for (i1 = 0; i1 < 6; i1++) {
      c_I[i1] = 0.0;
      vJ[i1] = 0.0;
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        Smod_tmp_tmp = 0.0;
        for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
          Smod_tmp_tmp += b_I_tmp[i1 + 6 * XJ_tmp] *
                          b_I[(XJ_tmp + 6 * b_XJ_tmp) + 36 * (i + 1)];
        }
        XJ_tmp = b_XJ_tmp + b_i;
        c_I[i1] += b_I[(i1 + 6 * b_XJ_tmp) + 36 * (i + 1)] * a[XJ_tmp];
        vJ[i1] += Smod_tmp_tmp * b_v[XJ_tmp];
      }
      f[i1 + b_i] = c_I[i1] + vJ[i1];
    }
  }
  /*  Composite Rigid Body Algorithm to calculate M */
  /*  composite inertia calculation */
  for (i = 0; i < 4; i++) {
    if (4 - i == 1) {
      del = 0.0;
      for (b_i = 0; b_i < 6; b_i++) {
        del += (double)c_a[b_i] * f[b_i];
      }
      C[0] = del;
    } else {
      XJ_tmp = (2 - i) * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        qi_data[b_i] = (signed char)((XJ_tmp + b_i) + 2);
      }
      for (b_i = 0; b_i < 3; b_i++) {
        tmp_data[b_i] = (signed char)(qi_data[b_i] - 1);
      }
      for (b_i = 0; b_i < 3; b_i++) {
        Smod_tmp_tmp = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          Smod_tmp_tmp +=
              2.0 * b_Smod[(i1 + 6 * b_i) + 18 * (2 - i)] * f[i1 + 6 * (3 - i)];
        }
        b_dq[b_i] = Smod_tmp_tmp;
      }
      for (b_i = 0; b_i < 3; b_i++) {
        C[tmp_data[b_i]] = b_dq[b_i];
      }
      for (b_i = 0; b_i < 6; b_i++) {
        Smod_tmp_tmp = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          Smod_tmp_tmp +=
              Xup[(i1 + 6 * b_i) + 36 * (3 - i)] * f[i1 + 6 * (3 - i)];
        }
        c_I[b_i] = f[b_i + 6 * (2 - i)] + Smod_tmp_tmp;
      }
      for (b_i = 0; b_i < 6; b_i++) {
        f[b_i + 6 * (2 - i)] = c_I[b_i];
      }
    }
    if (4 - i != 1) {
      b_i = 36 * (3 - i);
      for (i1 = 0; i1 < 6; i1++) {
        for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
          Smod_tmp_tmp = 0.0;
          for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
            Smod_tmp_tmp += Xup[(XJ_tmp + 6 * i1) + b_i] *
                            b_I[(XJ_tmp + 6 * b_XJ_tmp) + b_i];
          }
          I_tmp[i1 + 6 * b_XJ_tmp] = Smod_tmp_tmp;
        }
        for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
          Smod_tmp_tmp = 0.0;
          for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
            Smod_tmp_tmp += I_tmp[i1 + 6 * XJ_tmp] *
                            Xup[(XJ_tmp + 6 * b_XJ_tmp) + 36 * (3 - i)];
          }
          XJ_tmp = i1 + 6 * b_XJ_tmp;
          b_I_tmp[XJ_tmp] = b_I[XJ_tmp + 36 * (2 - i)] + Smod_tmp_tmp;
        }
      }
      memcpy(&b_I[i * -36 + 72], &b_I_tmp[0], 36U * sizeof(double));
    }
  }
  memset(&M[0], 0, 100U * sizeof(double));
  /*  fh3 = zeros(6,3); */
  /*  fh1 = zeros(6,1); */
  for (i = 0; i < 4; i++) {
    if (i + 1 == 1) {
      Smod_tmp_tmp = 0.0;
      for (b_i = 0; b_i < 6; b_i++) {
        fh3_tmp = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          fh3_tmp += b_I[b_i + 6 * i1] * (double)b_a[i1];
        }
        Smod_tmp_tmp += (double)iv[b_i] * fh3_tmp;
      }
      M[0] = Smod_tmp_tmp;
    } else {
      signed char b_tmp_data[9];
      XJ_tmp = (i - 1) * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        qi_data[b_i] = (signed char)((XJ_tmp + b_i) + 2);
      }
      for (b_i = 0; b_i < 6; b_i++) {
        for (i1 = 0; i1 < 3; i1++) {
          Smod_tmp_tmp = 0.0;
          for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
            Smod_tmp_tmp += b_I[(b_i + 6 * b_XJ_tmp) + 36 * i] *
                            b_Smod[(b_XJ_tmp + 6 * i1) + 18 * (i - 1)];
          }
          fh3[b_i + 6 * i1] = Smod_tmp_tmp;
        }
      }
      for (b_i = 0; b_i < 3; b_i++) {
        i1 = qi_data[b_i];
        tmp_data[b_i] = (signed char)(i1 - 1);
        b_tmp_data[b_i] = (signed char)(i1 - 1);
      }
      for (b_i = 0; b_i < 3; b_i++) {
        for (i1 = 0; i1 < 3; i1++) {
          Smod_tmp_tmp = 0.0;
          for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
            Smod_tmp_tmp += b_Smod[(b_XJ_tmp + 6 * b_i) + 18 * (i - 1)] *
                            fh3[b_XJ_tmp + 6 * i1];
          }
          Smod[b_i + 3 * i1] = Smod_tmp_tmp;
        }
      }
      for (b_i = 0; b_i < 3; b_i++) {
        for (i1 = 0; i1 < 3; i1++) {
          M[tmp_data[i1] + 10 * b_tmp_data[b_i]] = Smod[i1 + 3 * b_i];
        }
      }
      j = i - 1;
      while (j + 2 > 1) {
        for (b_i = 0; b_i < 6; b_i++) {
          for (i1 = 0; i1 < 3; i1++) {
            Smod_tmp_tmp = 0.0;
            for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
              Smod_tmp_tmp += Xup[(b_XJ_tmp + 6 * b_i) + 36 * (j + 1)] *
                              fh3[b_XJ_tmp + 6 * i1];
            }
            c_Xup[b_i + 6 * i1] = Smod_tmp_tmp;
          }
        }
        memcpy(&fh3[0], &c_Xup[0], 18U * sizeof(double));
        j--;
        if (j + 2 == 1) {
          for (b_i = 0; b_i < 3; b_i++) {
            Smod_tmp_tmp = 0.0;
            for (i1 = 0; i1 < 6; i1++) {
              Smod_tmp_tmp += fh3[i1 + 6 * b_i] * (double)b_a[i1];
            }
            b_dq[b_i] = Smod_tmp_tmp;
          }
          for (b_i = 0; b_i < 3; b_i++) {
            M[b_tmp_data[b_i]] = b_dq[b_i];
          }
          /* (S{j}' * fh).'; */
          for (b_i = 0; b_i < 3; b_i++) {
            Smod_tmp_tmp = 0.0;
            for (i1 = 0; i1 < 6; i1++) {
              Smod_tmp_tmp += (double)iv[i1] * fh3[i1 + 6 * b_i];
            }
            b_dq[b_i] = Smod_tmp_tmp;
          }
          for (b_i = 0; b_i < 3; b_i++) {
            M[10 * b_tmp_data[b_i]] = b_dq[b_i];
          }
        } else {
          signed char c_tmp_data[6];
          signed char d_tmp_data[6];
          b_XJ_tmp = j * 3 + 2;
          XJ_tmp = j * 3 + 4;
          if (XJ_tmp < b_XJ_tmp) {
            c_XJ_tmp = 0;
          } else {
            XJ_tmp = (signed char)XJ_tmp - (signed char)b_XJ_tmp;
            c_XJ_tmp = XJ_tmp + 1;
            for (b_i = 0; b_i <= XJ_tmp; b_i++) {
              qi_data[b_i] =
                  (signed char)((signed char)b_XJ_tmp + (signed char)b_i);
            }
          }
          for (b_i = 0; b_i < c_XJ_tmp; b_i++) {
            i1 = qi_data[b_i];
            c_tmp_data[b_i] = (signed char)i1;
            d_tmp_data[b_i] = (signed char)(i1 - 1);
          }
          for (b_i = 0; b_i < 3; b_i++) {
            for (i1 = 0; i1 < 3; i1++) {
              Smod_tmp_tmp = 0.0;
              for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
                Smod_tmp_tmp += fh3[b_XJ_tmp + 6 * b_i] *
                                b_Smod[(b_XJ_tmp + 6 * i1) + 18 * j];
              }
              Smod[b_i + 3 * i1] = Smod_tmp_tmp;
            }
          }
          /* (S{j}' * fh).'; */
          for (b_i = 0; b_i < c_XJ_tmp; b_i++) {
            for (i1 = 0; i1 < 3; i1++) {
              M[b_tmp_data[i1] + 10 * d_tmp_data[b_i]] = Smod[i1 + 3 * b_i];
            }
            d_tmp_data[b_i] = (signed char)(c_tmp_data[b_i] - 1);
          }
          for (b_i = 0; b_i < 3; b_i++) {
            for (i1 = 0; i1 < 3; i1++) {
              Smod_tmp_tmp = 0.0;
              for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
                Smod_tmp_tmp += b_Smod[(b_XJ_tmp + 6 * b_i) + 18 * j] *
                                fh3[b_XJ_tmp + 6 * i1];
              }
              Smod[b_i + 3 * i1] = Smod_tmp_tmp;
            }
          }
          for (b_i = 0; b_i < 3; b_i++) {
            for (i1 = 0; i1 < c_XJ_tmp; i1++) {
              M[d_tmp_data[i1] + 10 * b_tmp_data[b_i]] =
                  Smod[i1 + c_XJ_tmp * b_i];
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
  b_C[0] = 0.0;
  memcpy(&b_C[1], &skw[0], 9U * sizeof(double));
  memset(&K[0], 0, 100U * sizeof(double));
  for (j = 0; j < 10; j++) {
    K[j + 10 * j] = b_C[j];
  }
  b_dq[0] = bb;
  b_dq[1] = bb;
  b_dq[2] = bs;
  repmat(b_dq, skw);
  b_C[0] = bm;
  memcpy(&b_C[1], &skw[0], 9U * sizeof(double));
  memset(&D[0], 0, 100U * sizeof(double));
  for (j = 0; j < 10; j++) {
    D[j + 10 * j] = b_C[j];
  }
  /*  Fci = cell(1,4); */
  /*  Ai = cell(1,4); */
  for (i = 0; i < 3; i++) {
    Smod_tmp_tmp = qp[3 * i];
    fh3_tmp = qp[3 * i + 1];
    del_tmp = Smod_tmp_tmp * Smod_tmp_tmp;
    b_del_tmp = fh3_tmp * fh3_tmp;
    del = sqrt(del_tmp + b_del_tmp);
    theta = del / d;
    b_a_tmp = sin(del);
    Dq = del - b_a_tmp;
    g_tmp = qp[3 * i + 2] + L0;
    if ((Smod_tmp_tmp < 1.0E-5) && (fh3_tmp < 1.0E-5)) {
      cq[i] = g_tmp / 3.0;
      for (b_i = 0; b_i < 9; b_i++) {
        skw[b_i] = iv1[b_i];
      }
    } else {
      cq[i] = 2.0 * (g_tmp / theta - d) * sin(theta / 6.0);
      theta = rt_powd_snf(del, 3.0);
      skw[0] = Smod_tmp_tmp * fh3_tmp * Dq / theta;
      skw[3] = (-del_tmp * del - b_del_tmp * b_a_tmp) / theta;
      skw[6] = Smod_tmp_tmp * Dq * g_tmp / theta;
      skw[1] = (b_del_tmp * del + del_tmp * b_a_tmp) / theta;
      skw[4] = -Smod_tmp_tmp * fh3_tmp * Dq / theta;
      skw[7] = fh3_tmp * Dq * g_tmp / theta;
      skw[2] = 0.0;
      skw[5] = 0.0;
      skw[8] = b_a_tmp / del;
    }
    /*  Al = [d*cosd(60) d*cosd(60) -d; */
    /*  -d*cosd(30) d*cosd(30) 0; */
    /*  1 1 1]; */
    if (i + 1 == 2) {
      Smod[0] = -d;
      Smod[3] = d * 0.86602540378443871;
      Smod[6] = d * 0.86602540378443871;
      Smod[1] = 0.0;
      Smod[4] = d * 0.49999999999999994;
      Smod[7] = -d * 0.49999999999999994;
      Smod[2] = 1.0;
      Smod[5] = 1.0;
      Smod[8] = 1.0;
      for (b_i = 0; b_i < 3; b_i++) {
        signed char i2;
        signed char i3;
        signed char i4;
        i2 = d_a[b_i];
        i3 = d_a[b_i + 3];
        i4 = d_a[b_i + 6];
        for (i1 = 0; i1 < 3; i1++) {
          Al[b_i + 3 * i1] =
              ((double)i2 * Smod[3 * i1] + (double)i3 * Smod[3 * i1 + 1]) +
              (double)i4 * Smod[3 * i1 + 2];
        }
      }
    } else {
      Al[0] = -d;
      Al[3] = d * 0.86602540378443871;
      Al[6] = d * 0.86602540378443871;
      Al[1] = 0.0;
      Al[4] = d * 0.49999999999999994;
      Al[7] = -d * 0.49999999999999994;
      Al[2] = 1.0;
      Al[5] = 1.0;
      Al[8] = 1.0;
    }
    for (b_i = 0; b_i < 3; b_i++) {
      Smod_tmp_tmp = skw[b_i];
      fh3_tmp = skw[b_i + 3];
      del = skw[b_i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        b_skw[b_i + 3 * i1] = (Smod_tmp_tmp * (double)c_b[3 * i1] +
                               fh3_tmp * (double)c_b[3 * i1 + 1]) +
                              del * (double)c_b[3 * i1 + 2];
      }
      Smod_tmp_tmp = b_skw[b_i];
      fh3_tmp = b_skw[b_i + 3];
      del = b_skw[b_i + 6];
      for (i1 = 0; i1 < 3; i1++) {
        Ai[(b_i + 3 * i1) + 9 * i] =
            (Smod_tmp_tmp * Al[3 * i1] + fh3_tmp * Al[3 * i1 + 1]) +
            del * Al[3 * i1 + 2];
      }
    }
  }
  /*  A = blkdiag(1,Ai(:,:,1),Ai(:,:,2),Ai(:,:,3)); */
  memset(&A[0], 0, 100U * sizeof(double));
  A[0] = 1.0;
  b_dq[0] = conv_pcc;
  b_dq[1] = conv_pcc;
  b_dq[2] = conv_pcc;
  for (i = 0; i < 3; i++) {
    b_i = 3 * i;
    for (XJ_tmp = 0; XJ_tmp < 3; XJ_tmp++) {
      b_XJ_tmp = 3 * XJ_tmp + 9 * i;
      c_XJ_tmp = XJ_tmp + b_i;
      j = b_i + 10 * (c_XJ_tmp + 1);
      A[j + 1] = Ai[b_XJ_tmp];
      A[j + 2] = Ai[b_XJ_tmp + 1];
      A[j + 3] = Ai[b_XJ_tmp + 2];
      skw[c_XJ_tmp] = b_dq[XJ_tmp];
    }
  }
  b_C[0] = conv_motor;
  memcpy(&b_C[1], &skw[0], 9U * sizeof(double));
  memset(&conversion[0], 0, 100U * sizeof(double));
  for (j = 0; j < 10; j++) {
    conversion[j + 10 * j] = b_C[j];
    Smod_tmp_tmp = 0.0;
    b_qd[j] = qd[j] - q[j];
    fh3_tmp = 0.0;
    del = 0.0;
    for (b_i = 0; b_i < 10; b_i++) {
      i1 = j + 10 * b_i;
      Smod_tmp_tmp += M[i1] * ddqd[b_i];
      del += K[i1] * qd[b_i];
      fh3_tmp += D[i1] * dqd[b_i];
    }
    Smod_tmp_tmp += C[j];
    b_C[j] = Smod_tmp_tmp;
    c_C[j] = (Smod_tmp_tmp + del) + fh3_tmp;
  }
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      Smod_tmp_tmp += Kp[b_i + 10 * i1] * b_qd[i1];
    }
    b_K[b_i] = Smod_tmp_tmp;
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_qd[b_i] = dqd[b_i] - dq[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      Smod_tmp_tmp += KD[b_i + 10 * i1] * b_qd[i1];
    }
    b_C[b_i] = (c_C[b_i] + b_K[b_i]) + Smod_tmp_tmp;
  }
  mldivide(A, b_C);
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      Smod_tmp_tmp += conversion[b_i + 10 * i1] * b_C[i1];
    }
    tau[b_i] = Smod_tmp_tmp;
  }
  memcpy(&Kpr[0], &Kp[0], 100U * sizeof(double));
  Kpr[0] = 0.0;
  memcpy(&KDr[0], &KD[0], 100U * sizeof(double));
  KDr[0] = 0.0;
  memcpy(&c_C[0], &q[0], 10U * sizeof(double));
  J_r(c_C, L0, d, J, x);
  for (b_i = 0; b_i < 3; b_i++) {
    for (i1 = 0; i1 < 10; i1++) {
      b_tmp_tmp[i1 + 10 * b_i] = J[b_i + 3 * i1];
    }
  }
  memcpy(&b_tmp[0], &b_tmp_tmp[0], 30U * sizeof(double));
  b_mldivide(M, b_tmp);
  for (b_i = 0; b_i < 3; b_i++) {
    for (i1 = 0; i1 < 3; i1++) {
      Smod_tmp_tmp = 0.0;
      for (b_XJ_tmp = 0; b_XJ_tmp < 10; b_XJ_tmp++) {
        Smod_tmp_tmp += J[b_i + 3 * b_XJ_tmp] * b_tmp[b_XJ_tmp + 10 * i1];
      }
      Smod[b_i + 3 * i1] = Smod_tmp_tmp;
    }
  }
  mrdiv(b_A, Smod, skw);
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = b_tmp[b_i];
    fh3_tmp = b_tmp[b_i + 10];
    del = b_tmp[b_i + 20];
    for (i1 = 0; i1 < 3; i1++) {
      b_b_tmp[b_i + 10 * i1] =
          (Smod_tmp_tmp * skw[3 * i1] + fh3_tmp * skw[3 * i1 + 1]) +
          del * skw[3 * i1 + 2];
    }
  }
  mtimes(J, b_b_tmp, y_tmp);
  memset(&d_I[0], 0, 100U * sizeof(signed char));
  for (XJ_tmp = 0; XJ_tmp < 10; XJ_tmp++) {
    d_I[XJ_tmp + 10 * XJ_tmp] = 1;
    Smod_tmp_tmp = 0.0;
    fh3_tmp = 0.0;
    for (b_i = 0; b_i < 10; b_i++) {
      i1 = XJ_tmp + 10 * b_i;
      Smod_tmp_tmp += K[i1] * q[b_i];
      fh3_tmp += D[i1] * dq[b_i];
    }
    b_K[XJ_tmp] = Smod_tmp_tmp + fh3_tmp;
    Smod_tmp_tmp = b_tmp_tmp[XJ_tmp];
    fh3_tmp = b_tmp_tmp[XJ_tmp + 10];
    del = b_tmp_tmp[XJ_tmp + 20];
    for (b_i = 0; b_i < 3; b_i++) {
      b_tmp[XJ_tmp + 10 * b_i] =
          (Smod_tmp_tmp * skw[3 * b_i] + fh3_tmp * skw[3 * b_i + 1]) +
          del * skw[3 * b_i + 2];
    }
  }
  b_dq[0] = Kpx * (xd[0] - x[0]) + KDx * (dxd[0] - dxr[0]);
  b_dq[1] = Kpx * (xd[1] - x[1]) + KDx * (dxd[1] - dxr[1]);
  b_dq[2] = Kpx * (xd[2] - x[2]) + KDx * (dxd[2] - dxr[2]);
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      Smod_tmp_tmp += y_tmp[b_i + 10 * i1] * b_K[i1];
    }
    b_C[b_i] = (C[b_i] + Smod_tmp_tmp) +
               ((b_tmp[b_i] * b_dq[0] + b_tmp[b_i + 10] * b_dq[1]) +
                b_tmp[b_i + 20] * b_dq[2]);
  }
  mldivide(A, b_C);
  for (b_i = 0; b_i < 100; b_i++) {
    Kpr[b_i] = -Kpr[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = 0.0;
    fh3_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      b_XJ_tmp = b_i + 10 * i1;
      Smod_tmp_tmp += Kpr[b_XJ_tmp] * q[i1];
      fh3_tmp += KDr[b_XJ_tmp] * dq[i1];
    }
    b_qd[b_i] = fh3_tmp;
    b_K[b_i] = Smod_tmp_tmp;
  }
  for (b_i = 0; b_i < 100; b_i++) {
    y_tmp[b_i] = (double)d_I[b_i] - y_tmp[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    b_K[b_i] -= b_qd[b_i];
  }
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      Smod_tmp_tmp += y_tmp[b_i + 10 * i1] * b_K[i1];
    }
    b_qd[b_i] = b_C[b_i] + Smod_tmp_tmp;
  }
  for (b_i = 0; b_i < 10; b_i++) {
    Smod_tmp_tmp = 0.0;
    for (i1 = 0; i1 < 10; i1++) {
      Smod_tmp_tmp += conversion[b_i + 10 * i1] * b_qd[i1];
    }
    tau_r[b_i] = Smod_tmp_tmp;
  }
}

/*
 * File trailer for helix_controller.c
 *
 * [EOF]
 */
