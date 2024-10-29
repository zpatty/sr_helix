/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: MC_3_cg.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 29-Oct-2024 18:41:09
 */

/* Include Files */
#include "MC_3_cg.h"
#include "helix_controller_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * q = state
 *  qd = velocity
 *  qdd = acceleration
 *  m = mass of a module;
 *  mm = mass of a motor;
 *  hm = height of a motor;
 *  rm = "radius" a motor;
 *  r = radius of the hex plate;
 *  L0 = Length of a module;
 *  d = distance to cable;
 *  N = number of links (number of modules plus number of motors)
 *
 * Arguments    : const double q[10]
 *                const double qd[10]
 *                double m
 *                double r
 *                double L0
 *                double M[100]
 *                double C[10]
 * Return Type  : void
 */
void MC_3_cg(const double q[10], const double qd[10], double m, double r,
             double L0, double M[100], double C[10])
{
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
  static const signed char c_a[6] = {0, 0, 0, 0, 0, 1};
  static const signed char d_a[6] = {0, 0, 0, 0, 0, 2};
  static const signed char iv[6] = {0, 0, 0, 0, 0, 1};
  double XJ[144];
  double Xup[144];
  double b_I[144];
  double Smod[54];
  double I_tmp[36];
  double b_I_tmp[36];
  double a[24];
  double b_v[24];
  double f[24];
  double c_Xup[18];
  double fh3[18];
  double g[16];
  double b_skw[9];
  double c_skw[9];
  double dv[9];
  double skw[9];
  double c_I[6];
  double vJ[6];
  double b_qd[3];
  double a_tmp;
  double b_a;
  double b_g_tmp;
  double g_tmp;
  double g_tmp_tmp;
  double vJ_tmp;
  double vJ_tmp_tmp;
  int XJ_tmp;
  int b_XJ_tmp;
  int b_i;
  int c_XJ_tmp;
  int i;
  int i1;
  int j;
  signed char Xtree[36];
  signed char qi_data[9];
  signed char tmp_data[9];
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
  vJ_tmp = 0.25 * m * vJ_tmp_tmp;
  vJ[3] = vJ_tmp;
  vJ[4] = vJ_tmp;
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
  dv[0] = 0.0;
  dv[3] = -0.0;
  dv[6] = 0.0;
  dv[1] = 0.0;
  dv[4] = 0.0;
  dv[7] = -0.0;
  dv[2] = -0.0;
  dv[5] = 0.0;
  dv[8] = 0.0;
  for (b_i = 0; b_i < 3; b_i++) {
    b_a = b_skw[b_i];
    g_tmp_tmp = b_skw[b_i + 3];
    a_tmp = b_skw[b_i + 6];
    for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
      c_skw[b_i + 3 * b_XJ_tmp] =
          (b_a * dv[3 * b_XJ_tmp] + g_tmp_tmp * dv[3 * b_XJ_tmp + 1]) +
          a_tmp * dv[3 * b_XJ_tmp + 2];
      XJ[b_XJ_tmp + 6 * b_i] = skw[b_XJ_tmp + 3 * b_i];
    }
  }
  for (b_i = 0; b_i < 3; b_i++) {
    XJ_tmp = 6 * (b_i + 3);
    XJ[XJ_tmp] = c_skw[3 * b_i];
    XJ[6 * b_i + 3] = 0.0;
    XJ[XJ_tmp + 3] = skw[3 * b_i];
    b_XJ_tmp = 3 * b_i + 1;
    XJ[XJ_tmp + 1] = c_skw[b_XJ_tmp];
    XJ[6 * b_i + 4] = 0.0;
    XJ[XJ_tmp + 4] = skw[b_XJ_tmp];
    b_XJ_tmp = 3 * b_i + 2;
    XJ[XJ_tmp + 2] = c_skw[b_XJ_tmp];
    XJ[6 * b_i + 5] = 0.0;
    XJ[XJ_tmp + 5] = skw[b_XJ_tmp];
  }
  for (i = 0; i < 6; i++) {
    for (b_i = 0; b_i < 6; b_i++) {
      b_a = 0.0;
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        b_a += XJ[i + 6 * b_XJ_tmp] * (double)b[b_XJ_tmp + 6 * b_i];
      }
      Xup[i + 6 * b_i] = b_a;
    }
    b_a = (double)c_a[i] * qd[0];
    vJ[i] = b_a;
    b_v[i] = b_a;
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
    vJ_tmp = skw[3 * b_i];
    I_tmp[6 * b_i] = vJ_tmp;
    XJ_tmp = 6 * (b_i + 3);
    I_tmp[XJ_tmp + 3] = vJ_tmp;
    vJ_tmp = skw[3 * b_i + 1];
    I_tmp[6 * b_i + 1] = vJ_tmp;
    I_tmp[XJ_tmp + 4] = vJ_tmp;
    vJ_tmp = skw[3 * b_i + 2];
    I_tmp[6 * b_i + 2] = vJ_tmp;
    I_tmp[XJ_tmp + 5] = vJ_tmp;
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
    b_a = 0.0;
    g_tmp_tmp = 0.0;
    for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
      i1 = b_i + 6 * b_XJ_tmp;
      b_a += Xup[i1] * b_b[b_XJ_tmp];
      a_tmp = I_tmp[i1];
      g_tmp_tmp += a_tmp * vJ[b_XJ_tmp];
      b_I_tmp[b_XJ_tmp + 6 * b_i] = -a_tmp;
    }
    a[b_i] = b_a + g_tmp_tmp;
  }
  for (b_i = 0; b_i < 6; b_i++) {
    c_I[b_i] = 0.0;
    vJ[b_i] = 0.0;
    for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
      b_a = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        b_a += b_I_tmp[b_i + 6 * i1] * b_I[i1 + 6 * b_XJ_tmp];
      }
      c_I[b_i] += b_I[b_i + 6 * b_XJ_tmp] * a[b_XJ_tmp];
      vJ[b_i] += b_a * b_v[b_XJ_tmp];
    }
    f[b_i] = c_I[b_i] + vJ[b_i];
  }
  /*  Recursive Newton Euler to Calculate C+G */
  dv[0] = 0.0;
  dv[4] = 0.0;
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    double b_Xup[6];
    double b_del_tmp;
    double c_del_tmp;
    double d_Smod_tmp_tmp;
    double d_del_tmp;
    double d_fh3_tmp;
    double del_tmp;
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
      double fh3_tmp;
      g[12] = 0.0;
      g[13] = 0.0;
      g_tmp_tmp = q[qi_data[2] - 1];
      g_tmp = L0 + g_tmp_tmp;
      g[14] = g_tmp;
      g[3] = 0.0;
      g[7] = 0.0;
      g[11] = 0.0;
      g[15] = 1.0;
      Smod_tmp = g_tmp / (2.0 * r);
      Smod[18 * i] = Smod_tmp;
      Smod[18 * i + 6] = 0.0;
      Smod[18 * i + 12] = 0.0;
      Smod[18 * i + 1] = 0.0;
      Smod[18 * i + 7] = Smod_tmp;
      Smod[18 * i + 13] = 0.0;
      Smod[18 * i + 3] = 0.0;
      Smod[18 * i + 9] = -1.0 / r;
      Smod[18 * i + 15] = 0.0;
      Smod[18 * i + 4] = 1.0 / r;
      Smod[18 * i + 10] = 0.0;
      Smod[18 * i + 16] = 0.0;
      fh3_tmp = qd[qi_data[2] - 1] / (2.0 * r);
      fh3[0] = fh3_tmp;
      fh3[6] = 0.0;
      b_fh3_tmp = qd[qi_data[0] - 1];
      fh3[12] = -b_fh3_tmp / (2.0 * r);
      fh3[1] = 0.0;
      fh3[7] = fh3_tmp;
      fh3_tmp = qd[qi_data[1] - 1];
      fh3[13] = -fh3_tmp / (2.0 * r);
      c_fh3_tmp = 6.0 * vJ_tmp_tmp;
      fh3[2] = b_fh3_tmp * g_tmp / c_fh3_tmp;
      fh3[8] = (L0 * fh3_tmp + fh3_tmp * g_tmp_tmp) / c_fh3_tmp;
      fh3[14] = 0.0;
      Smod[18 * i + 2] = 0.0;
      Smod[18 * i + 5] = 0.0;
      fh3[3] = 0.0;
      fh3[4] = 0.0;
      Smod[18 * i + 8] = 0.0;
      Smod[18 * i + 11] = 0.0;
      fh3[9] = 0.0;
      fh3[10] = 0.0;
      Smod[18 * i + 14] = 1.0;
      Smod[18 * i + 17] = 0.0;
      fh3[15] = 0.0;
      fh3[16] = 0.0;
      c_fh3_tmp = 2.0 * vJ_tmp_tmp;
      fh3[5] = fh3_tmp / c_fh3_tmp;
      fh3[11] = -b_fh3_tmp / c_fh3_tmp;
      fh3[17] = 0.0;
    } else {
      double Smod_tmp;
      double Smod_tmp_tmp;
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
      double c_g_tmp_tmp;
      double d_Smod_tmp;
      double d_fh3_tmp_tmp;
      double e_Smod_tmp;
      double e_Smod_tmp_tmp;
      double e_fh3_tmp;
      double e_fh3_tmp_tmp;
      double f_Smod_tmp;
      double f_Smod_tmp_tmp;
      double f_fh3_tmp;
      double fh3_tmp;
      double fh3_tmp_tmp;
      double fh3_tmp_tmp_tmp;
      double g_Smod_tmp;
      double g_Smod_tmp_tmp;
      double g_fh3_tmp;
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
      vJ_tmp = f_del_tmp * f_del_tmp;
      b_a_tmp = L0 + q[qi_data[2] - 1];
      a_tmp = r * b_a_tmp;
      b_a = a_tmp / vJ_tmp;
      g_tmp_tmp = f_del_tmp / r;
      b_g_tmp_tmp = cos(g_tmp_tmp);
      g[0] = c_del_tmp / vJ_tmp * (b_g_tmp_tmp - 1.0) + 1.0;
      c_g_tmp_tmp = del_tmp * b_del_tmp;
      g_tmp = c_g_tmp_tmp / vJ_tmp * (b_g_tmp_tmp - 1.0);
      g[1] = g_tmp;
      b_g_tmp = sin(g_tmp_tmp);
      g[2] = -del_tmp / f_del_tmp * b_g_tmp;
      g[4] = g_tmp;
      g[5] = d_del_tmp / vJ_tmp * (b_g_tmp_tmp - 1.0) + 1.0;
      g[6] = -b_del_tmp / f_del_tmp * b_g_tmp;
      g[8] = del_tmp / f_del_tmp * b_g_tmp;
      g[9] = b_del_tmp / f_del_tmp * b_g_tmp;
      g[10] = b_g_tmp_tmp;
      g[12] = b_a * (del_tmp * (1.0 - b_g_tmp_tmp));
      g[13] = b_a * (b_del_tmp * (1.0 - b_g_tmp_tmp));
      g_tmp = f_del_tmp * b_g_tmp;
      g[14] = b_a * g_tmp;
      g[3] = 0.0;
      g[7] = 0.0;
      g[11] = 0.0;
      g[15] = 1.0;
      Smod_tmp_tmp = c_del_tmp * b_g_tmp;
      Smod_tmp = Smod_tmp_tmp * b_a_tmp;
      b_Smod_tmp = e_del_tmp * e_del_tmp;
      c_Smod_tmp = rt_powd_snf(e_del_tmp, 1.5);
      vJ_tmp = a_tmp * (b_g_tmp_tmp - 1.0);
      d_Smod_tmp = vJ_tmp / e_del_tmp;
      b_Smod_tmp_tmp = r * b_g_tmp;
      e_Smod_tmp = b_g_tmp_tmp * f_del_tmp - b_Smod_tmp_tmp;
      c_Smod_tmp_tmp = rt_powd_snf(e_del_tmp, 3.0);
      Smod_tmp_tmp_tmp = c_del_tmp * d_del_tmp;
      d_Smod_tmp_tmp = Smod_tmp_tmp_tmp * b_a_tmp * (b_g_tmp_tmp - 1.0);
      g_tmp_tmp = (g_tmp - 2.0 * r) + 2.0 * r * b_g_tmp_tmp;
      f_Smod_tmp = d_Smod_tmp_tmp * g_tmp_tmp / c_Smod_tmp_tmp;
      e_Smod_tmp_tmp = c_del_tmp * b_g_tmp_tmp;
      g_Smod_tmp = d_del_tmp + e_Smod_tmp_tmp;
      h_Smod_tmp =
          (Smod_tmp / c_Smod_tmp - d_Smod_tmp) +
          2.0 * r * c_del_tmp * b_a_tmp * (b_g_tmp_tmp - 1.0) / b_Smod_tmp;
      Smod[18 * i] = (g_Smod_tmp * h_Smod_tmp / e_del_tmp -
                      Smod_tmp * e_Smod_tmp / b_Smod_tmp) +
                     f_Smod_tmp;
      Smod[18 * i + 6] = 0.0;
      Smod_tmp = r * del_tmp * (b_g_tmp_tmp - 1.0);
      Smod[18 * i + 12] = Smod_tmp / e_del_tmp;
      Smod[18 * i + 1] = 0.0;
      f_Smod_tmp_tmp = d_del_tmp * b_g_tmp;
      i_Smod_tmp = f_Smod_tmp_tmp * b_a_tmp;
      g_Smod_tmp_tmp = d_del_tmp * b_g_tmp_tmp;
      j_Smod_tmp = c_del_tmp + g_Smod_tmp_tmp;
      d_Smod_tmp =
          (i_Smod_tmp / c_Smod_tmp - d_Smod_tmp) +
          2.0 * r * d_del_tmp * b_a_tmp * (b_g_tmp_tmp - 1.0) / b_Smod_tmp;
      Smod[18 * i + 7] = (j_Smod_tmp * d_Smod_tmp / e_del_tmp -
                          i_Smod_tmp * e_Smod_tmp / b_Smod_tmp) +
                         f_Smod_tmp;
      f_Smod_tmp = r * b_del_tmp * (b_g_tmp_tmp - 1.0);
      Smod[18 * i + 13] = f_Smod_tmp / e_del_tmp;
      i_Smod_tmp = c_Smod_tmp - b_Smod_tmp_tmp * e_del_tmp;
      k_Smod_tmp = rt_powd_snf(e_del_tmp, 2.5);
      l_Smod_tmp = del_tmp * b_a_tmp;
      Smod[18 * i + 2] = l_Smod_tmp * i_Smod_tmp / k_Smod_tmp;
      m_Smod_tmp = b_del_tmp * b_a_tmp;
      Smod[18 * i + 8] = m_Smod_tmp * i_Smod_tmp / k_Smod_tmp;
      Smod[18 * i + 14] = b_Smod_tmp_tmp / f_del_tmp;
      n_Smod_tmp = r * c_Smod_tmp;
      h_Smod_tmp_tmp = f_del_tmp - b_Smod_tmp_tmp;
      o_Smod_tmp = c_g_tmp_tmp * h_Smod_tmp_tmp;
      Smod[18 * i + 3] = -o_Smod_tmp / n_Smod_tmp;
      p_Smod_tmp = d_del_tmp * f_del_tmp + r * c_del_tmp * b_g_tmp;
      Smod[18 * i + 9] = -p_Smod_tmp / n_Smod_tmp;
      Smod[18 * i + 15] = 0.0;
      q_Smod_tmp = c_del_tmp * f_del_tmp + r * d_del_tmp * b_g_tmp;
      Smod[18 * i + 4] = q_Smod_tmp / n_Smod_tmp;
      Smod[18 * i + 10] = o_Smod_tmp / n_Smod_tmp;
      Smod[18 * i + 16] = 0.0;
      o_Smod_tmp = b_del_tmp * (b_g_tmp_tmp - 1.0);
      Smod[18 * i + 5] = -o_Smod_tmp / e_del_tmp;
      r_Smod_tmp = del_tmp * (b_g_tmp_tmp - 1.0);
      Smod[18 * i + 11] = r_Smod_tmp / e_del_tmp;
      Smod[18 * i + 17] = 0.0;
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
      fh3_tmp = qd[qi_data[2] - 1];
      fh3_tmp_tmp_tmp = qd[qi_data[0] - 1];
      fh3_tmp_tmp = 2.0 * del_tmp * fh3_tmp_tmp_tmp;
      b_fh3_tmp_tmp_tmp = qd[qi_data[1] - 1];
      b_fh3_tmp_tmp = 2.0 * b_del_tmp * b_fh3_tmp_tmp_tmp;
      b_fh3_tmp = fh3_tmp_tmp + b_fh3_tmp_tmp;
      c_fh3_tmp_tmp = fh3_tmp * c_del_tmp;
      c_fh3_tmp = c_fh3_tmp_tmp * b_g_tmp;
      e_fh3_tmp = 2.0 * r * b_Smod_tmp;
      f_fh3_tmp = fh3_tmp_tmp * b_g_tmp * b_a_tmp;
      g_fh3_tmp = e_Smod_tmp_tmp * b_a_tmp;
      h_fh3_tmp = 2.0 * c_del_tmp;
      i_fh3_tmp = 2.0 * c_Smod_tmp;
      j_fh3_tmp = b_g_tmp * b_a_tmp * b_fh3_tmp / i_fh3_tmp -
                  r * fh3_tmp * (b_g_tmp_tmp - 1.0) / e_del_tmp;
      k_fh3_tmp = 2.0 * k_Smod_tmp;
      l_fh3_tmp = 2.0 * r * fh3_tmp;
      a_tmp = vJ_tmp * b_fh3_tmp / b_Smod_tmp;
      g_tmp = 2.0 * r * f_del_tmp;
      m_fh3_tmp = b_g_tmp * b_g_tmp;
      d_fh3_tmp_tmp = b_g_tmp_tmp * b_fh3_tmp;
      e_fh3_tmp_tmp = 2.0 * f_del_tmp;
      d_fh3_tmp =
          d_Smod_tmp_tmp *
          (b_g_tmp * b_fh3_tmp / e_fh3_tmp_tmp - d_fh3_tmp_tmp / (2.0 * r)) /
          c_Smod_tmp_tmp;
      n_fh3_tmp = c_fh3_tmp_tmp * d_del_tmp * (b_g_tmp_tmp - 1.0) * g_tmp_tmp /
                  c_Smod_tmp_tmp;
      d_Smod_tmp_tmp = 2.0 * r * k_Smod_tmp;
      o_fh3_tmp = 3.0 * c_del_tmp * d_del_tmp * b_a_tmp * (b_g_tmp_tmp - 1.0) *
                  b_fh3_tmp * g_tmp_tmp / rt_powd_snf(e_del_tmp, 4.0);
      p_fh3_tmp = fh3_tmp_tmp * d_del_tmp * b_a_tmp * (b_g_tmp_tmp - 1.0) *
                  g_tmp_tmp / c_Smod_tmp_tmp;
      b_a = h_fh3_tmp * b_del_tmp * b_fh3_tmp_tmp_tmp * b_a_tmp *
            (b_g_tmp_tmp - 1.0) * g_tmp_tmp / c_Smod_tmp_tmp;
      c_fh3_tmp_tmp = rt_powd_snf(e_del_tmp, 3.5);
      vJ_tmp = Smod_tmp_tmp_tmp * b_g_tmp * b_a_tmp * b_fh3_tmp * g_tmp_tmp /
               (2.0 * r * c_fh3_tmp_tmp);
      fh3[0] =
          ((((((((((((g_Smod_tmp *
                          ((((((((j_fh3_tmp + c_fh3_tmp / c_Smod_tmp) -
                                 5.0 * c_del_tmp * b_g_tmp * b_a_tmp *
                                     b_fh3_tmp / k_fh3_tmp) +
                                f_fh3_tmp / c_Smod_tmp) +
                               l_fh3_tmp * c_del_tmp * (b_g_tmp_tmp - 1.0) /
                                   b_Smod_tmp) +
                              a_tmp) +
                             4.0 * r * del_tmp * fh3_tmp_tmp_tmp * b_a_tmp *
                                 (b_g_tmp_tmp - 1.0) / b_Smod_tmp) -
                            4.0 * r * c_del_tmp * b_a_tmp *
                                (b_g_tmp_tmp - 1.0) * b_fh3_tmp /
                                c_Smod_tmp_tmp) +
                           g_fh3_tmp * b_fh3_tmp / e_fh3_tmp) /
                          e_del_tmp +
                      h_Smod_tmp *
                          ((b_fh3_tmp_tmp + fh3_tmp_tmp * b_g_tmp_tmp) -
                           Smod_tmp_tmp * b_fh3_tmp / g_tmp) /
                          e_del_tmp) -
                     g_Smod_tmp * b_fh3_tmp * h_Smod_tmp / b_Smod_tmp) -
                    c_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
                   c_del_tmp * m_fh3_tmp * b_a_tmp * b_fh3_tmp / e_fh3_tmp) +
                  h_fh3_tmp * b_g_tmp * b_a_tmp * e_Smod_tmp * b_fh3_tmp /
                      c_Smod_tmp_tmp) -
                 d_fh3_tmp) -
                f_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
               n_fh3_tmp) -
              g_fh3_tmp * e_Smod_tmp * b_fh3_tmp / d_Smod_tmp_tmp) -
             o_fh3_tmp) +
            p_fh3_tmp) +
           b_a) -
          vJ_tmp;
      fh3[6] = 0.0;
      c_fh3_tmp = del_tmp * b_g_tmp * b_fh3_tmp;
      fh3[12] = (r * fh3_tmp_tmp_tmp * (b_g_tmp_tmp - 1.0) / e_del_tmp -
                 c_fh3_tmp / i_fh3_tmp) -
                Smod_tmp * b_fh3_tmp / b_Smod_tmp;
      fh3[1] = 0.0;
      f_fh3_tmp = fh3_tmp * d_del_tmp * b_g_tmp;
      g_fh3_tmp = b_fh3_tmp_tmp * b_g_tmp * b_a_tmp;
      h_fh3_tmp = g_Smod_tmp_tmp * b_a_tmp;
      fh3[7] =
          ((((((((((((j_Smod_tmp *
                          ((((((((j_fh3_tmp + f_fh3_tmp / c_Smod_tmp) -
                                 5.0 * d_del_tmp * b_g_tmp * b_a_tmp *
                                     b_fh3_tmp / k_fh3_tmp) +
                                g_fh3_tmp / c_Smod_tmp) +
                               l_fh3_tmp * d_del_tmp * (b_g_tmp_tmp - 1.0) /
                                   b_Smod_tmp) +
                              a_tmp) +
                             4.0 * r * b_del_tmp * b_fh3_tmp_tmp_tmp * b_a_tmp *
                                 (b_g_tmp_tmp - 1.0) / b_Smod_tmp) -
                            4.0 * r * d_del_tmp * b_a_tmp *
                                (b_g_tmp_tmp - 1.0) * b_fh3_tmp /
                                c_Smod_tmp_tmp) +
                           h_fh3_tmp * b_fh3_tmp / e_fh3_tmp) /
                          e_del_tmp +
                      d_Smod_tmp *
                          ((fh3_tmp_tmp + b_fh3_tmp_tmp * b_g_tmp_tmp) -
                           f_Smod_tmp_tmp * b_fh3_tmp / g_tmp) /
                          e_del_tmp) -
                     j_Smod_tmp * b_fh3_tmp * d_Smod_tmp / b_Smod_tmp) -
                    f_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
                   d_del_tmp * m_fh3_tmp * b_a_tmp * b_fh3_tmp / e_fh3_tmp) +
                  2.0 * d_del_tmp * b_g_tmp * b_a_tmp * e_Smod_tmp * b_fh3_tmp /
                      c_Smod_tmp_tmp) -
                 d_fh3_tmp) -
                g_fh3_tmp * e_Smod_tmp / b_Smod_tmp) +
               n_fh3_tmp) -
              h_fh3_tmp * e_Smod_tmp * b_fh3_tmp / d_Smod_tmp_tmp) -
             o_fh3_tmp) +
            p_fh3_tmp) +
           b_a) -
          vJ_tmp;
      e_fh3_tmp = b_del_tmp * b_g_tmp * b_fh3_tmp;
      fh3[13] = (r * b_fh3_tmp_tmp_tmp * (b_g_tmp_tmp - 1.0) / e_del_tmp -
                 e_fh3_tmp / i_fh3_tmp) -
                f_Smod_tmp * b_fh3_tmp / b_Smod_tmp;
      vJ_tmp = b_Smod_tmp_tmp * b_fh3_tmp;
      f_fh3_tmp = (d_fh3_tmp_tmp * f_del_tmp / 2.0 -
                   3.0 * b_fh3_tmp * f_del_tmp / 2.0) +
                  vJ_tmp;
      g_fh3_tmp = 2.0 * c_fh3_tmp_tmp;
      fh3[2] = ((fh3_tmp_tmp_tmp * b_a_tmp * i_Smod_tmp / k_Smod_tmp -
                 l_Smod_tmp * f_fh3_tmp / k_Smod_tmp) +
                del_tmp * fh3_tmp * i_Smod_tmp / k_Smod_tmp) -
               5.0 * del_tmp * b_a_tmp * i_Smod_tmp * b_fh3_tmp / g_fh3_tmp;
      fh3[8] = ((b_fh3_tmp_tmp_tmp * b_a_tmp * i_Smod_tmp / k_Smod_tmp -
                 m_Smod_tmp * f_fh3_tmp / k_Smod_tmp) +
                b_del_tmp * fh3_tmp * i_Smod_tmp / k_Smod_tmp) -
               5.0 * b_del_tmp * b_a_tmp * i_Smod_tmp * b_fh3_tmp / g_fh3_tmp;
      fh3[14] = d_fh3_tmp_tmp / (2.0 * e_del_tmp) - vJ_tmp / i_fh3_tmp;
      fh3_tmp = c_g_tmp_tmp *
                (b_fh3_tmp / e_fh3_tmp_tmp - d_fh3_tmp_tmp / e_fh3_tmp_tmp) /
                n_Smod_tmp;
      f_fh3_tmp = del_tmp * b_fh3_tmp_tmp_tmp * h_Smod_tmp_tmp / n_Smod_tmp;
      g_fh3_tmp = fh3_tmp_tmp_tmp * b_del_tmp * h_Smod_tmp_tmp / n_Smod_tmp;
      h_fh3_tmp = 3.0 * del_tmp * b_del_tmp * h_Smod_tmp_tmp * b_fh3_tmp /
                  d_Smod_tmp_tmp;
      fh3[3] = ((h_fh3_tmp - f_fh3_tmp) - g_fh3_tmp) - fh3_tmp;
      fh3[9] = 3.0 * p_Smod_tmp * b_fh3_tmp / d_Smod_tmp_tmp -
               (((d_del_tmp * b_fh3_tmp / e_fh3_tmp_tmp +
                  b_fh3_tmp_tmp * f_del_tmp) +
                 2.0 * r * del_tmp * fh3_tmp_tmp_tmp * b_g_tmp) +
                e_Smod_tmp_tmp * b_fh3_tmp / e_fh3_tmp_tmp) /
                   n_Smod_tmp;
      fh3[15] = 0.0;
      fh3[4] =
          (((c_del_tmp * b_fh3_tmp / e_fh3_tmp_tmp + fh3_tmp_tmp * f_del_tmp) +
            2.0 * r * b_del_tmp * b_fh3_tmp_tmp_tmp * b_g_tmp) +
           g_Smod_tmp_tmp * b_fh3_tmp / e_fh3_tmp_tmp) /
              n_Smod_tmp -
          3.0 * q_Smod_tmp * b_fh3_tmp / d_Smod_tmp_tmp;
      fh3[10] = ((fh3_tmp + f_fh3_tmp) + g_fh3_tmp) - h_fh3_tmp;
      fh3[16] = 0.0;
      fh3_tmp = 2.0 * r * c_Smod_tmp;
      fh3[5] = (o_Smod_tmp * b_fh3_tmp / b_Smod_tmp -
                b_fh3_tmp_tmp_tmp * (b_g_tmp_tmp - 1.0) / e_del_tmp) +
               e_fh3_tmp / fh3_tmp;
      fh3[11] = (fh3_tmp_tmp_tmp * (b_g_tmp_tmp - 1.0) / e_del_tmp -
                 r_Smod_tmp * b_fh3_tmp / b_Smod_tmp) -
                c_fh3_tmp / fh3_tmp;
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
    dv[3] = -g[14];
    dv[6] = g[13];
    dv[1] = g[14];
    dv[7] = -g[12];
    dv[2] = -g[13];
    dv[5] = g[12];
    for (b_i = 0; b_i < 3; b_i++) {
      b_a = b_skw[b_i];
      g_tmp_tmp = b_skw[b_i + 3];
      a_tmp = b_skw[b_i + 6];
      for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
        c_skw[b_i + 3 * b_XJ_tmp] =
            (b_a * dv[3 * b_XJ_tmp] + g_tmp_tmp * dv[3 * b_XJ_tmp + 1]) +
            a_tmp * dv[3 * b_XJ_tmp + 2];
        XJ[(b_XJ_tmp + 6 * b_i) + 36 * (i + 1)] = skw[b_XJ_tmp + 3 * b_i];
      }
    }
    XJ_tmp = 36 * (i + 1);
    for (b_i = 0; b_i < 3; b_i++) {
      b_XJ_tmp = 6 * (b_i + 3) + XJ_tmp;
      XJ[b_XJ_tmp] = c_skw[3 * b_i];
      c_XJ_tmp = 6 * b_i + XJ_tmp;
      XJ[c_XJ_tmp + 3] = 0.0;
      XJ[b_XJ_tmp + 3] = skw[3 * b_i];
      j = 3 * b_i + 1;
      XJ[b_XJ_tmp + 1] = c_skw[j];
      XJ[c_XJ_tmp + 4] = 0.0;
      XJ[b_XJ_tmp + 4] = skw[j];
      j = 3 * b_i + 2;
      XJ[b_XJ_tmp + 2] = c_skw[j];
      XJ[c_XJ_tmp + 5] = 0.0;
      XJ[b_XJ_tmp + 5] = skw[j];
    }
    for (b_i = 0; b_i < 6; b_i++) {
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        b_a = 0.0;
        for (i1 = 0; i1 < 6; i1++) {
          b_a += XJ[(b_i + 6 * i1) + 36 * (i + 1)] *
                 (double)Xtree[i1 + 6 * b_XJ_tmp];
        }
        Xup[(b_i + 6 * b_XJ_tmp) + 36 * (i + 1)] = b_a;
      }
    }
    for (b_i = 0; b_i < 36; b_i++) {
      Xtree[b_i] = 0;
    }
    g_tmp = qd[qi_data[0] - 1];
    d_Smod_tmp_tmp = qd[qi_data[1] - 1];
    d_fh3_tmp = qd[qi_data[2] - 1];
    for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
      Xtree[XJ_tmp + 6 * XJ_tmp] = 1;
      b_i = XJ_tmp + 18 * i;
      b_a = (Smod[b_i] * g_tmp + Smod[b_i + 6] * d_Smod_tmp_tmp) +
            Smod[b_i + 12] * d_fh3_tmp;
      vJ[XJ_tmp] = b_a;
      g_tmp_tmp = 0.0;
      for (b_i = 0; b_i < 6; b_i++) {
        g_tmp_tmp += Xup[(XJ_tmp + 6 * b_i) + 36 * (i + 1)] * b_v[b_i + 6 * i];
      }
      b_Xup[XJ_tmp] = g_tmp_tmp + b_a;
    }
    for (b_i = 0; b_i < 6; b_i++) {
      b_v[b_i + 6 * (i + 1)] = b_Xup[b_i];
    }
    memset(&I_tmp[0], 0, 36U * sizeof(double));
    skw[0] = 0.0;
    b_i = 6 * (i + 1);
    b_a = b_v[b_i + 5];
    skw[3] = -b_a;
    g_tmp_tmp = b_v[b_i + 4];
    skw[6] = g_tmp_tmp;
    skw[1] = b_a;
    skw[4] = 0.0;
    a_tmp = b_v[b_i + 3];
    skw[7] = -a_tmp;
    skw[2] = -g_tmp_tmp;
    skw[5] = a_tmp;
    skw[8] = 0.0;
    for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
      vJ_tmp = skw[3 * b_XJ_tmp];
      I_tmp[6 * b_XJ_tmp] = vJ_tmp;
      XJ_tmp = 6 * (b_XJ_tmp + 3);
      I_tmp[XJ_tmp + 3] = vJ_tmp;
      vJ_tmp = skw[3 * b_XJ_tmp + 1];
      I_tmp[6 * b_XJ_tmp + 1] = vJ_tmp;
      I_tmp[XJ_tmp + 4] = vJ_tmp;
      vJ_tmp = skw[3 * b_XJ_tmp + 2];
      I_tmp[6 * b_XJ_tmp + 2] = vJ_tmp;
      I_tmp[XJ_tmp + 5] = vJ_tmp;
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
    for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
      vJ_tmp = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        vJ_tmp += Xup[(b_XJ_tmp + 6 * i1) + 36 * (i + 1)] * a[i1 + 6 * i];
      }
      b_Xup[b_XJ_tmp] = vJ_tmp;
      c_I[b_XJ_tmp] =
          (fh3[b_XJ_tmp] * g_tmp + fh3[b_XJ_tmp + 6] * d_Smod_tmp_tmp) +
          fh3[b_XJ_tmp + 12] * d_fh3_tmp;
    }
    for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
      vJ_tmp = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        vJ_tmp += I_tmp[b_XJ_tmp + 6 * i1] * vJ[i1];
      }
      a[b_XJ_tmp + b_i] = (b_Xup[b_XJ_tmp] + c_I[b_XJ_tmp]) + vJ_tmp;
    }
    memset(&I_tmp[0], 0, 36U * sizeof(double));
    skw[0] = 0.0;
    skw[3] = -b_a;
    skw[6] = g_tmp_tmp;
    skw[1] = b_a;
    skw[4] = 0.0;
    skw[7] = -a_tmp;
    skw[2] = -g_tmp_tmp;
    skw[5] = a_tmp;
    skw[8] = 0.0;
    for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
      vJ_tmp = skw[3 * b_XJ_tmp];
      I_tmp[6 * b_XJ_tmp] = vJ_tmp;
      XJ_tmp = 6 * (b_XJ_tmp + 3);
      I_tmp[XJ_tmp + 3] = vJ_tmp;
      vJ_tmp = skw[3 * b_XJ_tmp + 1];
      I_tmp[6 * b_XJ_tmp + 1] = vJ_tmp;
      I_tmp[XJ_tmp + 4] = vJ_tmp;
      vJ_tmp = skw[3 * b_XJ_tmp + 2];
      I_tmp[6 * b_XJ_tmp + 2] = vJ_tmp;
      I_tmp[XJ_tmp + 5] = vJ_tmp;
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
    for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
      for (i1 = 0; i1 < 6; i1++) {
        b_I_tmp[i1 + 6 * b_XJ_tmp] = -I_tmp[b_XJ_tmp + 6 * i1];
      }
    }
    for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
      c_I[b_XJ_tmp] = 0.0;
      vJ[b_XJ_tmp] = 0.0;
      for (i1 = 0; i1 < 6; i1++) {
        b_a = 0.0;
        for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
          b_a += b_I_tmp[b_XJ_tmp + 6 * XJ_tmp] *
                 b_I[(XJ_tmp + 6 * i1) + 36 * (i + 1)];
        }
        XJ_tmp = i1 + b_i;
        c_I[b_XJ_tmp] += b_I[(b_XJ_tmp + 6 * i1) + 36 * (i + 1)] * a[XJ_tmp];
        vJ[b_XJ_tmp] += b_a * b_v[XJ_tmp];
      }
      f[b_XJ_tmp + b_i] = c_I[b_XJ_tmp] + vJ[b_XJ_tmp];
    }
  }
  /*  Composite Rigid Body Algorithm to calculate M */
  /*  composite inertia calculation */
  for (i = 0; i < 4; i++) {
    if (4 - i == 1) {
      b_a = 0.0;
      for (b_i = 0; b_i < 6; b_i++) {
        b_a += (double)d_a[b_i] * f[b_i];
      }
      C[0] = b_a;
    } else {
      XJ_tmp = (2 - i) * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        qi_data[b_i] = (signed char)((XJ_tmp + b_i) + 2);
      }
      for (b_i = 0; b_i < 3; b_i++) {
        tmp_data[b_i] = (signed char)(qi_data[b_i] - 1);
      }
      for (b_i = 0; b_i < 3; b_i++) {
        b_a = 0.0;
        for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
          b_a += 2.0 * Smod[(b_XJ_tmp + 6 * b_i) + 18 * (2 - i)] *
                 f[b_XJ_tmp + 6 * (3 - i)];
        }
        b_qd[b_i] = b_a;
      }
      for (b_i = 0; b_i < 3; b_i++) {
        C[tmp_data[b_i]] = b_qd[b_i];
      }
      for (b_i = 0; b_i < 6; b_i++) {
        b_a = 0.0;
        for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
          b_a += Xup[(b_XJ_tmp + 6 * b_i) + 36 * (3 - i)] *
                 f[b_XJ_tmp + 6 * (3 - i)];
        }
        c_I[b_i] = f[b_i + 6 * (2 - i)] + b_a;
      }
      for (b_i = 0; b_i < 6; b_i++) {
        f[b_i + 6 * (2 - i)] = c_I[b_i];
      }
    }
    if (4 - i != 1) {
      b_i = 36 * (3 - i);
      for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
        for (i1 = 0; i1 < 6; i1++) {
          b_a = 0.0;
          for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
            b_a += Xup[(XJ_tmp + 6 * b_XJ_tmp) + b_i] *
                   b_I[(XJ_tmp + 6 * i1) + b_i];
          }
          I_tmp[b_XJ_tmp + 6 * i1] = b_a;
        }
        for (i1 = 0; i1 < 6; i1++) {
          b_a = 0.0;
          for (XJ_tmp = 0; XJ_tmp < 6; XJ_tmp++) {
            b_a += I_tmp[b_XJ_tmp + 6 * XJ_tmp] *
                   Xup[(XJ_tmp + 6 * i1) + 36 * (3 - i)];
          }
          XJ_tmp = b_XJ_tmp + 6 * i1;
          b_I_tmp[XJ_tmp] = b_I[XJ_tmp + 36 * (2 - i)] + b_a;
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
      b_a = 0.0;
      for (b_i = 0; b_i < 6; b_i++) {
        g_tmp_tmp = 0.0;
        for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
          g_tmp_tmp += b_I[b_i + 6 * b_XJ_tmp] * (double)c_a[b_XJ_tmp];
        }
        b_a += (double)iv[b_i] * g_tmp_tmp;
      }
      M[0] = b_a;
    } else {
      signed char b_tmp_data[9];
      XJ_tmp = (i - 1) * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        qi_data[b_i] = (signed char)((XJ_tmp + b_i) + 2);
      }
      for (b_i = 0; b_i < 6; b_i++) {
        for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
          b_a = 0.0;
          for (i1 = 0; i1 < 6; i1++) {
            b_a += b_I[(b_i + 6 * i1) + 36 * i] *
                   Smod[(i1 + 6 * b_XJ_tmp) + 18 * (i - 1)];
          }
          fh3[b_i + 6 * b_XJ_tmp] = b_a;
        }
      }
      for (b_i = 0; b_i < 3; b_i++) {
        b_XJ_tmp = qi_data[b_i];
        tmp_data[b_i] = (signed char)(b_XJ_tmp - 1);
        b_tmp_data[b_i] = (signed char)(b_XJ_tmp - 1);
      }
      for (b_i = 0; b_i < 3; b_i++) {
        for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
          b_a = 0.0;
          for (i1 = 0; i1 < 6; i1++) {
            b_a += Smod[(i1 + 6 * b_i) + 18 * (i - 1)] * fh3[i1 + 6 * b_XJ_tmp];
          }
          skw[b_i + 3 * b_XJ_tmp] = b_a;
        }
      }
      for (b_i = 0; b_i < 3; b_i++) {
        for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
          M[tmp_data[b_XJ_tmp] + 10 * b_tmp_data[b_i]] =
              skw[b_XJ_tmp + 3 * b_i];
        }
      }
      j = i - 1;
      while (j + 2 > 1) {
        for (b_i = 0; b_i < 6; b_i++) {
          for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
            b_a = 0.0;
            for (i1 = 0; i1 < 6; i1++) {
              b_a +=
                  Xup[(i1 + 6 * b_i) + 36 * (j + 1)] * fh3[i1 + 6 * b_XJ_tmp];
            }
            c_Xup[b_i + 6 * b_XJ_tmp] = b_a;
          }
        }
        memcpy(&fh3[0], &c_Xup[0], 18U * sizeof(double));
        j--;
        if (j + 2 == 1) {
          for (b_i = 0; b_i < 3; b_i++) {
            b_a = 0.0;
            for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
              b_a += fh3[b_XJ_tmp + 6 * b_i] * (double)c_a[b_XJ_tmp];
            }
            b_qd[b_i] = b_a;
          }
          for (b_i = 0; b_i < 3; b_i++) {
            M[b_tmp_data[b_i]] = b_qd[b_i];
          }
          /* (S{j}' * fh).'; */
          for (b_i = 0; b_i < 3; b_i++) {
            b_a = 0.0;
            for (b_XJ_tmp = 0; b_XJ_tmp < 6; b_XJ_tmp++) {
              b_a += (double)iv[b_XJ_tmp] * fh3[b_XJ_tmp + 6 * b_i];
            }
            b_qd[b_i] = b_a;
          }
          for (b_i = 0; b_i < 3; b_i++) {
            M[10 * b_tmp_data[b_i]] = b_qd[b_i];
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
            b_XJ_tmp = qi_data[b_i];
            c_tmp_data[b_i] = (signed char)b_XJ_tmp;
            d_tmp_data[b_i] = (signed char)(b_XJ_tmp - 1);
          }
          for (b_i = 0; b_i < 3; b_i++) {
            for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
              b_a = 0.0;
              for (i1 = 0; i1 < 6; i1++) {
                b_a += fh3[i1 + 6 * b_i] * Smod[(i1 + 6 * b_XJ_tmp) + 18 * j];
              }
              skw[b_i + 3 * b_XJ_tmp] = b_a;
            }
          }
          /* (S{j}' * fh).'; */
          for (b_i = 0; b_i < c_XJ_tmp; b_i++) {
            for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
              M[b_tmp_data[b_XJ_tmp] + 10 * d_tmp_data[b_i]] =
                  skw[b_XJ_tmp + 3 * b_i];
            }
            d_tmp_data[b_i] = (signed char)(c_tmp_data[b_i] - 1);
          }
          for (b_i = 0; b_i < 3; b_i++) {
            for (b_XJ_tmp = 0; b_XJ_tmp < 3; b_XJ_tmp++) {
              b_a = 0.0;
              for (i1 = 0; i1 < 6; i1++) {
                b_a += Smod[(i1 + 6 * b_i) + 18 * j] * fh3[i1 + 6 * b_XJ_tmp];
              }
              skw[b_i + 3 * b_XJ_tmp] = b_a;
            }
          }
          for (b_i = 0; b_i < 3; b_i++) {
            for (b_XJ_tmp = 0; b_XJ_tmp < c_XJ_tmp; b_XJ_tmp++) {
              M[d_tmp_data[b_XJ_tmp] + 10 * b_tmp_data[b_i]] =
                  skw[b_XJ_tmp + c_XJ_tmp * b_i];
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
}

/*
 * File trailer for MC_3_cg.c
 *
 * [EOF]
 */
