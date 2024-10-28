/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: MC_3_cg.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 17:40:15
 */

/* Include Files */
#include "MC_3_cg.h"
#include "PCC_jacobian.h"
#include "helix_controller_emxutil.h"
#include "helix_controller_types.h"
#include "rt_nonfinite.h"
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
 *                double N
 *                emxArray_real_T *M
 *                emxArray_real_T *C
 *                emxArray_real_T *J
 *                double X[36]
 * Return Type  : void
 */
void MC_3_cg(const double q[10], const double qd[10], double m, double r,
             double L0, double N, emxArray_real_T *M, emxArray_real_T *C,
             emxArray_real_T *J, double X[36])
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
  static const signed char b_a[6] = {0, 0, 0, 0, 0, 1};
  static const signed char c_a[6] = {0, 0, 0, 0, 0, 2};
  static const signed char iv1[6] = {0, 0, 0, 0, 0, 1};
  emxArray_int32_T *b_r;
  emxArray_int32_T *r1;
  emxArray_real_T *Smod;
  emxArray_real_T *XJ;
  emxArray_real_T *Xup;
  emxArray_real_T *a;
  emxArray_real_T *b_I;
  emxArray_real_T *b_q;
  emxArray_real_T *b_qd;
  emxArray_real_T *b_v;
  emxArray_real_T *f;
  emxArray_real_T *qi;
  double I_tmp[36];
  double b_I_tmp[36];
  double c_Xup[18];
  double fh3[18];
  double g[16];
  double b_skw[9];
  double c_skw[9];
  double skw[9];
  double c_I[6];
  double vJ[6];
  double c_qd[3];
  double Nq_tmp;
  double d;
  double d1;
  double vJ_tmp;
  double vJ_tmp_tmp;
  double *C_data;
  double *I_data;
  double *J_data;
  double *Smod_data;
  double *XJ_data;
  double *Xup_data;
  double *a_data;
  double *f_data;
  double *q_data;
  double *qd_data;
  double *qi_data;
  double *v_data;
  int b_i;
  int i;
  int i1;
  int i2;
  int i3;
  int j;
  int loop_ub_tmp;
  int *r2;
  int *r3;
  signed char Xtree[36];
  signed char iv[9];
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  This function calculates the Mass and Coriolis + Gravity Matrices */
  /*  given the current state and geometric parameters of the pushpuppet robot
   */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  qi = {1, 2:4, 5, 6:8, 9, 10:12, 13, 14:16}; */
  Nq_tmp = (N - 1.0) * 3.0 + 1.0;
  emxInit_real_T(&b_I, 3);
  i = b_I->size[0] * b_I->size[1] * b_I->size[2];
  b_I->size[0] = 6;
  b_I->size[1] = 6;
  i1 = (int)N;
  b_I->size[2] = (int)N;
  emxEnsureCapacity_real_T(b_I, i);
  I_data = b_I->data;
  loop_ub_tmp = 36 * (int)N;
  for (i = 0; i < loop_ub_tmp; i++) {
    I_data[i] = 0.0;
  }
  for (i = 0; i < 36; i++) {
    I_data[i] = 0.0;
  }
  for (j = 0; j < 6; j++) {
    I_data[j + 6 * j] = v[j];
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
  for (i = 0; i < 36; i++) {
    I_data[i + 36] = I_tmp[i];
  }
  i = (int)(N - 2.0);
  for (b_i = 0; b_i < i; b_i++) {
    for (i2 = 0; i2 < 36; i2++) {
      I_data[(i2 + b_i * 36) + 72] = I_tmp[i2];
    }
  }
  emxInit_real_T(&Smod, 3);
  i = Smod->size[0] * Smod->size[1] * Smod->size[2];
  Smod->size[0] = 6;
  Smod->size[1] = 3;
  i2 = (int)(N - 1.0);
  Smod->size[2] = (int)(N - 1.0);
  emxEnsureCapacity_real_T(Smod, i);
  Smod_data = Smod->data;
  j = 18 * (int)(N - 1.0);
  for (i = 0; i < j; i++) {
    Smod_data[i] = 0.0;
  }
  emxInit_real_T(&XJ, 3);
  i = XJ->size[0] * XJ->size[1] * XJ->size[2];
  XJ->size[0] = 6;
  XJ->size[1] = 6;
  XJ->size[2] = (int)N;
  emxEnsureCapacity_real_T(XJ, i);
  XJ_data = XJ->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    XJ_data[i] = 0.0;
  }
  emxInit_real_T(&Xup, 3);
  i = Xup->size[0] * Xup->size[1] * Xup->size[2];
  Xup->size[0] = 6;
  Xup->size[1] = 6;
  Xup->size[2] = (int)N;
  emxEnsureCapacity_real_T(Xup, i);
  Xup_data = Xup->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    Xup_data[i] = 0.0;
  }
  emxInit_real_T(&b_v, 2);
  i = b_v->size[0] * b_v->size[1];
  b_v->size[0] = 6;
  b_v->size[1] = (int)N;
  emxEnsureCapacity_real_T(b_v, i);
  v_data = b_v->data;
  loop_ub_tmp = 6 * (int)N;
  for (i = 0; i < loop_ub_tmp; i++) {
    v_data[i] = 0.0;
  }
  emxInit_real_T(&a, 2);
  i = a->size[0] * a->size[1];
  a->size[0] = 6;
  a->size[1] = (int)N;
  emxEnsureCapacity_real_T(a, i);
  a_data = a->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    a_data[i] = 0.0;
  }
  emxInit_real_T(&f, 2);
  i = f->size[0] * f->size[1];
  f->size[0] = 6;
  f->size[1] = (int)N;
  emxEnsureCapacity_real_T(f, i);
  f_data = f->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    f_data[i] = 0.0;
  }
  loop_ub_tmp = (int)Nq_tmp;
  i = C->size[0];
  C->size[0] = (int)Nq_tmp;
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    C_data[i] = 0.0;
  }
  i = J->size[0] * J->size[1];
  J->size[0] = 6;
  J->size[1] = (int)Nq_tmp;
  emxEnsureCapacity_real_T(J, i);
  J_data = J->data;
  j = 6 * (int)Nq_tmp;
  for (i = 0; i < j; i++) {
    J_data[i] = 0.0;
  }
  for (i = 0; i < 36; i++) {
    Xtree[i] = b[i];
  }
  /* [diag(ones(3,1)) [0;0;0]; 0 0 0 1]; */
  vJ_tmp_tmp = sin(q[0]);
  vJ_tmp = cos(q[0]);
  g[0] = vJ_tmp;
  g[4] = -vJ_tmp_tmp;
  g[8] = 0.0;
  g[12] = 0.0;
  g[1] = vJ_tmp_tmp;
  g[5] = vJ_tmp;
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
  iv[0] = 0;
  iv[3] = 0;
  iv[6] = 0;
  iv[1] = 0;
  iv[4] = 0;
  iv[7] = 0;
  iv[2] = 0;
  iv[5] = 0;
  iv[8] = 0;
  for (i = 0; i < 3; i++) {
    d = b_skw[i];
    vJ_tmp = b_skw[i + 3];
    d1 = b_skw[i + 6];
    for (i3 = 0; i3 < 3; i3++) {
      c_skw[i + 3 * i3] =
          (d * (double)iv[3 * i3] + vJ_tmp * (double)iv[3 * i3 + 1]) +
          d1 * (double)iv[3 * i3 + 2];
      XJ_data[i3 + 6 * i] = skw[i3 + 3 * i];
    }
  }
  for (i = 0; i < 3; i++) {
    i3 = 6 * (i + 3);
    XJ_data[i3] = c_skw[3 * i];
    XJ_data[6 * i + 3] = 0.0;
    XJ_data[i3 + 3] = skw[3 * i];
    loop_ub_tmp = 3 * i + 1;
    XJ_data[i3 + 1] = c_skw[loop_ub_tmp];
    XJ_data[6 * i + 4] = 0.0;
    XJ_data[i3 + 4] = skw[loop_ub_tmp];
    loop_ub_tmp = 3 * i + 2;
    XJ_data[i3 + 2] = c_skw[loop_ub_tmp];
    XJ_data[6 * i + 5] = 0.0;
    XJ_data[i3 + 5] = skw[loop_ub_tmp];
  }
  for (b_i = 0; b_i < 6; b_i++) {
    for (i = 0; i < 6; i++) {
      d = 0.0;
      for (i3 = 0; i3 < 6; i3++) {
        d += XJ_data[b_i + 6 * i3] * (double)b[i3 + 6 * i];
      }
      Xup_data[b_i + 6 * i] = d;
    }
    d = (double)b_a[b_i] * qd[0];
    vJ[b_i] = d;
    v_data[b_i] = d;
  }
  memset(&I_tmp[0], 0, 36U * sizeof(double));
  skw[0] = 0.0;
  skw[3] = -v_data[5];
  skw[6] = v_data[4];
  skw[1] = v_data[5];
  skw[4] = 0.0;
  skw[7] = -v_data[3];
  skw[2] = -v_data[4];
  skw[5] = v_data[3];
  skw[8] = 0.0;
  for (i = 0; i < 3; i++) {
    vJ_tmp_tmp = skw[3 * i];
    I_tmp[6 * i] = vJ_tmp_tmp;
    loop_ub_tmp = 6 * (i + 3);
    I_tmp[loop_ub_tmp + 3] = vJ_tmp_tmp;
    vJ_tmp_tmp = skw[3 * i + 1];
    I_tmp[6 * i + 1] = vJ_tmp_tmp;
    I_tmp[loop_ub_tmp + 4] = vJ_tmp_tmp;
    vJ_tmp_tmp = skw[3 * i + 2];
    I_tmp[6 * i + 2] = vJ_tmp_tmp;
    I_tmp[loop_ub_tmp + 5] = vJ_tmp_tmp;
  }
  I_tmp[18] = 0.0;
  I_tmp[24] = -v_data[2];
  I_tmp[30] = v_data[1];
  I_tmp[19] = v_data[2];
  I_tmp[25] = 0.0;
  I_tmp[31] = -v_data[0];
  I_tmp[20] = -v_data[1];
  I_tmp[26] = v_data[0];
  I_tmp[32] = 0.0;
  for (i = 0; i < 6; i++) {
    d = 0.0;
    vJ_tmp = 0.0;
    for (i3 = 0; i3 < 6; i3++) {
      j = i + 6 * i3;
      d += Xup_data[j] * b_b[i3];
      d1 = I_tmp[j];
      vJ_tmp += d1 * vJ[i3];
      b_I_tmp[i3 + 6 * i] = -d1;
    }
    a_data[i] = d + vJ_tmp;
  }
  for (i = 0; i < 6; i++) {
    c_I[i] = 0.0;
    vJ[i] = 0.0;
    for (i3 = 0; i3 < 6; i3++) {
      d = 0.0;
      for (j = 0; j < 6; j++) {
        d += b_I_tmp[i + 6 * j] * I_data[j + 6 * i3];
      }
      c_I[i] += I_data[i + 6 * i3] * a_data[i3];
      vJ[i] += d * v_data[i3];
      j = i3 + 6 * i;
      X[j] = Xup_data[j];
    }
    f_data[i] = c_I[i] + vJ[i];
  }
  /*  Recursive Newton Euler to Calculate C+G */
  emxInit_real_T(&qi, 2);
  qi_data = qi->data;
  emxInit_real_T(&b_q, 1);
  emxInit_real_T(&b_qd, 1);
  for (b_i = 0; b_i < i2; b_i++) {
    double b_Xup[6];
    double b_qd_tmp;
    double c_qd_tmp;
    double d2;
    double d3;
    double d4;
    double qd_tmp;
    vJ_tmp_tmp = (((double)b_i + 2.0) - 2.0) * 3.0;
    if (vJ_tmp_tmp + 4.0 < vJ_tmp_tmp + 2.0) {
      qi->size[0] = 1;
      qi->size[1] = 0;
    } else {
      i = qi->size[0] * qi->size[1];
      qi->size[0] = 1;
      j = (int)((vJ_tmp_tmp + 4.0) - (vJ_tmp_tmp + 2.0));
      qi->size[1] = j + 1;
      emxEnsureCapacity_real_T(qi, i);
      qi_data = qi->data;
      for (i = 0; i <= j; i++) {
        qi_data[i] = (vJ_tmp_tmp + 2.0) + (double)i;
      }
    }
    i = b_q->size[0];
    b_q->size[0] = qi->size[1];
    emxEnsureCapacity_real_T(b_q, i);
    q_data = b_q->data;
    j = qi->size[1];
    i = b_qd->size[0];
    b_qd->size[0] = qi->size[1];
    emxEnsureCapacity_real_T(b_qd, i);
    qd_data = b_qd->data;
    for (i = 0; i < j; i++) {
      d = qi_data[i];
      q_data[i] = q[(int)d - 1];
      qd_data[i] = qd[(int)d - 1];
    }
    PCC_jacobian(b_q, r, L0, b_qd, I_tmp, c_Xup, g, fh3);
    memcpy(&(*(double(*)[18]) & Smod_data[18 * b_i])[0], &c_Xup[0],
           18U * sizeof(double));
    memcpy(&(*(double(*)[36]) & XJ_data[36 * (b_i + 1)])[0], &I_tmp[0],
           36U * sizeof(double));
    for (i = 0; i < 6; i++) {
      for (i3 = 0; i3 < 6; i3++) {
        d = 0.0;
        for (j = 0; j < 6; j++) {
          d +=
              XJ_data[(i + 6 * j) + 36 * (b_i + 1)] * (double)Xtree[j + 6 * i3];
        }
        Xup_data[(i + 6 * i3) + 36 * (b_i + 1)] = d;
      }
      for (i3 = 0; i3 < 6; i3++) {
        d = 0.0;
        for (j = 0; j < 6; j++) {
          d += Xup_data[(i + 6 * j) + 36 * (b_i + 1)] * X[j + 6 * i3];
        }
        I_tmp[i + 6 * i3] = d;
      }
    }
    for (i = 0; i < 36; i++) {
      X[i] = I_tmp[i];
      Xtree[i] = 0;
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < 6; loop_ub_tmp++) {
      Xtree[loop_ub_tmp + 6 * loop_ub_tmp] = 1;
    }
    qd_tmp = qd[(int)qi_data[0] - 1];
    b_qd_tmp = qd[(int)qi_data[1] - 1];
    c_qd_tmp = qd[(int)qi_data[2] - 1];
    for (i = 0; i < 6; i++) {
      i3 = i + 18 * b_i;
      d = (Smod_data[i3] * qd_tmp + Smod_data[i3 + 6] * b_qd_tmp) +
          Smod_data[i3 + 12] * c_qd_tmp;
      vJ[i] = d;
      vJ_tmp = 0.0;
      for (i3 = 0; i3 < 6; i3++) {
        vJ_tmp +=
            Xup_data[(i + 6 * i3) + 36 * (b_i + 1)] * v_data[i3 + 6 * b_i];
      }
      b_Xup[i] = vJ_tmp + d;
    }
    for (i = 0; i < 6; i++) {
      v_data[i + 6 * (b_i + 1)] = b_Xup[i];
    }
    memset(&I_tmp[0], 0, 36U * sizeof(double));
    skw[0] = 0.0;
    i = 6 * (b_i + 1);
    d = v_data[i + 5];
    skw[3] = -d;
    vJ_tmp = v_data[i + 4];
    skw[6] = vJ_tmp;
    skw[1] = d;
    skw[4] = 0.0;
    d1 = v_data[i + 3];
    skw[7] = -d1;
    skw[2] = -vJ_tmp;
    skw[5] = d1;
    skw[8] = 0.0;
    for (i3 = 0; i3 < 3; i3++) {
      vJ_tmp_tmp = skw[3 * i3];
      I_tmp[6 * i3] = vJ_tmp_tmp;
      loop_ub_tmp = 6 * (i3 + 3);
      I_tmp[loop_ub_tmp + 3] = vJ_tmp_tmp;
      vJ_tmp_tmp = skw[3 * i3 + 1];
      I_tmp[6 * i3 + 1] = vJ_tmp_tmp;
      I_tmp[loop_ub_tmp + 4] = vJ_tmp_tmp;
      vJ_tmp_tmp = skw[3 * i3 + 2];
      I_tmp[6 * i3 + 2] = vJ_tmp_tmp;
      I_tmp[loop_ub_tmp + 5] = vJ_tmp_tmp;
    }
    I_tmp[18] = 0.0;
    d2 = v_data[i + 2];
    I_tmp[24] = -d2;
    d3 = v_data[i + 1];
    I_tmp[30] = d3;
    I_tmp[19] = d2;
    I_tmp[25] = 0.0;
    d4 = v_data[i];
    I_tmp[31] = -d4;
    I_tmp[20] = -d3;
    I_tmp[26] = d4;
    I_tmp[32] = 0.0;
    for (i3 = 0; i3 < 6; i3++) {
      vJ_tmp_tmp = 0.0;
      for (j = 0; j < 6; j++) {
        vJ_tmp_tmp +=
            Xup_data[(i3 + 6 * j) + 36 * (b_i + 1)] * a_data[j + 6 * b_i];
      }
      b_Xup[i3] = vJ_tmp_tmp;
      c_I[i3] =
          (fh3[i3] * qd_tmp + fh3[i3 + 6] * b_qd_tmp) + fh3[i3 + 12] * c_qd_tmp;
    }
    for (i3 = 0; i3 < 6; i3++) {
      vJ_tmp_tmp = 0.0;
      for (j = 0; j < 6; j++) {
        vJ_tmp_tmp += I_tmp[i3 + 6 * j] * vJ[j];
      }
      a_data[i3 + i] = (b_Xup[i3] + c_I[i3]) + vJ_tmp_tmp;
    }
    memset(&I_tmp[0], 0, 36U * sizeof(double));
    skw[0] = 0.0;
    skw[3] = -d;
    skw[6] = vJ_tmp;
    skw[1] = d;
    skw[4] = 0.0;
    skw[7] = -d1;
    skw[2] = -vJ_tmp;
    skw[5] = d1;
    skw[8] = 0.0;
    for (i3 = 0; i3 < 3; i3++) {
      vJ_tmp_tmp = skw[3 * i3];
      I_tmp[6 * i3] = vJ_tmp_tmp;
      loop_ub_tmp = 6 * (i3 + 3);
      I_tmp[loop_ub_tmp + 3] = vJ_tmp_tmp;
      vJ_tmp_tmp = skw[3 * i3 + 1];
      I_tmp[6 * i3 + 1] = vJ_tmp_tmp;
      I_tmp[loop_ub_tmp + 4] = vJ_tmp_tmp;
      vJ_tmp_tmp = skw[3 * i3 + 2];
      I_tmp[6 * i3 + 2] = vJ_tmp_tmp;
      I_tmp[loop_ub_tmp + 5] = vJ_tmp_tmp;
    }
    I_tmp[18] = 0.0;
    I_tmp[24] = -d2;
    I_tmp[30] = d3;
    I_tmp[19] = d2;
    I_tmp[25] = 0.0;
    I_tmp[31] = -d4;
    I_tmp[20] = -d3;
    I_tmp[26] = d4;
    I_tmp[32] = 0.0;
    for (i3 = 0; i3 < 6; i3++) {
      for (j = 0; j < 6; j++) {
        b_I_tmp[j + 6 * i3] = -I_tmp[i3 + 6 * j];
      }
    }
    for (i3 = 0; i3 < 6; i3++) {
      c_I[i3] = 0.0;
      vJ[i3] = 0.0;
      for (j = 0; j < 6; j++) {
        d = 0.0;
        for (loop_ub_tmp = 0; loop_ub_tmp < 6; loop_ub_tmp++) {
          d += b_I_tmp[i3 + 6 * loop_ub_tmp] *
               I_data[(loop_ub_tmp + 6 * j) + 36 * (b_i + 1)];
        }
        loop_ub_tmp = j + i;
        c_I[i3] += I_data[(i3 + 6 * j) + 36 * (b_i + 1)] * a_data[loop_ub_tmp];
        vJ[i3] += d * v_data[loop_ub_tmp];
      }
      f_data[i3 + i] = c_I[i3] + vJ[i3];
    }
  }
  emxFree_real_T(&a);
  emxFree_real_T(&b_v);
  i = (int)-((-1.0 - N) + 1.0);
  emxInit_int32_T(&b_r, 1);
  for (b_i = 0; b_i < i; b_i++) {
    vJ_tmp = N - (double)b_i;
    if (vJ_tmp == 1.0) {
      vJ_tmp_tmp = 0.0;
      for (i2 = 0; i2 < 6; i2++) {
        vJ_tmp_tmp += (double)c_a[i2] * f_data[i2];
      }
      C_data[0] = vJ_tmp_tmp;
      for (i2 = 0; i2 < 6; i2++) {
        d = 0.0;
        for (i3 = 0; i3 < 6; i3++) {
          d += X[i2 + 6 * i3] * (double)b_a[i3];
        }
        J_data[i2] = d;
      }
    } else {
      vJ_tmp_tmp = (vJ_tmp - 2.0) * 3.0;
      if (vJ_tmp_tmp + 4.0 < vJ_tmp_tmp + 2.0) {
        qi->size[0] = 1;
        qi->size[1] = 0;
      } else if ((rtIsInf(vJ_tmp_tmp + 2.0) || rtIsInf(vJ_tmp_tmp + 4.0)) &&
                 (vJ_tmp_tmp + 2.0 == vJ_tmp_tmp + 4.0)) {
        i2 = qi->size[0] * qi->size[1];
        qi->size[0] = 1;
        qi->size[1] = 1;
        emxEnsureCapacity_real_T(qi, i2);
        qi_data = qi->data;
        qi_data[0] = rtNaN;
      } else {
        i2 = qi->size[0] * qi->size[1];
        qi->size[0] = 1;
        j = (int)((vJ_tmp_tmp + 4.0) - (vJ_tmp_tmp + 2.0));
        qi->size[1] = j + 1;
        emxEnsureCapacity_real_T(qi, i2);
        qi_data = qi->data;
        for (i2 = 0; i2 <= j; i2++) {
          qi_data[i2] = (vJ_tmp_tmp + 2.0) + (double)i2;
        }
      }
      i2 = b_q->size[0];
      b_q->size[0] = qi->size[1];
      emxEnsureCapacity_real_T(b_q, i2);
      q_data = b_q->data;
      j = qi->size[1];
      for (i2 = 0; i2 < j; i2++) {
        q_data[i2] = qi_data[i2];
      }
      i2 = b_r->size[0];
      b_r->size[0] = b_q->size[0];
      emxEnsureCapacity_int32_T(b_r, i2);
      r2 = b_r->data;
      j = b_q->size[0];
      for (i2 = 0; i2 < j; i2++) {
        r2[i2] = (int)q_data[i2] - 1;
      }
      for (i2 = 0; i2 < 3; i2++) {
        d = 0.0;
        for (i3 = 0; i3 < 6; i3++) {
          d += 2.0 * Smod_data[(i3 + 6 * i2) + 18 * ((int)vJ_tmp - 2)] *
               f_data[i3 + 6 * ((int)vJ_tmp - 1)];
        }
        c_qd[i2] = d;
      }
      j = b_r->size[0];
      for (i2 = 0; i2 < j; i2++) {
        C_data[r2[i2]] = c_qd[i2];
      }
      for (i2 = 0; i2 < 6; i2++) {
        d = 0.0;
        for (i3 = 0; i3 < 6; i3++) {
          d += Xup_data[(i3 + 6 * i2) + 36 * ((int)vJ_tmp - 1)] *
               f_data[i3 + 6 * ((int)vJ_tmp - 1)];
        }
        c_I[i2] = f_data[i2 + 6 * ((int)vJ_tmp - 2)] + d;
      }
      for (i2 = 0; i2 < 6; i2++) {
        f_data[i2 + 6 * ((int)vJ_tmp - 2)] = c_I[i2];
      }
      i2 = b_r->size[0];
      b_r->size[0] = b_q->size[0];
      emxEnsureCapacity_int32_T(b_r, i2);
      r2 = b_r->data;
      j = b_q->size[0];
      for (i2 = 0; i2 < j; i2++) {
        r2[i2] = (int)q_data[i2] - 1;
      }
      for (i2 = 0; i2 < 6; i2++) {
        for (i3 = 0; i3 < 3; i3++) {
          d = 0.0;
          for (j = 0; j < 6; j++) {
            d += X[i2 + 6 * j] *
                 Smod_data[(j + 6 * i3) + 18 * ((int)vJ_tmp - 2)];
          }
          fh3[i2 + 6 * i3] = d;
        }
      }
      j = b_r->size[0];
      for (i2 = 0; i2 < j; i2++) {
        for (i3 = 0; i3 < 6; i3++) {
          J_data[i3 + 6 * r2[i2]] = fh3[i3 + 6 * i2];
        }
      }
      for (i2 = 0; i2 < 6; i2++) {
        for (i3 = 0; i3 < 6; i3++) {
          d = 0.0;
          for (j = 0; j < 6; j++) {
            d += X[i2 + 6 * j] * XJ_data[(j + 6 * i3) + 36 * ((int)vJ_tmp - 1)];
          }
          I_tmp[i2 + 6 * i3] = d;
        }
        for (i3 = 0; i3 < 6; i3++) {
          d = 0.0;
          for (j = 0; j < 6; j++) {
            d += I_tmp[i2 + 6 * j] * (double)Xtree[j + 6 * i3];
          }
          X[i2 + 6 * i3] = d;
        }
      }
    }
  }
  emxFree_real_T(&f);
  emxFree_real_T(&XJ);
  /*  Composite Rigid Body Algorithm to calculate M */
  /*  composite inertia calculation */
  for (b_i = 0; b_i < i; b_i++) {
    vJ_tmp = N - (double)b_i;
    if (vJ_tmp != 1.0) {
      i2 = 36 * ((int)vJ_tmp - 1);
      for (i3 = 0; i3 < 6; i3++) {
        for (j = 0; j < 6; j++) {
          d = 0.0;
          for (loop_ub_tmp = 0; loop_ub_tmp < 6; loop_ub_tmp++) {
            d += Xup_data[(loop_ub_tmp + 6 * i3) + i2] *
                 I_data[(loop_ub_tmp + 6 * j) + i2];
          }
          I_tmp[i3 + 6 * j] = d;
        }
        for (j = 0; j < 6; j++) {
          d = 0.0;
          for (loop_ub_tmp = 0; loop_ub_tmp < 6; loop_ub_tmp++) {
            d += I_tmp[i3 + 6 * loop_ub_tmp] *
                 Xup_data[(loop_ub_tmp + 6 * j) + 36 * ((int)vJ_tmp - 1)];
          }
          loop_ub_tmp = i3 + 6 * j;
          b_I_tmp[loop_ub_tmp] =
              I_data[loop_ub_tmp + 36 * ((int)vJ_tmp - 2)] + d;
        }
      }
      for (i2 = 0; i2 < 6; i2++) {
        for (i3 = 0; i3 < 6; i3++) {
          j = i3 + 6 * i2;
          I_data[j + 36 * ((int)vJ_tmp - 2)] = b_I_tmp[j];
        }
      }
    }
  }
  i = M->size[0] * M->size[1];
  M->size[0] = (int)Nq_tmp;
  M->size[1] = (int)Nq_tmp;
  emxEnsureCapacity_real_T(M, i);
  XJ_data = M->data;
  j = (int)Nq_tmp * (int)Nq_tmp;
  for (i = 0; i < j; i++) {
    XJ_data[i] = 0.0;
  }
  /*  fh3 = zeros(6,3); */
  /*  fh1 = zeros(6,1); */
  emxInit_int32_T(&r1, 1);
  for (b_i = 0; b_i < i1; b_i++) {
    if ((unsigned int)b_i + 1U == 1U) {
      d = 0.0;
      for (i = 0; i < 6; i++) {
        vJ_tmp = 0.0;
        for (i2 = 0; i2 < 6; i2++) {
          vJ_tmp += I_data[i + 6 * i2] * (double)b_a[i2];
        }
        d += (double)iv1[i] * vJ_tmp;
      }
      XJ_data[0] = d;
    } else {
      unsigned int b_j;
      vJ_tmp_tmp = (((double)b_i + 1.0) - 2.0) * 3.0;
      if (vJ_tmp_tmp + 4.0 < vJ_tmp_tmp + 2.0) {
        qi->size[0] = 1;
        qi->size[1] = 0;
      } else {
        i = qi->size[0] * qi->size[1];
        qi->size[0] = 1;
        j = (int)((vJ_tmp_tmp + 4.0) - (vJ_tmp_tmp + 2.0));
        qi->size[1] = j + 1;
        emxEnsureCapacity_real_T(qi, i);
        qi_data = qi->data;
        for (i = 0; i <= j; i++) {
          qi_data[i] = (vJ_tmp_tmp + 2.0) + (double)i;
        }
      }
      for (i = 0; i < 6; i++) {
        for (i2 = 0; i2 < 3; i2++) {
          d = 0.0;
          for (i3 = 0; i3 < 6; i3++) {
            d += I_data[(i + 6 * i3) + 36 * b_i] *
                 Smod_data[(i3 + 6 * i2) + 18 * (b_i - 1)];
          }
          fh3[i + 6 * i2] = d;
        }
      }
      i = b_q->size[0];
      b_q->size[0] = qi->size[1];
      emxEnsureCapacity_real_T(b_q, i);
      q_data = b_q->data;
      j = qi->size[1];
      for (i = 0; i < j; i++) {
        q_data[i] = qi_data[i];
      }
      i = b_r->size[0];
      b_r->size[0] = b_q->size[0];
      emxEnsureCapacity_int32_T(b_r, i);
      r2 = b_r->data;
      j = b_q->size[0];
      i = r1->size[0];
      r1->size[0] = b_q->size[0];
      emxEnsureCapacity_int32_T(r1, i);
      r3 = r1->data;
      for (i = 0; i < j; i++) {
        i2 = (int)q_data[i] - 1;
        r2[i] = i2;
        r3[i] = i2;
      }
      for (i = 0; i < 3; i++) {
        for (i2 = 0; i2 < 3; i2++) {
          d = 0.0;
          for (i3 = 0; i3 < 6; i3++) {
            d += Smod_data[(i3 + 6 * i) + 18 * (b_i - 1)] * fh3[i3 + 6 * i2];
          }
          skw[i + 3 * i2] = d;
        }
      }
      loop_ub_tmp = b_r->size[0];
      j = r1->size[0];
      for (i = 0; i < j; i++) {
        for (i2 = 0; i2 < loop_ub_tmp; i2++) {
          XJ_data[r2[i2] + M->size[0] * r3[i]] = skw[i2 + loop_ub_tmp * i];
        }
      }
      b_j = (unsigned int)b_i + 1U;
      while ((b_j < N + 1.0) && (b_j > 1U)) {
        for (i = 0; i < 6; i++) {
          for (i2 = 0; i2 < 3; i2++) {
            d = 0.0;
            for (i3 = 0; i3 < 6; i3++) {
              d += Xup_data[(i3 + 6 * i) + 36 * ((int)b_j - 1)] *
                   fh3[i3 + 6 * i2];
            }
            c_Xup[i + 6 * i2] = d;
          }
        }
        memcpy(&fh3[0], &c_Xup[0], 18U * sizeof(double));
        b_j = (unsigned int)((int)b_j - 1);
        if ((int)b_j == 1) {
          i = b_r->size[0];
          b_r->size[0] = b_q->size[0];
          emxEnsureCapacity_int32_T(b_r, i);
          r2 = b_r->data;
          j = b_q->size[0];
          for (i = 0; i < j; i++) {
            r2[i] = (int)q_data[i] - 1;
          }
          for (i = 0; i < 3; i++) {
            d = 0.0;
            for (i2 = 0; i2 < 6; i2++) {
              d += fh3[i2 + 6 * i] * (double)b_a[i2];
            }
            c_qd[i] = d;
          }
          j = b_r->size[0];
          for (i = 0; i < j; i++) {
            XJ_data[r2[i]] = c_qd[i];
          }
          /* (S{j}' * fh).'; */
          i = b_r->size[0];
          b_r->size[0] = b_q->size[0];
          emxEnsureCapacity_int32_T(b_r, i);
          r2 = b_r->data;
          j = b_q->size[0];
          for (i = 0; i < j; i++) {
            r2[i] = (int)q_data[i] - 1;
          }
          for (i = 0; i < 3; i++) {
            d = 0.0;
            for (i2 = 0; i2 < 6; i2++) {
              d += (double)iv1[i2] * fh3[i2 + 6 * i];
            }
            c_qd[i] = d;
          }
          j = b_r->size[0];
          for (i = 0; i < j; i++) {
            XJ_data[M->size[0] * r2[i]] = c_qd[i];
          }
        } else {
          vJ_tmp_tmp = ((double)b_j - 2.0) * 3.0;
          if (vJ_tmp_tmp + 4.0 < vJ_tmp_tmp + 2.0) {
            qi->size[0] = 1;
            qi->size[1] = 0;
          } else {
            i = qi->size[0] * qi->size[1];
            qi->size[0] = 1;
            j = (int)((vJ_tmp_tmp + 4.0) - (vJ_tmp_tmp + 2.0));
            qi->size[1] = j + 1;
            emxEnsureCapacity_real_T(qi, i);
            qi_data = qi->data;
            for (i = 0; i <= j; i++) {
              qi_data[i] = (vJ_tmp_tmp + 2.0) + (double)i;
            }
          }
          i = b_r->size[0];
          b_r->size[0] = b_q->size[0];
          emxEnsureCapacity_int32_T(b_r, i);
          r2 = b_r->data;
          j = b_q->size[0];
          for (i = 0; i < j; i++) {
            r2[i] = (int)q_data[i] - 1;
          }
          i = b_qd->size[0];
          b_qd->size[0] = qi->size[1];
          emxEnsureCapacity_real_T(b_qd, i);
          qd_data = b_qd->data;
          j = qi->size[1];
          for (i = 0; i < j; i++) {
            qd_data[i] = qi_data[i];
          }
          i = r1->size[0];
          r1->size[0] = b_qd->size[0];
          emxEnsureCapacity_int32_T(r1, i);
          r3 = r1->data;
          j = b_qd->size[0];
          for (i = 0; i < j; i++) {
            r3[i] = (int)qd_data[i] - 1;
          }
          for (i = 0; i < 3; i++) {
            for (i2 = 0; i2 < 3; i2++) {
              d = 0.0;
              for (i3 = 0; i3 < 6; i3++) {
                d += fh3[i3 + 6 * i] *
                     Smod_data[(i3 + 6 * i2) + 18 * ((int)b_j - 2)];
              }
              skw[i + 3 * i2] = d;
            }
          }
          loop_ub_tmp = b_r->size[0];
          j = r1->size[0];
          for (i = 0; i < j; i++) {
            for (i2 = 0; i2 < loop_ub_tmp; i2++) {
              XJ_data[r2[i2] + M->size[0] * r3[i]] = skw[i2 + loop_ub_tmp * i];
            }
          }
          /* (S{j}' * fh).'; */
          i = b_r->size[0];
          b_r->size[0] = b_qd->size[0];
          emxEnsureCapacity_int32_T(b_r, i);
          r2 = b_r->data;
          j = b_qd->size[0];
          for (i = 0; i < j; i++) {
            r2[i] = (int)qd_data[i] - 1;
          }
          i = r1->size[0];
          r1->size[0] = b_q->size[0];
          emxEnsureCapacity_int32_T(r1, i);
          r3 = r1->data;
          j = b_q->size[0];
          for (i = 0; i < j; i++) {
            r3[i] = (int)q_data[i] - 1;
          }
          for (i = 0; i < 3; i++) {
            for (i2 = 0; i2 < 3; i2++) {
              d = 0.0;
              for (i3 = 0; i3 < 6; i3++) {
                d += Smod_data[(i3 + 6 * i) + 18 * ((int)b_j - 2)] *
                     fh3[i3 + 6 * i2];
              }
              skw[i + 3 * i2] = d;
            }
          }
          loop_ub_tmp = b_r->size[0];
          j = r1->size[0];
          for (i = 0; i < j; i++) {
            for (i2 = 0; i2 < loop_ub_tmp; i2++) {
              XJ_data[r2[i2] + M->size[0] * r3[i]] = skw[i2 + loop_ub_tmp * i];
            }
          }
        }
      }
    }
  }
  emxFree_real_T(&b_qd);
  emxFree_real_T(&b_q);
  emxFree_int32_T(&b_r);
  emxFree_int32_T(&r1);
  emxFree_real_T(&qi);
  emxFree_real_T(&Xup);
  emxFree_real_T(&Smod);
  emxFree_real_T(&b_I);
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
