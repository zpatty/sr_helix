/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_api.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 24-Oct-2024 11:42:25
 */

/* Include Files */
#include "_coder_helix_controller_api.h"
#include "_coder_helix_controller_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131643U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "helix_controller",                                   /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[10];

static const mxArray *b_emlrt_marshallOut(const real_T u[3]);

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier);

static const mxArray *c_emlrt_marshallOut(const real_T u[100]);

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[100];

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier))[10];

static const mxArray *emlrt_marshallOut(const real_T u[10]);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[100];

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3];

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3];

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[10];

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[100];

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3];

/* Function Definitions */
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[10]
 */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[10]
{
  real_T(*y)[10];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const real_T u[3]
 * Return Type  : const mxArray *
 */
static const mxArray *b_emlrt_marshallOut(const real_T u[3])
{
  static const int32_T i = 0;
  static const int32_T i1 = 3;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const real_T u[100]
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrt_marshallOut(const real_T u[100])
{
  static const int32_T iv[2] = {0, 0};
  static const int32_T iv1[2] = {10, 10};
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, &iv1[0], 2);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real_T (*)[100]
 */
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[100]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[100];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real_T (*)[10]
 */
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier))[10]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[10];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const real_T u[10]
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real_T u[10])
{
  static const int32_T i = 0;
  static const int32_T i1 = 10;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[100]
 */
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[100]
{
  real_T(*y)[100];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real_T (*)[3]
 */
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[3];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[3]
 */
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3]
{
  real_T(*y)[3];
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[10]
 */
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[10]
{
  static const int32_T dims = 10;
  real_T(*ret)[10];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[10])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[100]
 */
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[100]
{
  static const int32_T dims[2] = {10, 10};
  real_T(*ret)[100];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[100])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[3]
 */
static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3]
{
  static const int32_T dims = 3;
  real_T(*ret)[3];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[3])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const mxArray * const prhs[23]
 *                int32_T nlhs
 *                const mxArray *plhs[7]
 * Return Type  : void
 */
void helix_controller_api(const mxArray *const prhs[23], int32_T nlhs,
                          const mxArray *plhs[7])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  real_T(*A)[100];
  real_T(*KD)[100];
  real_T(*Kp)[100];
  real_T(*M)[100];
  real_T(*C)[10];
  real_T(*ddqd)[10];
  real_T(*dq)[10];
  real_T(*dqd)[10];
  real_T(*q)[10];
  real_T(*qd)[10];
  real_T(*tau)[10];
  real_T(*tau_r)[10];
  real_T(*cq)[3];
  real_T(*dxd)[3];
  real_T(*dxr)[3];
  real_T(*x)[3];
  real_T(*xd)[3];
  real_T KDx;
  real_T Kpx;
  real_T L0;
  real_T bb;
  real_T bm;
  real_T bs;
  real_T conv_motor;
  real_T conv_pcc;
  real_T d;
  real_T kb;
  real_T ks;
  real_T m;
  real_T r;
  st.tls = emlrtRootTLSGlobal;
  tau = (real_T(*)[10])mxMalloc(sizeof(real_T[10]));
  tau_r = (real_T(*)[10])mxMalloc(sizeof(real_T[10]));
  x = (real_T(*)[3])mxMalloc(sizeof(real_T[3]));
  M = (real_T(*)[100])mxMalloc(sizeof(real_T[100]));
  C = (real_T(*)[10])mxMalloc(sizeof(real_T[10]));
  A = (real_T(*)[100])mxMalloc(sizeof(real_T[100]));
  cq = (real_T(*)[3])mxMalloc(sizeof(real_T[3]));
  /* Marshall function inputs */
  q = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "q");
  dq = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "dq");
  qd = emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "qd");
  dqd = emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "dqd");
  ddqd = emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "ddqd");
  d = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "d");
  m = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "m");
  r = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "r");
  kb = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "kb");
  ks = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "ks");
  bb = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "bb");
  bs = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "bs");
  bm = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "bm");
  L0 = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "L0");
  Kp = e_emlrt_marshallIn(&st, emlrtAlias(prhs[14]), "Kp");
  KD = e_emlrt_marshallIn(&st, emlrtAlias(prhs[15]), "KD");
  Kpx = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[16]), "Kpx");
  KDx = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[17]), "KDx");
  xd = g_emlrt_marshallIn(&st, emlrtAlias(prhs[18]), "xd");
  dxd = g_emlrt_marshallIn(&st, emlrtAlias(prhs[19]), "dxd");
  dxr = g_emlrt_marshallIn(&st, emlrtAlias(prhs[20]), "dxr");
  conv_pcc = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[21]), "conv_pcc");
  conv_motor = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[22]), "conv_motor");
  /* Invoke the target function */
  helix_controller(*q, *dq, *qd, *dqd, *ddqd, d, m, r, kb, ks, bb, bs, bm, L0,
                   *Kp, *KD, Kpx, KDx, *xd, *dxd, *dxr, conv_pcc, conv_motor,
                   *tau, *tau_r, *x, *M, *C, *A, *cq);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*tau);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(*tau_r);
  }
  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(*x);
  }
  if (nlhs > 3) {
    plhs[3] = c_emlrt_marshallOut(*M);
  }
  if (nlhs > 4) {
    plhs[4] = emlrt_marshallOut(*C);
  }
  if (nlhs > 5) {
    plhs[5] = c_emlrt_marshallOut(*A);
  }
  if (nlhs > 6) {
    plhs[6] = b_emlrt_marshallOut(*cq);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void helix_controller_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  helix_controller_xil_terminate();
  helix_controller_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void helix_controller_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void helix_controller_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_helix_controller_api.c
 *
 * [EOF]
 */
