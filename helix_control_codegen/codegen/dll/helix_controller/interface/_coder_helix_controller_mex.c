/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_mex.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 20:27:52
 */

/* Include Files */
#include "_coder_helix_controller_mex.h"
#include "_coder_helix_controller_api.h"

/* Function Definitions */
/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[]
 *                int32_T nrhs
 *                const mxArray *prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&helix_controller_atexit);
  /* Module initialization. */
  helix_controller_initialize();
  /* Dispatch the entry-point. */
  unsafe_helix_controller_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  helix_controller_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[7]
 *                int32_T nrhs
 *                const mxArray *prhs[23]
 * Return Type  : void
 */
void unsafe_helix_controller_mexFunction(int32_T nlhs, mxArray *plhs[7],
                                         int32_T nrhs, const mxArray *prhs[23])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *b_prhs[23];
  const mxArray *outputs[7];
  int32_T i;
  int32_T i1;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 23) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 23, 4,
                        16, "helix_controller");
  }
  if (nlhs > 7) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 16,
                        "helix_controller");
  }
  /* Call the function. */
  for (i = 0; i < 23; i++) {
    b_prhs[i] = prhs[i];
  }
  helix_controller_api(b_prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i1 = 1;
  } else {
    i1 = nlhs;
  }
  emlrtReturnArrays(i1, &plhs[0], &outputs[0]);
}

/*
 * File trailer for _coder_helix_controller_mex.c
 *
 * [EOF]
 */
