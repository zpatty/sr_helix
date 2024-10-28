/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_info.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 17:40:15
 */

/* Include Files */
#include "_coder_helix_controller_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[7] = {
      "789ced594f73d240140f0e6a0ffee1a4170ff6a0176da3d0dae9c19902a5800eb150c019"
      "1d275d926d1248b231018a1c3c7bf32b70f4e857f12b78f3e647b021"
      "594232b343244c5ab679071e8f1fbbbf97f7d81f6f122655ada51886b9c73876fad4f177"
      "dd38e3fa1b8cdf8278caf5b70331b69b4cdab70ee3df5d2f20bd0f47",
      "7d27d08106672b45a4293ad0fbcd2f06644c68217508c52972a6a8b0a968f0643ee0ec48"
      "3b9a8366810dd9ef8b32147a27038d3165cbcb509d0f66f5f849b8de"
      "74c87a5409f5c804f08fa54fac8c34c88e8120b392d2df32a1812cd6327919aaca889dbe"
      "f276994ca45e78114a50f77faa42735b0be47f1a317f523f332e522b",
      "f2395e9056c6778bc8e720221a7454e85ddf2422df0191cf8f2fd91fb73cdbdaa2fadc0f"
      "996fd07bdfdf98fa67e9c792ede3e2fbf6f7f77e9c7cd82e8b6f44d8"
      "2fecefed01812f13c0cd6c7ddcdca98d7365f0863b94f27abedd7a59f2f2385ec0b3280f"
      "8610c7b5ff84b07eddce6dd87c3702b197af8300b1cb0b40152e4b47",
      "7f44e4cb13f9fcf892fdc0e509fdffb6aaf3ffe7490bc6c987edbaf0c5a5a75aa3081b35"
      "aea27cf85ccaeadd7abb5d13db657af494b6f31b36ef7420f6f27610"
      "ab07cf57c917f77cfa9ac8e7c797ec8b5d9e694f68d5d3cd5f5fb7e2e4c346bb9eee4ad9"
      "d20b137067ad8e7568bcddcb0f55ae5ca0474f2784f557753e3522e6",
      "1bbc2f13cc17e39601fa0ab848c444967dbf645de7d523229f1f5f5657e7cb640b6ca2af"
      "abe1c346bbbecac6bb3c1a1a85bab6d7e272fbbb92d1787f4e91bed2"
      "767e9379d5b1645e8dc697e8a963c9bcfa7ffb4f08eb699d57ef2cc817e3c7c522df0502"
      "ea2840b7e3759d574b443e3fbe647fe6cb44b5beb2cf1f3e8a930f1b",
      "edfa2aec18bd6e55977b0354100b66b9f56ad4ae54e8d157dace6ff2dccab1ab76df7b5d"
      "7515db75e14b9e5b45dbff1fdc113a0e",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 9848U, &nameCaptureInfo);
  return nameCaptureInfo;
}

/*
 * Arguments    : void
 * Return Type  : mxArray *
 */
mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9] = {"Version",
                                    "ResolvedFunctions",
                                    "Checksum",
                                    "EntryPoints",
                                    "CoverageInfo",
                                    "IsPolymorphic",
                                    "PropertyList",
                                    "UUID",
                                    "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8] = {
      "Name",     "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "FullPath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 24);
  emlrtSetField(xEntryPoints, 0, "Name",
                emlrtMxCreateString("helix_controller"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(24.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "FullPath",
      emlrtMxCreateString("/home/zach/git-repos/sr_helix/helix_control_codegen/"
                          "helix_controller.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739553.72140046291));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("23.2.0.2428915 (R2023b) Update 4"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("8wg3TCGbIGNOOz86czPxp"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_helix_controller_info.c
 *
 * [EOF]
 */
