/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_info.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 29-Oct-2024 18:41:09
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
      "789ced5acd72d250140e8a3f0babacd4850b5db8d23608c5ea8cce14520ab4060b053aad"
      "e3c44b729b04f2671228ea8c2fe0c295ae19572ebbf4357c0c677c04"
      "17129200c9cc1d109820979c0587c3c7bddfe15cee9733372122053a4210c475c2b61f0f"
      "6cbfe6c431c75f20bce6c7238ebfea8b5dbb44443de35cfcb3e35955",
      "3161c7b40305c8703092536551018a5979a7414287862ab521d7474e45095644191e8e06"
      "452b927747a0416041d67b4a806cf3b02513ba600c33944683413dce"
      "11bf373a613d0a887ac47cf8abec6b52506548be07ac40f2a2b9ae434d35484367042889"
      "1db2ffca5865d255a9e739c843c5fba904f50dd997ff9b19f3bf82cc",
      "df46688a49322c3f37becb483e1be1d4565d82c3dfd79d916f1bc9e7c5a75c1fa73c1bf2"
      "b8fadc98305fbf1f7edfde795f5fdce52d1f14df874fdf9e04c9e7da"
      "a2f83a88f926fdbfdd44f0c57c78833d3ece7372219b3472d5c44e217d04f6a9cc308f83"
      "313ce3f220107150f37711e3976ddf4e9a2feaba187310c035181648",
      "eca274f4fb8c7c69249f179f723ddcf24c7c7d9bd7feff75bf0a83e4736d55f882d253b9"
      "4cc1325dcc8b276fb309a551aad568ae96c3474f71dbbf93e61df5c5"
      "c3bc6dc468c2b379f205dd9f3e47f279f129d7c52a4f7f4d70d5d37b3f3fae07c9e71aee"
      "7a9ae213d9b80e8aa7d5bab1a3ed6fa5db523117f6a70beb4fb519f3",
      "f59fcbf8f375714303a6087a89e8aa619d972c6bbfba8be4f3e2d3eaea68992c810df575"
      "3e7caee1aeaf82f632adb6b54c49deaa16934f53bc563e3ac3485f71"
      "dbbf61bf6a5bd8afcec617eaa96d61bffa6ff37711e371ed57af8dc9d7c50f288a690056"
      "ad8b40b1e265ed57b3483e2f3ee5fa8c96096b7d251fdeba13249f6b",
      "b8eb2bbba9351b054568b6d40c97d173d5c79d5a3e8f8fbee2b67fc3fb56b6fd6fe7decb"
      "aaabaead0a5f78df6a3ef39f23c6afcaf3551791f9dbc81ea3cf952f"
      "e8f38067483e2f3ee5faf4ca632f49503a90fd733bd0e7abb6bffc5e0b92cf35dc75558a"
      "efe54b347894aa98954d7a538a8b152051cbafab7f01764c859e",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 11216U, &nameCaptureInfo);
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
  xInputs = emlrtCreateLogicalMatrix(1, 23);
  emlrtSetField(xEntryPoints, 0, "Name",
                emlrtMxCreateString("helix_controller"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(23.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "FullPath",
      emlrtMxCreateString("/home/zach/git-repos/sr_helix/helix_control_codegen/"
                          "helix_controller.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739554.77636574069));
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
                emlrtMxCreateString("AIJ61yIlNCnrQJhg6bEiGF"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_helix_controller_info.c
 *
 * [EOF]
 */
