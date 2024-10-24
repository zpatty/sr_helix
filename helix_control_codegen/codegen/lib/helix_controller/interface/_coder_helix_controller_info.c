/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_info.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 23-Oct-2024 10:17:18
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
      "789ced5acb92d240140d233e168eb2f251e5c6851bcb210c0fd172e3c0f074c0e1a93596"
      "9569921e12e83c4c02c3b8f0177433552ee797dccd2f58fa07ba7030"
      "09905475110aaa8190bba0b975489fc3bdf44993840a144a018aa2ee50465c3e35c66d33"
      "0f99e316650f271e30c75b8edc8aeb54d0769c857f354756967438d0",
      "8d4402221c1dc9c9a2200149af9f299052a126a33ee4fe23270282754184b5c9a43cccc4"
      "ec04344a86d0f07d9a876cb7d6132995d7c60ad16432aac705e6fb06"
      "5dd6e315a61e2107fe21f391e66511d29f01cbd36d41df51a1226bb40874045af47ef580"
      "e6211206c62b332c972a2304d5b068d37b3ca7de9b58bd06524a3331",
      "866d2f8cef0696cf4038b9d74270fcfdbecfc917c5f2d971d7fd300b62b6616a5deebad4"
      "e91cc79f3756d8e5d6c33649be5f7fce7749f259b12cbe01663eb7bf"
      "b37b18be9003ef4a07f16c2772a6d59bb542315f7e9f3d6aa5f2631d875378a6e9a03039"
      "a9f9d77dbdbad5893bef854c04701d8605885d57df8c61f9ecb8eb3e",
      "58057179fe5ad43afff9a40149f259b1297ca47c53aca661b554ce0b479f3251a9536936"
      "4b5c33e7fbe6aaac57b73a838e7cacd340b42e3c5d24dfacbef96d4e"
      "3e1acb67c75df7615890ab1e78d52f1ffff8b24392cf0aaffb65a21dcd4454503e69b4b4"
      "7de54d72af8fcab994ef97a4f799ca9c3a9dd74f9c3a2d5c53802e00",
      "c4b0aaac69d4f2fc73debebcc0f2d971f7fe395998b0e8fbe862f8acf0ba8ff2cadb3db9"
      "afa42a62b2518ebd4cb495eabb53df475766bdfafb4e23fc7de76c7c"
      "be5f1ae1ef3b679b7f53f69db7a7e8b4f0c3749ae900566e09409ae43f9e939ff4be3389"
      "e5b3e3aefb3259180ffb28fdecfe23927c5678dd47d9b8d2ed1424be",
      "db93535c4acd359e0f9a79ff3ed1caac57ff3e9111cbbeeebcaebe69c5a6f0f9f7891633"
      "ff05e678af3e9f740dabd7408a8cba503ed2ffdbc3583e3beeba1f57"
      "0519b680d43acffc7d40f4b9a4d7e7bfb749f259e175df449162be5202bb89ba5e8f97e2"
      "2822d4014aafbf6ffe0364221d19",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 10992U, &nameCaptureInfo);
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
      emlrtMxCreateString(
          "/home/zach/git-repos/matlab/DRL/helix/helix_controller.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739548.42464120372));
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
                emlrtMxCreateString("5KhMDmcLi57OLkh3TREVt"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_helix_controller_info.c
 *
 * [EOF]
 */
