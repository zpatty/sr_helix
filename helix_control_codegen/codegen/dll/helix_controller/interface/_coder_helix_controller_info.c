/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_info.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 20:27:52
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
      "789ced5acd72d250140e8ad6855556d6850b5db8aa6d68018b333a5348f9ab92f2ef8c8e"
      "1343724b02373f4d80a20b5fc1957b962ebbf4357c0c677c035d4808"
      "0192993b44c8a4e53667c1e1f071ef773897fbe5cc4d8850a1182208e21e61da8f6dd36f"
      "4ee2c8c4df20ece6c443137fc7115b768b08dbc659f8d789e714b90b",
      "065d339059094c47f28a24caacdcad7d5201a1015d817dc08f915311829a2881ea7c401b"
      "91949d83a6810119ef2901709d6a4f2234419f6508e783693d2e10bf"
      "37ecb21e05443d220efc7de603292812203fb39c40b6c4ee8e06544527758d11001407e4"
      "f89531caa42970e479d002b2fd5308b45dc991ffc715f3df40e66f22",
      "458a89315ccb33bedb483e13e1955e1382d9ef1baec87788e4b3e34baecfa43cbbd2a2fa"
      "dc7799afd3cfbe6feebc3fdb8f5b86f78baf5adedaf093cfb2cbe21b"
      "20e673fb7f7b80e08b38f03dbacf71fa9b52f52cd948b68f7ac958f324999de5515ac0b3"
      "280f0211fb35ff10317eddf6addb7c51d7c5c80461f936c3b190bb2c",
      "1dfdbe225f0ac967c7975c0fab3caeaf6f5eedff5f4febc04f3ecbae0b9f5f7a2a552850"
      "29d279f1dd59665f6e971b8d22dfc8e1a3a7b8ed5fb779871df12c6f"
      "13d13be0dc4b3ebffbd357483e3bbee4ba18e519af09ae7afae4e7971d3ff92cc35d4f13"
      "adfd4c5463e9d37a533f525f1fa4fa90cea5f1d1d32162fc55ed4fd5",
      "15f3759ecb38f3b5705d65bb223b4a445374e3bc645dfbd52c92cf8e2fababf365320436"
      "d0576ff82cc35d5f05f524a5f4d574593aa8d3b11789965a797b8e91"
      "bee2b67f837ed5b4a05f5d8d2fd053d3827ef5ffe61f22c6e3daafde5d90af8597288a69"
      "b39cd21459d988d7b55fcd20f9ecf892eb335f26acf5957cb6f5c84f",
      "3ecb70d7572eae76da0559e8f494349fd672f5e783463e8f8fbee2b67f83fb56a65db573"
      "ef75d555cbae0b5f70dfca9bf92f10e3afcbf3553791f99bc831a379"
      "cae7f779c04b249f1d5f727d46e53197c42f1dc8fc7de8ebf35587df7e6ffac96719eeba"
      "0aa3c7f97291dd4bd4bab578310ea3628d85d4faebea3fcea984e5",
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
                emlrtMxCreateDoubleScalar(739553.8524768519));
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
                emlrtMxCreateString("zyVEBdzOInuUR50avsek2"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_helix_controller_info.c
 *
 * [EOF]
 */
