/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_info.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 24-Oct-2024 11:42:25
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
  const char_T *data[8] = {
      "789ced5a4b73d250140e151f0babac7cccb871e14a5b680b769cd19902f254b03c75da71"
      "62486e49e0e6610294baf02fb8ea9ea5cb2efd1beefc0b8efe035d48"
      "7ab9403273871898b45c72161c0e1ff77e8773b95fcedc8409e40a0186616e31c8be3d46"
      "7e7d1487467e8db19a1d0f8cfc0d5b8ced2a13b48cc3f89791e755a5",
      "03fa1d14289c0cc623055596144ee9544f34c0e8c050610f08e7c8910441559241653a28"
      "9a919c9e82c6810999ef9322e0db95aecce8a231c9104e07e37a9c11"
      "7e6fd0613d72847a846cf861ea7d58546510fec4f162b829753674a0a946d8d0591140a9"
      "1f3e7f65cd32e92a1c7a01348162fd14027d53b6e5ff61cefcaf13f3",
      "474821c9eeb07c73617cd7887c0811d46e0382c9ef1bccc9b747e4b3e22ed767549e4d79"
      "567d6e3bccd7ee27df473befc7dafda6e9bde2fbf5e774cb4b3e6c17"
      "c5d727cce7f4ff7687c017b2e16de57534dd8a9c18d57a2597cf16dfa50f1a89ec248ffd"
      "193cb3f26008b157f30f08e3976ddf3acd97745d0c8d104e68b13c07",
      "f98bd2d1af73f2c5897c56dce57ae0f238bebe2d6affff7c54035ef2615b153eaff4542e"
      "2741b950cc4a071f53db4aab54af17847a861e3da56dff3acd3b688b"
      "277923c46883e345f279dd9fbe20f2597197eb6296e77c4d68d5d387df3f6f78c9878d76"
      "3d8d35b753119d2b1ed51ac64bedd56ebc078b99043d7a3a208cbfac",
      "fda93667bef673197bbe183734ae2371c34474d530cf4b96b55f4d13f9acb85b5d9d2e93"
      "29b0bebe2e860f1bedfa2a6a6fe26a4f4b94e4dd5a71e759aca995df"
      "1e53a4afb4ed5fbf5f45e6f7abf3f1f97a8accef57ff6ffe01613cadfdeacd19f9627c3f"
      "99645b1caf36244e31e365ed5753443e2bee727da6cb44b5be869fdc",
      "7de0251f36daf5958f6aed564e11db5d352124f44ced69bf9ea5e87e156dfbd7bf6f85ec"
      "b29d7b2fabae625b153effbed562e63f238c5f95e7abae10f347489e"
      "d517cae7f579c073229f1577b93ec3f2a025f14a07527fef79fa7cd5dee9ef752ff9b0d1"
      "aeab3092cf960adc56acdaa9460b511891aa1c4c2ebfaefe030aa486",
      "9a",
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
