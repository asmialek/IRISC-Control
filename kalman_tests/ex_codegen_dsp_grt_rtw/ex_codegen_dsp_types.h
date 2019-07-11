/*
 * ex_codegen_dsp_types.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "ex_codegen_dsp".
 *
 * Model version              : 1.21
 * Simulink Coder version : 9.1 (R2019a) 23-Nov-2018
 * C source code generated on : Thu May 23 10:23:32 2019
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_ex_codegen_dsp_types_h_
#define RTW_HEADER_ex_codegen_dsp_types_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#ifndef struct_md128433f1ef248ca39112e85ac51f2a3
#define struct_md128433f1ef248ca39112e85ac51f2a3

struct md128433f1ef248ca39112e85ac51f2a3
{
  int32_T S0_isInitialized;
  real_T W0_states[38];
  real_T P0_InitialStates;
  real_T P1_Coefficients[39];
};

#endif                              /*struct_md128433f1ef248ca39112e85ac51f2a3*/

#ifndef typedef_dsp_FIRFilter_0_ex_codegen_ds_T
#define typedef_dsp_FIRFilter_0_ex_codegen_ds_T

typedef struct md128433f1ef248ca39112e85ac51f2a3 dsp_FIRFilter_0_ex_codegen_ds_T;

#endif                               /*typedef_dsp_FIRFilter_0_ex_codegen_ds_T*/

#ifndef struct_md5Fbh9n3EThqnsmmDPSZQHD
#define struct_md5Fbh9n3EThqnsmmDPSZQHD

struct md5Fbh9n3EThqnsmmDPSZQHD
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  dsp_FIRFilter_0_ex_codegen_ds_T cSFunObject;
};

#endif                                 /*struct_md5Fbh9n3EThqnsmmDPSZQHD*/

#ifndef typedef_dspcodegen_FIRFilter_ex_codeg_T
#define typedef_dspcodegen_FIRFilter_ex_codeg_T

typedef struct md5Fbh9n3EThqnsmmDPSZQHD dspcodegen_FIRFilter_ex_codeg_T;

#endif                               /*typedef_dspcodegen_FIRFilter_ex_codeg_T*/

#ifndef typedef_c_cell_wrap_ex_codegen_dsp_T
#define typedef_c_cell_wrap_ex_codegen_dsp_T

typedef struct {
  uint32_T f1[8];
} c_cell_wrap_ex_codegen_dsp_T;

#endif                                 /*typedef_c_cell_wrap_ex_codegen_dsp_T*/

#ifndef struct_mdipIcnF8zxxhFu1C3kGyN8D
#define struct_mdipIcnF8zxxhFu1C3kGyN8D

struct mdipIcnF8zxxhFu1C3kGyN8D
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  c_cell_wrap_ex_codegen_dsp_T inputVarSize;
  int32_T NumChannels;
  dspcodegen_FIRFilter_ex_codeg_T *FilterObj;
};

#endif                                 /*struct_mdipIcnF8zxxhFu1C3kGyN8D*/

#ifndef typedef_dsp_LowpassFilter_ex_codegen__T
#define typedef_dsp_LowpassFilter_ex_codegen__T

typedef struct mdipIcnF8zxxhFu1C3kGyN8D dsp_LowpassFilter_ex_codegen__T;

#endif                               /*typedef_dsp_LowpassFilter_ex_codegen__T*/

/* Parameters (default storage) */
typedef struct P_ex_codegen_dsp_T_ P_ex_codegen_dsp_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_ex_codegen_dsp_T RT_MODEL_ex_codegen_dsp_T;

#endif                                 /* RTW_HEADER_ex_codegen_dsp_types_h_ */
