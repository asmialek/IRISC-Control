/*
 * ex_codegen_dsp_data.c
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

#include "ex_codegen_dsp.h"
#include "ex_codegen_dsp_private.h"

/* Block parameters (default storage) */
P_ex_codegen_dsp_T ex_codegen_dsp_P = {
  /* Mask Parameter: SineWave_Amplitude
   * Referenced by: '<Root>/Sine Wave'
   */
  1.0,

  /* Mask Parameter: SineWave_Frequency
   * Referenced by: '<Root>/Sine Wave'
   */
  0.5,

  /* Mask Parameter: RandomSource_MeanVal
   * Referenced by: '<Root>/Random Source'
   */
  0.0,

  /* Mask Parameter: SineWave_Phase
   * Referenced by: '<Root>/Sine Wave'
   */
  0.0,

  /* Mask Parameter: LMSFilter_leakage
   * Referenced by: '<Root>/LMS Filter'
   */
  1.0,

  /* Mask Parameter: LMSFilter_mu
   * Referenced by: '<Root>/LMS Filter'
   */
  0.1,

  /* Mask Parameter: RandomSource_rawSeed
   * Referenced by: '<Root>/Random Source'
   */
  23341U,

  /* Computed Parameter: RandomSource_VarianceRTP
   * Referenced by: '<Root>/Random Source'
   */
  1.0
};
