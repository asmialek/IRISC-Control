/*
 * ex_codegen_dsp.c
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

/* Block states (default storage) */
DW_ex_codegen_dsp_T ex_codegen_dsp_DW;

/* Real-time model */
RT_MODEL_ex_codegen_dsp_T ex_codegen_dsp_M_;
RT_MODEL_ex_codegen_dsp_T *const ex_codegen_dsp_M = &ex_codegen_dsp_M_;

/* Forward declaration for local functions */
static void ex_codegen_d_SystemCore_release(dspcodegen_FIRFilter_ex_codeg_T *obj);
static void ex_c_LPHPFilterBase_releaseImpl(dsp_LowpassFilter_ex_codegen__T *obj);
static void ex_co_SystemCore_releaseWrapper(dsp_LowpassFilter_ex_codegen__T *obj);
static void ex_codegen_SystemCore_release_l(dsp_LowpassFilter_ex_codegen__T *obj);
static void ex_codegen__SystemCore_delete_l(dsp_LowpassFilter_ex_codegen__T *obj);
static void matlabCodegenHandle_matlabCod_l(dsp_LowpassFilter_ex_codegen__T *obj);
static void ex_codegen_ds_SystemCore_delete(dspcodegen_FIRFilter_ex_codeg_T *obj);
static void matlabCodegenHandle_matlabCodeg(dspcodegen_FIRFilter_ex_codeg_T *obj);
void RandSrcInitState_GZ(const uint32_T seed[], uint32_T state[], int32_T nChans)
{
  int32_T i;
  int32_T tmp;

  /* InitializeConditions for S-Function (sdsprandsrc2): '<Root>/Random Source' */
  /* RandSrcInitState_GZ */
  for (i = 0; i < nChans; i++) {
    tmp = i << 1;
    state[tmp] = 362436069U;
    state[tmp + 1] = seed[i] == 0U ? 521288629U : seed[i];
  }

  /* End of InitializeConditions for S-Function (sdsprandsrc2): '<Root>/Random Source' */
}

void RandSrc_GZ_D(real_T y[], const real_T mean[], int32_T meanLen, const real_T
                  xstd[], int32_T xstdLen, uint32_T state[], int32_T nChans,
                  int32_T nSamps)
{
  int32_T i;
  int32_T j;
  real_T r;
  real_T x;
  real_T s;
  real_T y_0;
  int32_T chan;
  real_T std;
  uint32_T icng;
  int32_T samp;
  static const real_T vt[65] = { 0.340945, 0.4573146, 0.5397793, 0.6062427,
    0.6631691, 0.7136975, 0.7596125, 0.8020356, 0.8417227, 0.8792102, 0.9148948,
    0.9490791, 0.9820005, 1.0138492, 1.044781, 1.0749254, 1.1043917, 1.1332738,
    1.161653, 1.189601, 1.2171815, 1.2444516, 1.2714635, 1.298265, 1.3249008,
    1.3514125, 1.3778399, 1.4042211, 1.4305929, 1.4569915, 1.4834527, 1.5100122,
    1.5367061, 1.5635712, 1.5906454, 1.617968, 1.6455802, 1.6735255, 1.7018503,
    1.7306045, 1.7598422, 1.7896223, 1.8200099, 1.851077, 1.8829044, 1.9155831,
    1.9492166, 1.9839239, 2.0198431, 2.0571356, 2.095993, 2.136645, 2.1793713,
    2.2245175, 2.2725186, 2.3239338, 2.3795008, 2.4402218, 2.5075117, 2.5834658,
    2.6713916, 2.7769942, 2.7769942, 2.7769942, 2.7769942 };

  int32_T icng_tmp;
  int32_T jsr_tmp;
  uint32_T jsr;

  /* S-Function (sdsprandsrc2): '<Root>/Random Source' */
  /* RandSrc_GZ_D */
  for (chan = 0; chan < nChans; chan++) {
    std = xstd[xstdLen > 1 ? chan : 0];
    icng_tmp = chan << 1;
    icng = state[icng_tmp];
    jsr_tmp = icng_tmp + 1;
    jsr = state[jsr_tmp];
    for (samp = 0; samp < nSamps; samp++) {
      icng = 69069U * icng + 1234567U;
      jsr ^= jsr << 13;
      jsr ^= jsr >> 17;
      jsr ^= jsr << 5;
      i = (int32_T)(icng + jsr);
      j = (i & 63) + 1;
      r = (real_T)i * 4.6566128730773926E-10 * vt[j];
      x = fabs(r);
      y_0 = vt[j - 1];
      if (!(x <= y_0)) {
        x = (x - y_0) / (vt[j] - y_0);
        icng = 69069U * icng + 1234567U;
        jsr ^= jsr << 13;
        jsr ^= jsr >> 17;
        jsr ^= jsr << 5;
        y_0 = (real_T)(int32_T)(icng + jsr) * 2.328306436538696E-10 + 0.5;
        s = x + y_0;
        if (s > 1.301198) {
          r = r < 0.0 ? 0.4878992 * x - 0.4878992 : 0.4878992 - 0.4878992 * x;
        } else {
          if (!(s <= 0.9689279)) {
            x = 0.4878992 - 0.4878992 * x;
            if (y_0 > 12.67706 - exp(-0.5 * x * x) * 12.37586) {
              r = r < 0.0 ? -x : x;
            } else {
              if (!(exp(-0.5 * vt[j] * vt[j]) + y_0 * 0.01958303 / vt[j] <= exp(
                    -0.5 * r * r))) {
                do {
                  icng = 69069U * icng + 1234567U;
                  jsr ^= jsr << 13;
                  jsr ^= jsr >> 17;
                  jsr ^= jsr << 5;
                  x = log((real_T)(int32_T)(icng + jsr) * 2.328306436538696E-10
                          + 0.5) / 2.776994;
                  icng = 69069U * icng + 1234567U;
                  jsr ^= jsr << 13;
                  jsr ^= jsr >> 17;
                  jsr ^= jsr << 5;
                } while (log((real_T)(int32_T)(icng + jsr) *
                             2.328306436538696E-10 + 0.5) * -2.0 <= x * x);

                r = r < 0.0 ? x - 2.776994 : 2.776994 - x;
              }
            }
          }
        }
      }

      y[chan * nSamps + samp] = mean[meanLen > 1 ? chan : 0] + std * r;
    }

    state[icng_tmp] = icng;
    state[jsr_tmp] = jsr;
  }

  /* End of S-Function (sdsprandsrc2): '<Root>/Random Source' */
}

void MWSPCGlmsnw_D(const real_T x[], const real_T d[], real_T mu, uint32_T
                   *startIdx, real_T xBuf[], real_T wBuf[], int32_T wLen, real_T
                   leakFac, int32_T xLen, real_T y[], real_T eY[], real_T wY[])
{
  int32_T i;
  int32_T j;
  real_T bufEnergy;
  int32_T j1;

  /* S-Function (sdsplms): '<Root>/LMS Filter' */
  for (i = 0; i < xLen; i++) {
    y[i] = 0.0;
  }

  for (i = 0; i < xLen; i++) {
    bufEnergy = 0.0;

    /* Copy the current sample at the END of the circular buffer and update BuffStartIdx
     */
    xBuf[*startIdx] = x[i];
    (*startIdx)++;
    if (*startIdx == (uint32_T)wLen) {
      *startIdx = 0U;
    }

    /* Multiply wgtBuff_vector (not yet updated) and inBuff_vector
     */
    /* Get the energy of the signal in updated buffer
     */
    j1 = 0;
    for (j = (int32_T)*startIdx; j < wLen; j++) {
      y[i] += wBuf[j1] * xBuf[j];
      bufEnergy += xBuf[j] * xBuf[j];
      j1++;
    }

    for (j = 0; j < (int32_T)*startIdx; j++) {
      y[i] += wBuf[j1] * xBuf[j];
      bufEnergy += xBuf[j] * xBuf[j];
      j1++;
    }

    /* Ger error for the current sample
     */
    eY[i] = d[i] - y[i];

    /* Update weight-vector for next input sample
     */
    j1 = 0;
    for (j = (int32_T)*startIdx; j < wLen; j++) {
      wBuf[j1] = xBuf[j] / (bufEnergy + 2.2204460492503131E-16) * eY[i] * mu +
        leakFac * wBuf[j1];
      j1++;
    }

    for (j = 0; j < (int32_T)*startIdx; j++) {
      wBuf[j1] = xBuf[j] / (bufEnergy + 2.2204460492503131E-16) * eY[i] * mu +
        leakFac * wBuf[j1];
      j1++;
    }
  }

  j1 = wLen;
  for (j = 0; j < wLen; j++) {
    wY[j] = wBuf[j1 - 1];
    j1--;
  }

  /* End of S-Function (sdsplms): '<Root>/LMS Filter' */
}

static void ex_codegen_d_SystemCore_release(dspcodegen_FIRFilter_ex_codeg_T *obj)
{
  if (obj->isInitialized == 1) {
    obj->isInitialized = 2;
  }
}

static void ex_c_LPHPFilterBase_releaseImpl(dsp_LowpassFilter_ex_codegen__T *obj)
{
  ex_codegen_d_SystemCore_release(obj->FilterObj);
  obj->NumChannels = -1;
}

static void ex_co_SystemCore_releaseWrapper(dsp_LowpassFilter_ex_codegen__T *obj)
{
  if (obj->isSetupComplete) {
    ex_c_LPHPFilterBase_releaseImpl(obj);
  }
}

static void ex_codegen_SystemCore_release_l(dsp_LowpassFilter_ex_codegen__T *obj)
{
  if (obj->isInitialized == 1) {
    ex_co_SystemCore_releaseWrapper(obj);
  }
}

static void ex_codegen__SystemCore_delete_l(dsp_LowpassFilter_ex_codegen__T *obj)
{
  ex_codegen_SystemCore_release_l(obj);
}

static void matlabCodegenHandle_matlabCod_l(dsp_LowpassFilter_ex_codegen__T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
    ex_codegen__SystemCore_delete_l(obj);
  }
}

static void ex_codegen_ds_SystemCore_delete(dspcodegen_FIRFilter_ex_codeg_T *obj)
{
  ex_codegen_d_SystemCore_release(obj);
}

static void matlabCodegenHandle_matlabCodeg(dspcodegen_FIRFilter_ex_codeg_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
    ex_codegen_ds_SystemCore_delete(obj);
  }
}

/* Model output function */
static void ex_codegen_dsp_output(void)
{
  /* local block i/o variables */
  real_T rtb_LMSFilter_o2;
  real_T rtb_LMSFilter_o3[32];
  real_T rtb_SineWave;
  real_T rtb_Sum;
  real_T rtb_Sum1;
  int32_T j;
  dsp_FIRFilter_0_ex_codegen_ds_T *obj;
  real_T acc2;
  real_T rtb_LowpassFilter;
  real_T rtb_RandomSource;

  /* S-Function (sdsprandsrc2): '<Root>/Random Source' */
  RandSrc_GZ_D(&rtb_RandomSource, &ex_codegen_dsp_P.RandomSource_MeanVal, 1,
               &ex_codegen_dsp_P.RandomSource_VarianceRTP, 1,
               ex_codegen_dsp_DW.RandomSource_STATE_DWORK, 1, 1);

  /* MATLABSystem: '<Root>/Lowpass Filter' */
  if (ex_codegen_dsp_DW.obj.FilterObj->isInitialized != 1) {
    ex_codegen_dsp_DW.obj.FilterObj->isSetupComplete = false;
    ex_codegen_dsp_DW.obj.FilterObj->isInitialized = 1;
    ex_codegen_dsp_DW.obj.FilterObj->isSetupComplete = true;

    /* System object Initialization function: dsp.FIRFilter */
    for (j = 0; j < 38; j++) {
      ex_codegen_dsp_DW.obj.FilterObj->cSFunObject.W0_states[j] =
        ex_codegen_dsp_DW.obj.FilterObj->cSFunObject.P0_InitialStates;
    }
  }

  obj = &ex_codegen_dsp_DW.obj.FilterObj->cSFunObject;

  /* System object Outputs function: dsp.FIRFilter */
  /* Consume delay line and beginning of input samples */
  acc2 = rtb_RandomSource *
    ex_codegen_dsp_DW.obj.FilterObj->cSFunObject.P1_Coefficients[0];
  rtb_LowpassFilter = acc2;
  for (j = 0; j < 38; j++) {
    acc2 = obj->P1_Coefficients[1 + j] * obj->W0_states[j];
    rtb_LowpassFilter += acc2;
  }

  /* Update delay line for next frame */
  for (j = 36; j >= 0; j--) {
    obj->W0_states[1 + j] = obj->W0_states[j];
  }

  ex_codegen_dsp_DW.obj.FilterObj->cSFunObject.W0_states[0] = rtb_RandomSource;

  /* End of MATLABSystem: '<Root>/Lowpass Filter' */

  /* S-Function (sdsplms): '<Root>/LMS Filter' */
  MWSPCGlmsnw_D(&rtb_RandomSource, &rtb_LowpassFilter,
                ex_codegen_dsp_P.LMSFilter_mu,
                &ex_codegen_dsp_DW.LMSFilter_BUFF_IDX_DWORK,
                &ex_codegen_dsp_DW.LMSFilter_IN_BUFFER_DWORK[0U],
                &ex_codegen_dsp_DW.LMSFilter_WGT_IC_DWORK[0U], 32,
                ex_codegen_dsp_P.LMSFilter_leakage, 1, &rtb_Sum1,
                &rtb_LMSFilter_o2, &rtb_LMSFilter_o3[0U]);

  /* SignalToWorkspace: '<Root>/Signal To Workspace' */
  rt_UpdateLogVar((LogVar *)(LogVar*)
                  (ex_codegen_dsp_DW.SignalToWorkspace_PWORK.LoggedData),
                  &rtb_LMSFilter_o3[0], 0);

  /* S-Function (sdspsine2): '<Root>/Sine Wave' */
  rtb_SineWave = ex_codegen_dsp_P.SineWave_Amplitude * sin
    (ex_codegen_dsp_DW.SineWave_AccFreqNorm);

  /* Update accumulated normalized freq value
     for next sample.  Keep in range [0 2*pi) */
  ex_codegen_dsp_DW.SineWave_AccFreqNorm += ex_codegen_dsp_P.SineWave_Frequency *
    0.31415926535897931;
  if (ex_codegen_dsp_DW.SineWave_AccFreqNorm >= 6.2831853071795862) {
    ex_codegen_dsp_DW.SineWave_AccFreqNorm -= 6.2831853071795862;
  } else {
    if (ex_codegen_dsp_DW.SineWave_AccFreqNorm < 0.0) {
      ex_codegen_dsp_DW.SineWave_AccFreqNorm += 6.2831853071795862;
    }
  }

  /* End of S-Function (sdspsine2): '<Root>/Sine Wave' */

  /* Sum: '<Root>/Sum' */
  rtb_Sum = rtb_SineWave + rtb_LowpassFilter;

  /* Sum: '<Root>/Sum1' */
  rtb_Sum1 = rtb_Sum - rtb_Sum1;
}

/* Model update function */
static void ex_codegen_dsp_update(void)
{
  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++ex_codegen_dsp_M->Timing.clockTick0)) {
    ++ex_codegen_dsp_M->Timing.clockTickH0;
  }

  ex_codegen_dsp_M->Timing.t[0] = ex_codegen_dsp_M->Timing.clockTick0 *
    ex_codegen_dsp_M->Timing.stepSize0 + ex_codegen_dsp_M->Timing.clockTickH0 *
    ex_codegen_dsp_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
static void ex_codegen_dsp_initialize(void)
{
  {
    real_T arg;
    dspcodegen_FIRFilter_ex_codeg_T *iobj_0;
    static const real_T tmp[39] = { -0.00059194361517663689,
      -0.0020231187609940074, -0.0021060229913033154, 0.0012560050342482782,
      0.0045836073548239592, 0.00081354750291903971, -0.0074945947911449614,
      -0.0059854409238031476, 0.0089175430449157777, 0.014763032770771055,
      -0.0061136921890367051, -0.026666288251129273, -0.00451612337795531,
      0.040069256336854729, 0.028624786731400381, -0.052515463472640085,
      -0.082387824369927212, 0.061352556551937866, 0.30991387818783289,
      0.43544828894513982, 0.30991387818783289, 0.061352556551937866,
      -0.082387824369927212, -0.052515463472640085, 0.028624786731400381,
      0.040069256336854729, -0.00451612337795531, -0.026666288251129273,
      -0.0061136921890367051, 0.014763032770771055, 0.0089175430449157777,
      -0.0059854409238031476, -0.0074945947911449614, 0.00081354750291903971,
      0.0045836073548239592, 0.0012560050342482782, -0.0021060229913033154,
      -0.0020231187609940074, -0.00059194361517663689 };

    int32_T i;

    /* SetupRuntimeResources for SignalToWorkspace: '<Root>/Signal To Workspace' */
    {
      int_T dimensions[2] = { 32, 1 };

      ex_codegen_dsp_DW.SignalToWorkspace_PWORK.LoggedData = rt_CreateLogVar(
        ex_codegen_dsp_M->rtwLogInfo,
        0.0,
        rtmGetTFinal(ex_codegen_dsp_M),
        ex_codegen_dsp_M->Timing.stepSize0,
        (&rtmGetErrorStatus(ex_codegen_dsp_M)),
        "filter_wts",
        SS_DOUBLE,
        0,
        0,
        0,
        32,
        2,
        dimensions,
        NO_LOGVALDIMS,
        (NULL),
        (NULL),
        0,
        1,
        0.05,
        1);
      if (ex_codegen_dsp_DW.SignalToWorkspace_PWORK.LoggedData == (NULL))
        return;
    }

    /* Start for MATLABSystem: '<Root>/Lowpass Filter' */
    ex_codegen_dsp_DW.gobj_1.matlabCodegenIsDeleted = true;
    ex_codegen_dsp_DW.gobj_0.matlabCodegenIsDeleted = true;
    ex_codegen_dsp_DW.obj.matlabCodegenIsDeleted = true;
    ex_codegen_dsp_DW.obj.isInitialized = 0;
    ex_codegen_dsp_DW.obj.NumChannels = -1;
    ex_codegen_dsp_DW.obj.matlabCodegenIsDeleted = false;
    ex_codegen_dsp_DW.objisempty = true;
    iobj_0 = &ex_codegen_dsp_DW.gobj_0;
    ex_codegen_dsp_DW.obj.isSetupComplete = false;
    ex_codegen_dsp_DW.obj.isInitialized = 1;
    ex_codegen_dsp_DW.gobj_0.isInitialized = 0;

    /* System object Constructor function: dsp.FIRFilter */
    iobj_0->cSFunObject.P0_InitialStates = 0.0;
    for (i = 0; i < 39; i++) {
      iobj_0->cSFunObject.P1_Coefficients[i] = tmp[i];
    }

    ex_codegen_dsp_DW.gobj_0.matlabCodegenIsDeleted = false;
    ex_codegen_dsp_DW.obj.FilterObj = &ex_codegen_dsp_DW.gobj_0;
    ex_codegen_dsp_DW.obj.NumChannels = 1;
    ex_codegen_dsp_DW.obj.isSetupComplete = true;

    /* End of Start for MATLABSystem: '<Root>/Lowpass Filter' */
    /* Start for S-Function (sdspsine2): '<Root>/Sine Wave' */
    /* Trigonometric mode: compute accumulated
       normalized trig fcn argument for each channel */
    /* Keep normalized value in range [0 2*pi) */
    arg = fmod(ex_codegen_dsp_P.SineWave_Phase, 6.2831853071795862);
    if (arg < 0.0) {
      arg += 6.2831853071795862;
    }

    ex_codegen_dsp_DW.SineWave_AccFreqNorm = arg;

    /* End of Start for S-Function (sdspsine2): '<Root>/Sine Wave' */
  }

  {
    real_T arg;
    int32_T i;

    /* InitializeConditions for S-Function (sdsprandsrc2): '<Root>/Random Source' */
    ex_codegen_dsp_DW.RandomSource_SEED_DWORK =
      ex_codegen_dsp_P.RandomSource_rawSeed;
    RandSrcInitState_GZ(&ex_codegen_dsp_DW.RandomSource_SEED_DWORK,
                        ex_codegen_dsp_DW.RandomSource_STATE_DWORK, 1);

    /* InitializeConditions for S-Function (sdsplms): '<Root>/LMS Filter' */
    ex_codegen_dsp_DW.LMSFilter_BUFF_IDX_DWORK = 0U;
    memset(&ex_codegen_dsp_DW.LMSFilter_WGT_IC_DWORK[0], 0, sizeof(real_T) << 5U);
    memset(&ex_codegen_dsp_DW.LMSFilter_IN_BUFFER_DWORK[0], 0, sizeof(real_T) <<
           5U);

    /* InitializeConditions for S-Function (sdspsine2): '<Root>/Sine Wave' */
    /* This code only executes when block is re-enabled in an
       enabled subsystem when the enabled subsystem states on
       re-enabling are set to 'Reset' */
    /* Reset to time zero on re-enable */
    /* Trigonometric mode: compute accumulated
       normalized trig fcn argument for each channel */
    /* Keep normalized value in range [0 2*pi) */
    arg = fmod(ex_codegen_dsp_P.SineWave_Phase, 6.2831853071795862);
    if (arg < 0.0) {
      arg += 6.2831853071795862;
    }

    ex_codegen_dsp_DW.SineWave_AccFreqNorm = arg;

    /* End of InitializeConditions for S-Function (sdspsine2): '<Root>/Sine Wave' */

    /* InitializeConditions for MATLABSystem: '<Root>/Lowpass Filter' */
    if (ex_codegen_dsp_DW.obj.FilterObj->isInitialized == 1) {
      /* System object Initialization function: dsp.FIRFilter */
      for (i = 0; i < 38; i++) {
        ex_codegen_dsp_DW.obj.FilterObj->cSFunObject.W0_states[i] =
          ex_codegen_dsp_DW.obj.FilterObj->cSFunObject.P0_InitialStates;
      }
    }

    /* End of InitializeConditions for MATLABSystem: '<Root>/Lowpass Filter' */
  }
}

/* Model terminate function */
static void ex_codegen_dsp_terminate(void)
{
  /* Terminate for MATLABSystem: '<Root>/Lowpass Filter' */
  matlabCodegenHandle_matlabCod_l(&ex_codegen_dsp_DW.obj);
  matlabCodegenHandle_matlabCodeg(&ex_codegen_dsp_DW.gobj_0);
  matlabCodegenHandle_matlabCodeg(&ex_codegen_dsp_DW.gobj_1);
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  ex_codegen_dsp_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  ex_codegen_dsp_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  ex_codegen_dsp_initialize();
}

void MdlTerminate(void)
{
  ex_codegen_dsp_terminate();
}

/* Registration function */
RT_MODEL_ex_codegen_dsp_T *ex_codegen_dsp(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)ex_codegen_dsp_M, 0,
                sizeof(RT_MODEL_ex_codegen_dsp_T));

  /* Initialize timing info */
  {
    int_T *mdlTsMap = ex_codegen_dsp_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    ex_codegen_dsp_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    ex_codegen_dsp_M->Timing.sampleTimes =
      (&ex_codegen_dsp_M->Timing.sampleTimesArray[0]);
    ex_codegen_dsp_M->Timing.offsetTimes =
      (&ex_codegen_dsp_M->Timing.offsetTimesArray[0]);

    /* task periods */
    ex_codegen_dsp_M->Timing.sampleTimes[0] = (0.05);

    /* task offsets */
    ex_codegen_dsp_M->Timing.offsetTimes[0] = (0.0);
  }

  rtmSetTPtr(ex_codegen_dsp_M, &ex_codegen_dsp_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = ex_codegen_dsp_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    ex_codegen_dsp_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(ex_codegen_dsp_M, 60.0);
  ex_codegen_dsp_M->Timing.stepSize0 = 0.05;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    ex_codegen_dsp_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(ex_codegen_dsp_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(ex_codegen_dsp_M->rtwLogInfo, (NULL));
    rtliSetLogT(ex_codegen_dsp_M->rtwLogInfo, "tout");
    rtliSetLogX(ex_codegen_dsp_M->rtwLogInfo, "");
    rtliSetLogXFinal(ex_codegen_dsp_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(ex_codegen_dsp_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(ex_codegen_dsp_M->rtwLogInfo, 0);
    rtliSetLogMaxRows(ex_codegen_dsp_M->rtwLogInfo, 1000);
    rtliSetLogDecimation(ex_codegen_dsp_M->rtwLogInfo, 1);
    rtliSetLogY(ex_codegen_dsp_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(ex_codegen_dsp_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(ex_codegen_dsp_M->rtwLogInfo, (NULL));
  }

  ex_codegen_dsp_M->solverInfoPtr = (&ex_codegen_dsp_M->solverInfo);
  ex_codegen_dsp_M->Timing.stepSize = (0.05);
  rtsiSetFixedStepSize(&ex_codegen_dsp_M->solverInfo, 0.05);
  rtsiSetSolverMode(&ex_codegen_dsp_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* parameters */
  ex_codegen_dsp_M->defaultParam = ((real_T *)&ex_codegen_dsp_P);

  /* states (dwork) */
  ex_codegen_dsp_M->dwork = ((void *) &ex_codegen_dsp_DW);
  (void) memset((void *)&ex_codegen_dsp_DW, 0,
                sizeof(DW_ex_codegen_dsp_T));

  /* Initialize Sizes */
  ex_codegen_dsp_M->Sizes.numContStates = (0);/* Number of continuous states */
  ex_codegen_dsp_M->Sizes.numY = (0);  /* Number of model outputs */
  ex_codegen_dsp_M->Sizes.numU = (0);  /* Number of model inputs */
  ex_codegen_dsp_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  ex_codegen_dsp_M->Sizes.numSampTimes = (1);/* Number of sample times */
  ex_codegen_dsp_M->Sizes.numBlocks = (10);/* Number of blocks */
  ex_codegen_dsp_M->Sizes.numBlockIO = (0);/* Number of block outputs */
  ex_codegen_dsp_M->Sizes.numBlockPrms = (8);/* Sum of parameter "widths" */
  return ex_codegen_dsp_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
