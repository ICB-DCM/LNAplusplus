/*
 * _coder_reactionFlux_api.h
 *
 * Code generation for function 'reactionFlux'
 *
 */

#ifndef ___CODER_REACTIONFLUX_API_H__
#define ___CODER_REACTIONFLUX_API_H__
/* Include files */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Function Declarations */
extern void reactionFlux_initialize(emlrtContext *aContext);
extern void reactionFlux_terminate(void);
extern void reactionFlux_atexit(void);
extern void reactionFlux_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void reactionFlux(real_T in1[2], real_T t, real_T in3[6], real_T F[4]);
extern void J_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void J(real_T in1[2], real_T t, real_T in3[6], real_T b_J[8]);
extern void dFdTheta_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void dFdTheta(real_T in1[2], real_T t, real_T in3[6], real_T b_dFdTheta[24]);
extern void d2fdTheta2_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2fdTheta2(real_T in1[2], real_T t, real_T in3[6], real_T b_d2fdTheta2[144]);
extern void Afunc_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void Afunc(real_T in1[2], real_T t, real_T in3[6], real_T A[4]);
extern void dAdTheta_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void dAdTheta(real_T in1[2], real_T t, real_T in3[6], real_T b_dAdTheta[24]);
extern void dAdPhi_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void dAdPhi(real_T in1[2], real_T t, real_T in3[6], real_T b_dAdPhi[8]);
extern void d2AdPhi2_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2AdPhi2(real_T in1[2], real_T t, real_T in3[6], real_T b_d2AdPhi2[16]);
extern void d2AdTheta2_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2AdTheta2(real_T in1[2], real_T t, real_T in3[6], real_T b_d2AdTheta2[144]);
extern void d2AdThetadPhi_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2AdThetadPhi(real_T in1[2], real_T t, real_T in3[6], real_T b_d2AdThetadPhi[48]);
extern void d2AdPhidTheta_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2AdPhidTheta(real_T in1[2], real_T t, real_T in3[6], real_T b_d2AdPhidTheta[48]);
extern void Efunc_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void Efunc(real_T in1[2], real_T t, real_T in3[6], real_T E[8]);
extern void dEdTheta_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void dEdTheta(real_T in1[2], real_T t, real_T in3[6], real_T b_dEdTheta[48]);
extern void d2EdTheta2_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2EdTheta2(real_T in1[2], real_T t, real_T in3[6], real_T b_d2EdTheta2[288]);
extern void d2EdThetadPhi_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2EdThetadPhi(real_T in1[2], real_T t, real_T in3[6], real_T b_d2EdThetadPhi[96]);
extern void d2EdPhidTheta_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2EdPhidTheta(real_T in1[2], real_T t, real_T in3[6], real_T b_d2EdPhidTheta[96]);
extern void dEdPhi_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void dEdPhi(real_T in1[2], real_T t, real_T in3[6], real_T b_dEdPhi[16]);
extern void d2EdPhi2_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void d2EdPhi2(real_T in1[2], real_T t, real_T in3[6], real_T b_d2EdPhi2[32]);
extern void MI_api(const mxArray *prhs[4], const mxArray *plhs[1]);
extern void MI(real_T in1[2], real_T t, real_T in3[6], real_T b_gamma, real_T b_MI[81]);
extern void S0_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void S0(real_T in1[6], real_T out1[12]);
extern void S20_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void S20(real_T in1[6], real_T out1[72]);
extern void SV0_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void SV0(real_T in1[6], real_T out1[18]);
extern void S2V0_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void S2V0(real_T in1[6], real_T out1[108]);
extern void Y0_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void Y0(real_T in1[6], real_T out1[2]);
extern void V0_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void V0(real_T in1[6], real_T out1[3]);
extern void systemJacobian_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void systemJacobian(real_T in1[2], real_T t, real_T in3[6], real_T b_systemJacobian[81]);
extern void reactionFlux_xil_terminate(void);

#endif
/* End of code generation (_coder_reactionFlux_api.h) */
