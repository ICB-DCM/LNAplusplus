/*
 * _coder_reactionFlux_api.c
 *
 * Code generation for function 'reactionFlux'
 *
 */

/* Include files */
#include "_coder_reactionFlux_api.h"

/* Function Declarations */
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *in1, const
  char_T *identifier))[2];
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2];
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *t, const
  char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *in3,
  const char_T *identifier))[6];
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[6];
static const mxArray *emlrt_marshallOut(const real_T u[4]);
static const mxArray *b_emlrt_marshallOut(const real_T u[8]);
static const mxArray *c_emlrt_marshallOut(const real_T u[24]);
static const mxArray *d_emlrt_marshallOut(const real_T u[144]);
static const mxArray *e_emlrt_marshallOut(const real_T u[4]);
static const mxArray *f_emlrt_marshallOut(const real_T u[16]);
static const mxArray *g_emlrt_marshallOut(const real_T u[48]);
static const mxArray *h_emlrt_marshallOut(const real_T u[48]);
static const mxArray *i_emlrt_marshallOut(const real_T u[8]);
static const mxArray *j_emlrt_marshallOut(const real_T u[288]);
static const mxArray *k_emlrt_marshallOut(const real_T u[96]);
static const mxArray *l_emlrt_marshallOut(const real_T u[96]);
static const mxArray *m_emlrt_marshallOut(const real_T u[32]);
static const mxArray *n_emlrt_marshallOut(const real_T u[81]);
static const mxArray *o_emlrt_marshallOut(const real_T u[12]);
static const mxArray *p_emlrt_marshallOut(const real_T u[72]);
static const mxArray *q_emlrt_marshallOut(const real_T u[18]);
static const mxArray *r_emlrt_marshallOut(const real_T u[108]);
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2];
static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6];

/* Function Definitions */
void reactionFlux_initialize(emlrtContext *aContext)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void reactionFlux_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void reactionFlux_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  reactionFlux_xil_terminate();
}

void reactionFlux_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*F)[4];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  F = (real_T (*)[4])mxMalloc(sizeof(real_T [4]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  reactionFlux(*in1, t, *in3, *F);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*F);
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *in1, const
  char_T *identifier))[2]
{
  real_T (*y)[2];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(sp, emlrtAlias(in1), &thisId);
  emlrtDestroyArray(&in1);
  return y;
}
  static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[2]
{
  real_T (*y)[2];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *t, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(sp, emlrtAlias(t), &thisId);
  emlrtDestroyArray(&t);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *in3,
  const char_T *identifier))[6]
{
  real_T (*y)[6];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = f_emlrt_marshallIn(sp, emlrtAlias(in3), &thisId);
  emlrtDestroyArray(&in3);
  return y;
}
  static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[6]
{
  real_T (*y)[6];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u[4])
{
  const mxArray *y;
  static const int32_T iv0[2] = { 0, 0 };

  const mxArray *m0;
  static const int32_T iv1[2] = { 1, 4 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u);
  emlrtSetDimensions((mxArray *)m0, iv1, 2);
  emlrtAssign(&y, m0);
  return y;
}

void J_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_J)[8];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_J = (real_T (*)[8])mxMalloc(sizeof(real_T [8]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  J(*in1, t, *in3, *b_J);

  /* Marshall function outputs */
  plhs[0] = b_emlrt_marshallOut(*b_J);
}

static const mxArray *b_emlrt_marshallOut(const real_T u[8])
{
  const mxArray *y;
  static const int32_T iv2[2] = { 0, 0 };

  const mxArray *m1;
  static const int32_T iv3[2] = { 4, 2 };

  y = NULL;
  m1 = emlrtCreateNumericArray(2, iv2, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m1, (void *)u);
  emlrtSetDimensions((mxArray *)m1, iv3, 2);
  emlrtAssign(&y, m1);
  return y;
}

void dFdTheta_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_dFdTheta)[24];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_dFdTheta = (real_T (*)[24])mxMalloc(sizeof(real_T [24]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  dFdTheta(*in1, t, *in3, *b_dFdTheta);

  /* Marshall function outputs */
  plhs[0] = c_emlrt_marshallOut(*b_dFdTheta);
}

static const mxArray *c_emlrt_marshallOut(const real_T u[24])
{
  const mxArray *y;
  static const int32_T iv4[2] = { 0, 0 };

  const mxArray *m2;
  static const int32_T iv5[2] = { 4, 6 };

  y = NULL;
  m2 = emlrtCreateNumericArray(2, iv4, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m2, (void *)u);
  emlrtSetDimensions((mxArray *)m2, iv5, 2);
  emlrtAssign(&y, m2);
  return y;
}

void d2fdTheta2_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2fdTheta2)[144];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2fdTheta2 = (real_T (*)[144])mxMalloc(sizeof(real_T [144]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2fdTheta2(*in1, t, *in3, *b_d2fdTheta2);

  /* Marshall function outputs */
  plhs[0] = d_emlrt_marshallOut(*b_d2fdTheta2);
}

static const mxArray *d_emlrt_marshallOut(const real_T u[144])
{
  const mxArray *y;
  static const int32_T iv6[2] = { 0, 0 };

  const mxArray *m3;
  static const int32_T iv7[2] = { 24, 6 };

  y = NULL;
  m3 = emlrtCreateNumericArray(2, iv6, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m3, (void *)u);
  emlrtSetDimensions((mxArray *)m3, iv7, 2);
  emlrtAssign(&y, m3);
  return y;
}

void Afunc_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*A)[4];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  A = (real_T (*)[4])mxMalloc(sizeof(real_T [4]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  Afunc(*in1, t, *in3, *A);

  /* Marshall function outputs */
  plhs[0] = e_emlrt_marshallOut(*A);
}

static const mxArray *e_emlrt_marshallOut(const real_T u[4])
{
  const mxArray *y;
  static const int32_T iv8[2] = { 0, 0 };

  const mxArray *m4;
  static const int32_T iv9[2] = { 2, 2 };

  y = NULL;
  m4 = emlrtCreateNumericArray(2, iv8, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m4, (void *)u);
  emlrtSetDimensions((mxArray *)m4, iv9, 2);
  emlrtAssign(&y, m4);
  return y;
}

void dAdTheta_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_dAdTheta)[24];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_dAdTheta = (real_T (*)[24])mxMalloc(sizeof(real_T [24]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  dAdTheta(*in1, t, *in3, *b_dAdTheta);

  /* Marshall function outputs */
  plhs[0] = c_emlrt_marshallOut(*b_dAdTheta);
}

void dAdPhi_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_dAdPhi)[8];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_dAdPhi = (real_T (*)[8])mxMalloc(sizeof(real_T [8]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  dAdPhi(*in1, t, *in3, *b_dAdPhi);

  /* Marshall function outputs */
  plhs[0] = b_emlrt_marshallOut(*b_dAdPhi);
}

void d2AdPhi2_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2AdPhi2)[16];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2AdPhi2 = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2AdPhi2(*in1, t, *in3, *b_d2AdPhi2);

  /* Marshall function outputs */
  plhs[0] = f_emlrt_marshallOut(*b_d2AdPhi2);
}

static const mxArray *f_emlrt_marshallOut(const real_T u[16])
{
  const mxArray *y;
  static const int32_T iv10[2] = { 0, 0 };

  const mxArray *m5;
  static const int32_T iv11[2] = { 8, 2 };

  y = NULL;
  m5 = emlrtCreateNumericArray(2, iv10, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m5, (void *)u);
  emlrtSetDimensions((mxArray *)m5, iv11, 2);
  emlrtAssign(&y, m5);
  return y;
}

void d2AdTheta2_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2AdTheta2)[144];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2AdTheta2 = (real_T (*)[144])mxMalloc(sizeof(real_T [144]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2AdTheta2(*in1, t, *in3, *b_d2AdTheta2);

  /* Marshall function outputs */
  plhs[0] = d_emlrt_marshallOut(*b_d2AdTheta2);
}

void d2AdThetadPhi_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2AdThetadPhi)[48];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2AdThetadPhi = (real_T (*)[48])mxMalloc(sizeof(real_T [48]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2AdThetadPhi(*in1, t, *in3, *b_d2AdThetadPhi);

  /* Marshall function outputs */
  plhs[0] = g_emlrt_marshallOut(*b_d2AdThetadPhi);
}

static const mxArray *g_emlrt_marshallOut(const real_T u[48])
{
  const mxArray *y;
  static const int32_T iv12[2] = { 0, 0 };

  const mxArray *m6;
  static const int32_T iv13[2] = { 24, 2 };

  y = NULL;
  m6 = emlrtCreateNumericArray(2, iv12, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m6, (void *)u);
  emlrtSetDimensions((mxArray *)m6, iv13, 2);
  emlrtAssign(&y, m6);
  return y;
}

void d2AdPhidTheta_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2AdPhidTheta)[48];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2AdPhidTheta = (real_T (*)[48])mxMalloc(sizeof(real_T [48]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2AdPhidTheta(*in1, t, *in3, *b_d2AdPhidTheta);

  /* Marshall function outputs */
  plhs[0] = h_emlrt_marshallOut(*b_d2AdPhidTheta);
}

static const mxArray *h_emlrt_marshallOut(const real_T u[48])
{
  const mxArray *y;
  static const int32_T iv14[2] = { 0, 0 };

  const mxArray *m7;
  static const int32_T iv15[2] = { 8, 6 };

  y = NULL;
  m7 = emlrtCreateNumericArray(2, iv14, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m7, (void *)u);
  emlrtSetDimensions((mxArray *)m7, iv15, 2);
  emlrtAssign(&y, m7);
  return y;
}

void Efunc_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*E)[8];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  E = (real_T (*)[8])mxMalloc(sizeof(real_T [8]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  Efunc(*in1, t, *in3, *E);

  /* Marshall function outputs */
  plhs[0] = i_emlrt_marshallOut(*E);
}

static const mxArray *i_emlrt_marshallOut(const real_T u[8])
{
  const mxArray *y;
  static const int32_T iv16[2] = { 0, 0 };

  const mxArray *m8;
  static const int32_T iv17[2] = { 2, 4 };

  y = NULL;
  m8 = emlrtCreateNumericArray(2, iv16, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m8, (void *)u);
  emlrtSetDimensions((mxArray *)m8, iv17, 2);
  emlrtAssign(&y, m8);
  return y;
}

void dEdTheta_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_dEdTheta)[48];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_dEdTheta = (real_T (*)[48])mxMalloc(sizeof(real_T [48]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  dEdTheta(*in1, t, *in3, *b_dEdTheta);

  /* Marshall function outputs */
  plhs[0] = h_emlrt_marshallOut(*b_dEdTheta);
}

void d2EdTheta2_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2EdTheta2)[288];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2EdTheta2 = (real_T (*)[288])mxMalloc(sizeof(real_T [288]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2EdTheta2(*in1, t, *in3, *b_d2EdTheta2);

  /* Marshall function outputs */
  plhs[0] = j_emlrt_marshallOut(*b_d2EdTheta2);
}

static const mxArray *j_emlrt_marshallOut(const real_T u[288])
{
  const mxArray *y;
  static const int32_T iv18[2] = { 0, 0 };

  const mxArray *m9;
  static const int32_T iv19[2] = { 48, 6 };

  y = NULL;
  m9 = emlrtCreateNumericArray(2, iv18, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m9, (void *)u);
  emlrtSetDimensions((mxArray *)m9, iv19, 2);
  emlrtAssign(&y, m9);
  return y;
}

void d2EdThetadPhi_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2EdThetadPhi)[96];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2EdThetadPhi = (real_T (*)[96])mxMalloc(sizeof(real_T [96]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2EdThetadPhi(*in1, t, *in3, *b_d2EdThetadPhi);

  /* Marshall function outputs */
  plhs[0] = k_emlrt_marshallOut(*b_d2EdThetadPhi);
}

static const mxArray *k_emlrt_marshallOut(const real_T u[96])
{
  const mxArray *y;
  static const int32_T iv20[2] = { 0, 0 };

  const mxArray *m10;
  static const int32_T iv21[2] = { 48, 2 };

  y = NULL;
  m10 = emlrtCreateNumericArray(2, iv20, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m10, (void *)u);
  emlrtSetDimensions((mxArray *)m10, iv21, 2);
  emlrtAssign(&y, m10);
  return y;
}

void d2EdPhidTheta_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2EdPhidTheta)[96];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2EdPhidTheta = (real_T (*)[96])mxMalloc(sizeof(real_T [96]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2EdPhidTheta(*in1, t, *in3, *b_d2EdPhidTheta);

  /* Marshall function outputs */
  plhs[0] = l_emlrt_marshallOut(*b_d2EdPhidTheta);
}

static const mxArray *l_emlrt_marshallOut(const real_T u[96])
{
  const mxArray *y;
  static const int32_T iv22[2] = { 0, 0 };

  const mxArray *m11;
  static const int32_T iv23[2] = { 16, 6 };

  y = NULL;
  m11 = emlrtCreateNumericArray(2, iv22, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m11, (void *)u);
  emlrtSetDimensions((mxArray *)m11, iv23, 2);
  emlrtAssign(&y, m11);
  return y;
}

void dEdPhi_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_dEdPhi)[16];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_dEdPhi = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  dEdPhi(*in1, t, *in3, *b_dEdPhi);

  /* Marshall function outputs */
  plhs[0] = f_emlrt_marshallOut(*b_dEdPhi);
}

void d2EdPhi2_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_d2EdPhi2)[32];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_d2EdPhi2 = (real_T (*)[32])mxMalloc(sizeof(real_T [32]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  d2EdPhi2(*in1, t, *in3, *b_d2EdPhi2);

  /* Marshall function outputs */
  plhs[0] = m_emlrt_marshallOut(*b_d2EdPhi2);
}

static const mxArray *m_emlrt_marshallOut(const real_T u[32])
{
  const mxArray *y;
  static const int32_T iv24[2] = { 0, 0 };

  const mxArray *m12;
  static const int32_T iv25[2] = { 16, 2 };

  y = NULL;
  m12 = emlrtCreateNumericArray(2, iv24, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m12, (void *)u);
  emlrtSetDimensions((mxArray *)m12, iv25, 2);
  emlrtAssign(&y, m12);
  return y;
}

void MI_api(const mxArray *prhs[4], const mxArray *plhs[1])
{
  real_T (*b_MI)[81];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  real_T b_gamma;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_MI = (real_T (*)[81])mxMalloc(sizeof(real_T [81]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");
  b_gamma = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "gamma");

  /* Invoke the target function */
  MI(*in1, t, *in3, b_gamma, *b_MI);

  /* Marshall function outputs */
  plhs[0] = n_emlrt_marshallOut(*b_MI);
}

static const mxArray *n_emlrt_marshallOut(const real_T u[81])
{
  const mxArray *y;
  static const int32_T iv26[2] = { 0, 0 };

  const mxArray *m13;
  static const int32_T iv27[2] = { 9, 9 };

  y = NULL;
  m13 = emlrtCreateNumericArray(2, iv26, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m13, (void *)u);
  emlrtSetDimensions((mxArray *)m13, iv27, 2);
  emlrtAssign(&y, m13);
  return y;
}

void S0_api(const mxArray *prhs[1], const mxArray *plhs[1])
{
  real_T (*out1)[12];
  real_T (*in1)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  out1 = (real_T (*)[12])mxMalloc(sizeof(real_T [12]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);

  /* Marshall function inputs */
  in1 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");

  /* Invoke the target function */
  S0(*in1, *out1);

  /* Marshall function outputs */
  plhs[0] = o_emlrt_marshallOut(*out1);
}

static const mxArray *o_emlrt_marshallOut(const real_T u[12])
{
  const mxArray *y;
  static const int32_T iv28[2] = { 0, 0 };

  const mxArray *m14;
  static const int32_T iv29[2] = { 2, 6 };

  y = NULL;
  m14 = emlrtCreateNumericArray(2, iv28, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m14, (void *)u);
  emlrtSetDimensions((mxArray *)m14, iv29, 2);
  emlrtAssign(&y, m14);
  return y;
}

void S20_api(const mxArray *prhs[1], const mxArray *plhs[1])
{
  real_T (*out1)[72];
  real_T (*in1)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  out1 = (real_T (*)[72])mxMalloc(sizeof(real_T [72]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);

  /* Marshall function inputs */
  in1 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");

  /* Invoke the target function */
  S20(*in1, *out1);

  /* Marshall function outputs */
  plhs[0] = p_emlrt_marshallOut(*out1);
}

static const mxArray *p_emlrt_marshallOut(const real_T u[72])
{
  const mxArray *y;
  static const int32_T iv30[2] = { 0, 0 };

  const mxArray *m15;
  static const int32_T iv31[2] = { 12, 6 };

  y = NULL;
  m15 = emlrtCreateNumericArray(2, iv30, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m15, (void *)u);
  emlrtSetDimensions((mxArray *)m15, iv31, 2);
  emlrtAssign(&y, m15);
  return y;
}

void SV0_api(const mxArray *prhs[1], const mxArray *plhs[1])
{
  real_T (*out1)[18];
  real_T (*in1)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  out1 = (real_T (*)[18])mxMalloc(sizeof(real_T [18]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);

  /* Marshall function inputs */
  in1 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");

  /* Invoke the target function */
  SV0(*in1, *out1);

  /* Marshall function outputs */
  plhs[0] = q_emlrt_marshallOut(*out1);
}

static const mxArray *q_emlrt_marshallOut(const real_T u[18])
{
  const mxArray *y;
  static const int32_T iv32[2] = { 0, 0 };

  const mxArray *m16;
  static const int32_T iv33[2] = { 3, 6 };

  y = NULL;
  m16 = emlrtCreateNumericArray(2, iv32, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m16, (void *)u);
  emlrtSetDimensions((mxArray *)m16, iv33, 2);
  emlrtAssign(&y, m16);
  return y;
}

void S2V0_api(const mxArray *prhs[1], const mxArray *plhs[1])
{
  real_T (*out1)[108];
  real_T (*in1)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  out1 = (real_T (*)[108])mxMalloc(sizeof(real_T [108]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);

  /* Marshall function inputs */
  in1 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");

  /* Invoke the target function */
  S2V0(*in1, *out1);

  /* Marshall function outputs */
  plhs[0] = r_emlrt_marshallOut(*out1);
}

static const mxArray *r_emlrt_marshallOut(const real_T u[108])
{
  const mxArray *y;
  static const int32_T iv34[2] = { 0, 0 };

  const mxArray *m17;
  static const int32_T iv35[2] = { 18, 6 };

  y = NULL;
  m17 = emlrtCreateNumericArray(2, iv34, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m17, (void *)u);
  emlrtSetDimensions((mxArray *)m17, iv35, 2);
  emlrtAssign(&y, m17);
  return y;
}

void systemJacobian_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*b_systemJacobian)[81];
  real_T (*in1)[2];
  real_T t;
  real_T (*in3)[6];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  b_systemJacobian = (real_T (*)[81])mxMalloc(sizeof(real_T [81]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  in1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "in1");
  t = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "t");
  in3 = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "in3");

  /* Invoke the target function */
  systemJacobian(*in1, t, *in3, *b_systemJacobian);

  /* Marshall function outputs */
  plhs[0] = n_emlrt_marshallOut(*b_systemJacobian);
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2]
{
  real_T (*ret)[2];
  int32_T iv36[2];
  int32_T i0;
  for (i0 = 0; i0 < 2; i0++) {
    iv36[i0] = 1 + i0;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv36);
  ret = (real_T (*)[2])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6]
{
  real_T (*ret)[6];
  int32_T iv37[2];
  int32_T i1;
  for (i1 = 0; i1 < 2; i1++) {
    iv37[i1] = 1 + 5 * i1;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv37);
  ret = (real_T (*)[6])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/* End of code generation (_coder_reactionFlux_api.c) */
