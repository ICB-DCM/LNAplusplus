/*
 * systemJacobian.h
 *
 * Code generation for function 'systemJacobian'
 *
 */

#ifndef __SYSTEMJACOBIAN_H__
#define __SYSTEMJACOBIAN_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "reactionFlux_types.h"

/* Function Declarations */
extern void systemJacobian(const double in1[2], double t, const double in3[6],
  double b_systemJacobian[81]);

#endif

/* End of code generation (systemJacobian.h) */
