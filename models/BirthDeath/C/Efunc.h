/*
 * Efunc.h
 *
 * Code generation for function 'Efunc'
 *
 */

#ifndef __EFUNC_H__
#define __EFUNC_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "reactionFlux_types.h"

/* Function Declarations */
extern void Efunc(const double in1[2], double t, const double in3[6], double E[8]);

#endif

/* End of code generation (Efunc.h) */
