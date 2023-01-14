/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lminvert.h
 *
 * Contents:  Declares square matrix inversion function.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2018-)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  https://jugit.fz-juelich.de/mlz/lmfit
 */

#include "lmdecls.h"

LM_DLL void lm_invert(double *const A, const int n, int *const P,
                      double *const ws, double *const IA, int*const failure);
