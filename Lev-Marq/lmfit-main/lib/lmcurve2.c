/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve2.c
 *
 * Contents:  Implements lmcurve2(), a variant of lmcurve() that weighs
 *            data points y(t) with the inverse of the standard deviations dy.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  https://jugit.fz-juelich.de/mlz/lmfit
 */

#include "lmcurve2.h"
#include <math.h>
#include "lmmin.h"

typedef struct {
    const double* t;
    const double* y;
    const double* dy;
    double (*f)(const double t, const double* par);
} lmcurve2_data_struct;

void lmcurve2_evaluate(
    const double* par, const int m_dat, const void* data, double* fvec,
    int* info)
{
    lmcurve2_data_struct* D = (lmcurve2_data_struct*)data;
    int i;
    for (i = 0; i < m_dat; i++)
        fvec[i] = ( D->y[i] - D->f(D->t[i], par) ) / D->dy[i];
}

void lmcurve2(
    const int n_par, double* par, double* parerr, double* covar,
    const int m_dat, const double* t, const double* y, const double* dy,
    double (*f)(const double t, const double* par),
    const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve2_data_struct data = { t, y, dy, f };

    lmmin2(n_par, par, NULL, covar, m_dat, NULL, (const void*)&data, lmcurve2_evaluate,
          control, status);
    int j;
    if (parerr)
        for (j = 0; j < n_par; j++)
            parerr[j] = sqrt(covar[j*n_par+j]);
}
