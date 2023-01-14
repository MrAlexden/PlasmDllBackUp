/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve2.h
 *
 * Contents:  Declares lmcurve2(), a variant of lmcurve() that weighs
 *            data points y(t) with the inverse of the standard deviations dy.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  https://jugit.fz-juelich.de/mlz/lmfit
 */

#ifndef LMCURVETYD_H
#define LMCURVETYD_H

#include "lmstruct.h"

__BEGIN_DECLS

LM_DLL void lmcurve2(
    const int n_par, double* par, double* parerr, double* covar,
    const int m_dat, const double* t, const double* y, const double* dy,
    double (*f)(double t, const double* par),
    const lm_control_struct* control, lm_status_struct* status);

__END_DECLS
#endif /* LMCURVETYD_H */
