/*
 * Library:  lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:     democpp/curve1.cpp
 *
 * Contents: Example for curve fitting with lmcurve():
 *           fit a data set y(x) by a curve f(x;p).
 *
 * Copyright: Janike Katter, Joachim Wuttke
 *            Forschungszentrum Juelich GmbH (2021)
 *
 * Licence:  see ../COPYING (FreeBSD)
 *
 * Homepage: https://jugit.fz-juelich.de/mlz/lmfit
 */

#include "lmfit.hpp"
#include <iostream>

double f(const double t, const double* p)
{
    return p[0] + p[1]*t + p[2]*t*t;
}

int main()
{
    std::vector<double> par0{ 100, 0, -10 };

    std::vector<double> t{ -4., -3., -2., -1.,  0., 1.,  2.,  3.,  4. };
    std::vector<double> y{ 16.6, 9.9, 4.4, 1.1, 0., 1.1, 4.2, 9.3, 16.4 };

    lm_control_struct control = lm_control_double;
    control.verbosity = 9;

    std::cout << "Fitting..." << '\n';
    auto result = lmfit::fit_curve(par0, t, y, &f, control);

    std::cout << "Results:" << '\n';
    std::cout << "status after " << result.status.nfev
              << " function evaluations:" << '\n';
    std::cout << lm_infmsg[result.status.outcome] << '\n';

    std::cout << "Obtained parameters:" << '\n';
    for (size_t j = 0; j < result.par.size(); ++j)
        std::cout << "par[" << j << "] = " << result.par[j] << '\n';
    std::cout << "Obtained norm:" << '\n';
    std::cout << result.status.fnorm << '\n';

    std::cout << "fitting data as follows:" << '\n';
    for (size_t i = 0; i < t.size(); ++i)
        std::cout << "t[" << i << "]= " << t[i] << '\n';

    if (result.status.outcome > 3) {
        std::cout << "FAILURE" << '\n';
        return 1;
    }
    std::cout << "SUCCESS" << '\n';
    return 0;
}
