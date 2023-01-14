/*
 * Library:  lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:     democpp/nonlin1.cpp
 *
 * Contents: Example for the solution of 2 nonlinear equations in 2 variables.
 *           Find the intersection of a circle and a parabola.
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

void evaluate_nonlin1(
    const double *p, int n, const void *data, double *f, int *info)
{
    f[0] = p[0]*p[0] + p[1]*p[1] - 1;
    f[1] = p[1] - p[0]*p[0];
}

int main( int argc, char **argv)
{
    std::vector<double> par0;

    lm_control_struct control = lm_control_double;
    control.verbosity = 31;

    /* get start values from command line */
    if (argc != 3 ){
        std::cout << "usage: nonlin1cpp x_start y_start" << '\n';
        exit(-1);
    }
    par0.push_back( atof(argv[1]) );
    par0.push_back( atof(argv[2]) );

    std::cout << "Minimization:" << '\n';
    auto result = lmfit::minimize(par0, nullptr, 2, &evaluate_nonlin1, control);

    std::cout << "\nlmmin status after " << result.status.nfev
              << " function evaluations:" << '\n';
    std::cout << lm_infmsg[result.status.outcome] << '\n';

    std::cout << '\n' << "Solution:\n";
    std::cout << "  x = " << result.par[0] << '\n';
    std::cout << "  y = " << result.par[1] << '\n';
    std::cout << "  d = " << result.status.fnorm << " => ";
    if (result.status.fnorm >= control.ftol)
        std::cout << "not a valid solution, try other starting values" << '\n';
    else
        std::cout << "valid, though not the only solution:"
                  << " try other starting values" << '\n';

    return 0;
}
