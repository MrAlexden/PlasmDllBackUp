#pragma once

#define WIN32_LEAN_AND_MEAN // Исключите редко используемые компоненты из заголовков Windows
// Файлы заголовков Windows
#include <windows.h>
#include <stdio.h>
#include <fstream>          // for file access
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <cstddef>          // for size_t
#include <algorithm>

using namespace std;

typedef float myflo;

////////////////////////////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////////
/*===================================== GAUSS =====================================*/
myflo fx_GAUSS(myflo x, vector <myflo>& vParams); // from MainDllFunc.cpp

/*===================================== STEP =====================================*/
myflo fx_STEP(myflo x, vector <myflo>& vParams); // from MainDllFunc.cpp

/*===================================== LINE =====================================*/
myflo fx_LINE(myflo, vector <myflo>&); // from Lev-Marq.cpp

/* make linear approximation of given data, returns vector(2) with A and B */
vector <myflo> linear_fit(vector <myflo>&, vector <myflo>&); // from Lev-Marq.cpp

/* multiply first and second vector, third elemnt is the result */
void vectorElementsProduct(vector <myflo>&, vector <myflo>&, vector <myflo>&); // from PeakFinder.cpp

/* multiply the given scalar on given vector */
void scalarProduct(int, vector <myflo>&); // from PeakFinder.cpp

namespace PeakFinder {
    const myflo EPS = 2.2204e-16f;

    /*
        Inputs
        x0: input signal
        extrema: 1 if maxima are desired, -1 if minima are desired
        includeEndpoints - If true the endpoints will be included as possible extrema otherwise they will not be included
        Output
        peakInds: Indices of peaks in x0
    */
    void findPeaks(vector <myflo>& x0, vector <int>& peakInds, bool includeEndpoints = true, int extrema = 1); // from PeakFinder.cpp
}
////////////////////////////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////////



////////////////////////////////////////////// GLOBAL VARIABLES //////////////////////////////////////////////

////////////////////////////////////////////// GLOBAL VARIABLES //////////////////////////////////////////////