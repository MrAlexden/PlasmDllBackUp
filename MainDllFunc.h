#pragma once

#define WIN32_LEAN_AND_MEAN // Исключите редко используемые компоненты из заголовков Windows

#define TOL 10E-30 /* smallest value allowed in cholesky_decomp() */
#define Pi 3.14159265
#define ef_koef 5 /* коэффициент эффективности, если == 1 то метод обрабатывает все точки, 
                если == 2 каждую вторую (в два раза быстре, но точность меньше) и тд */

// Файлы заголовков Windows
#include <windows.h>    // for windows interface
#include <iostream>     // for std::
#include <fstream>      // for file access
#include <sstream>      // for file access
#include <vector>       // for vector
#include <tuple>        // for tuple
#include <algorithm>    // for min_element|max_element
#include <omp.h>        // for multithreading

using namespace std;

typedef float myflo;

////////////////////////////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////////
/*===================================== GAUSS =====================================*/
myflo fx_GAUSS(myflo, vector <myflo>&); // from MainDllFunc.cpp

/*===================================== STEP =====================================*/
myflo fx_STEP(myflo, vector <myflo>&); // from MainDllFunc.cpp

/*===================================== LINE =====================================*/
myflo fx_LINE(myflo, vector <myflo>&); // from MainDllFunc.cpp

/* make linear approximation of given data, returns vector(2) with A and B */
vector <myflo> linear_fit(vector <myflo>&, vector <myflo>&); // from Lev-Marq.cpp

/* multiply first and second vector, third elemnt is the result */
void vectorElementsProduct(vector <myflo>&, vector <myflo>&, vector <myflo>&); // from PeakFinder.cpp

/* multiply the given scalar on given vector */
void vectormult(vector <myflo>&, int); // from PeakFinder.cpp
void vectormult(vector <myflo>&, myflo); // from PeakFinder.cpp
void vectormult(vector <int>&, int); // from PeakFinder.cpp
void vectormult(vector <int>&, myflo); // from PeakFinder.cpp

/* divide the given vector on given scalar */
void vectordiv(vector <myflo>&, int); // from PeakFinder.cpp
void vectordiv(vector <myflo>&, myflo); // from PeakFinder.cpp
void vectordiv(vector <int>&, int); // from PeakFinder.cpp
void vectordiv(vector <int>&, myflo); // from PeakFinder.cpp

void diff(_In_ vector <myflo>& in,
          _Out_ vector <myflo>& out,
          _In_ int k); // если -1 - то переворачивает функцию, если 1 - оставляет

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

class Matrix { // from matrix.cpp
private:
    unsigned m_rowSize;
    unsigned m_colSize;
    vector<vector<myflo> > m_matrix;
public:
    Matrix(unsigned, unsigned, myflo);
    Matrix(const char*);
    Matrix(const Matrix&);
    ~Matrix();

    // Matrix Operations
    Matrix operator+(Matrix&);
    Matrix operator-(Matrix&);
    Matrix operator*(Matrix&);
    Matrix transpose();

    // Scalar Operations
    Matrix operator+(myflo);
    Matrix operator-(myflo);
    Matrix operator*(myflo);
    Matrix operator/(myflo);

    // Aesthetic Methods
    myflo& operator()(const unsigned&, const unsigned&);
    void print() const;
    unsigned getRows() const;
    unsigned getCols() const;

    // Power Iteration
    tuple<Matrix, myflo, int> powerIter(unsigned, myflo);

    // Deflation
    Matrix deflation(Matrix&, myflo&);
};

/* calculate partial derivative of given func in particular point */
myflo partial_derivative(myflo, vector <myflo>&, myflo(*fx)(myflo, vector <myflo>&), int numofparam); // from Lev-Marq.cpp

void levmarq(vector <myflo>&,		        // independent data
    vector <myflo>&,				        // dependent data
    vector <myflo>&,			            // vector of parameters we are searching for
    vector <bool>&,				            // vector of param's fixed status 
    unsigned int,				            // number of max iterations this method does
    myflo (*fx)(myflo, vector <myflo>&));   // the function that the data will be approximated by

int find_signal_and_make_pila(vector <myflo> &, vector <myflo> &, vector <myflo> &, vector <int> &); // from SubFuncs.cpp

/* savitzky golay smoothing */
void sg_smooth(vector <myflo>&, vector <myflo>&, const int, const int); // from Savitsky-Golay.cpp

int make_one_segment(_In_ int,                      // diagnostics type (zond::0|setka::1|cilind::2)
                     _In_  vector <myflo>&,			// X data
                     _In_  vector <myflo>&,			// Y data
                     _Out_ vector <myflo>&,			// vector to be filled with the result
                     _Out_ vector <myflo>&,			// vector to be filled with the filtration
                     _Out_ vector <myflo>&);		// additional coeffs/results vector

extern "C" __declspec(dllexport) myflo *** Zond(_In_ vector <myflo> Pila,                // входной одномерный массив пилы
                                                _In_ vector <myflo> Signal,              // входной одномерный массив сигнала
                                                _In_ vector <myflo> AdditionalData);     // дополнительные данные по импульсу
extern "C" __declspec(dllexport) myflo *** Setka(_In_ vector <myflo> Pila,               // входной одномерный массив пилы
                                                 _In_ vector <myflo> Signal,             // входной одномерный массив сигнала
                                                 _In_ vector <myflo> AdditionalData);    // дополнительные данные по импульсу
extern "C" __declspec(dllexport) myflo *** Cilinder(_In_ vector <myflo> Pila,            // входной одномерный массив пилы
                                                    _In_ vector <myflo> Signal,          // входной одномерный массив сигнала
                                                    _In_ vector <myflo> AdditionalData); // дополнительные данные по импульсу

bool is_invalid(myflo val); // from SubFuncs.cpp
bool is_invalid(int val);   // from SubFuncs.cpp

bool is_signalpeakslookingdown(vector <myflo>& v); // from SubFuncs.cpp
////////////////////////////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////////



////////////////////////////////////////////// GLOBAL VARIABLES //////////////////////////////////////////////
extern int freqP,
           Num_iter,
           one_segment_width,
           fuel;

extern myflo leftP,
             rightP,
             linfitP,
             filtS,
             st_time_end_time[2],
             S,
             M_Ar,
             M_He;
////////////////////////////////////////////// GLOBAL VARIABLES //////////////////////////////////////////////