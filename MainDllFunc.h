#pragma once

#define WIN32_LEAN_AND_MEAN // Исключите редко используемые компоненты из заголовков Windows

#define TOL 10E-30 /* smallest value allowed in cholesky_decomp() */
#define Pi 3.14159265
#define ef_koef 10 /* коэффициент эффективности, если == 1 то метод обрабатывает все точки, 
                если == 2 каждую вторую (в два раза быстре, но точность меньше) и тд \
                ef_koef используется в **SubFuncs** **Lev-Marq** */

#define ERR(f) if (err = (f), err < 0) goto Error;

//#ifdef __APPLE__
//# include <OpenCL/opencl.h>      // для компьютеров на MacOsX
//#else
//# include <CL/cl.h>              // для компьютеров на Win\Linux указывайте путь к файлу cl.h 
//#endif
//#ifndef __cl_clang_function_pointers
//#error "Missing __cl_clang_function_pointers define"
//#endif
//#pragma OPENCL EXTENSION __cl_clang_function_pointers : enable

typedef float myflo;

// Файлы заголовков Windows
#include <windows.h>    // for windows interface
#include <iostream>     // for cin/cout
#include <fstream>      // for file access
#include <sstream>      // for file access
#include <vector>       // for vector
#include <tuple>        // for tuple
#include <algorithm>    // for min_element|max_element
#include <omp.h>        // for multithreading
//#define _SILENCE_AMP_DEPRECATION_WARNINGS // for GPU multithreading
//#include <amp.h>                          // for GPU multithreading

#include "Processing_Result_class/ProcessingResultClass.h"

//using namespace concurrency; // for GPU multithreading
using namespace std;         // for std::

////////////////////////////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////////
/*===================================== GAUSS =====================================*/
myflo fx_GAUSS(myflo, const vector <myflo> &); // from MainDllFunc.cpp

/*===================================== STEP =====================================*/
myflo fx_STEP(myflo, const vector <myflo> &); // from MainDllFunc.cpp

/*===================================== LINE =====================================*/
myflo fx_LINE(myflo, const vector <myflo> &); // from MainDllFunc.cpp

/* make linear approximation of given data, returns vector(2) with A and B */
vector <myflo> linear_fit(_In_ const vector <myflo> &,
                          _In_ const vector <myflo> &);  // from Lev-Marq.cpp

/* multiply first and second vector, third elemnt is the result */
inline void vectorElementsProduct(_In_ const vector <myflo> &,
                                  _In_ const vector <myflo> &,
                                  _Out_ vector <myflo> &); // from PeakFinder.cpp

/* multiply the given scalar on given vector */
template <typename v, typename s>
inline void vectormult(_Inout_ vector <v> &in, _In_ s scalar)
{
#pragma omp parallel for schedule(static, 1) 
    for (int i = 0; i < in.size(); ++i)
        in[i] *= scalar;
};

/* divide the given vector on given scalar */
template <typename v, typename s>
inline void vectordiv(_Inout_ vector <v> & in, _In_ s scalar)
{
    if (scalar == 0) return;
#pragma omp parallel for schedule(static, 1) 
    for (int i = 0; i < in.size(); ++i)
        in[i] /= scalar;
};

void diff(_In_ const vector <myflo> &,
          _Out_ vector <myflo> &,
          _In_opt_ int);  // если -1 - то переворачивает функцию, если 1 - оставляет

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
    int findPeaks(_In_ const vector <myflo>&,
                  _Out_ vector <int>&,
                  _In_opt_ bool = true,
                  _In_opt_ int = 1);    // from PeakFinder.cpp
};

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
myflo partial_derivative(_In_ myflo,
                         _In_ const vector <myflo> &,
                         _In_ myflo(*)(myflo, const vector <myflo> &),
                         _In_ int);     // from Lev-Marq.cpp

void levmarq(_In_ const vector <myflo> &,					 // independent data
             _In_ const vector <myflo> &,					 // dependent data
             _Inout_ vector <myflo> &,					     // vector of parameters we are searching for
             _In_ vector <bool> &,						     // vector of param's fixed status 
             _In_ unsigned int,							     // number of max iterations this method does
             _In_ myflo(*)(myflo, const vector <myflo> &));  // the function that the data will be approximated by

int find_signal_and_make_pila(_In_ const vector <myflo> &,
                              _In_ const vector <myflo> &,
                              _Out_ vector <myflo> &,
                              _Out_ vector <int> &);    // from SubFuncs.cpp

/* savitzky golay smoothing */
void sg_smooth(_In_ const vector <myflo> &,
               _Out_ vector <myflo> &,
               _In_ const int,
               _In_ int);  // from Savitsky-Golay.cpp

int make_one_segment(_In_ int,					     // diagnostics type (zond::0|setka::1|cilind::2)
                     _In_ const vector <myflo> &,	 // X data
                     _In_ vector <myflo> &,			 // Y data
                     _Out_ vector <myflo> &,		 // vector to be filled with the result
                     _Out_ vector <myflo> &,		 // vector to be filled with the filtration
                     _Inout_ vector <myflo> &);		 // additional coeffs/results vector

extern "C" __declspec(dllexport) int Zond(_In_ vector <myflo> & Pila,                      // входной одномерный массив пилы
                                          _In_ vector <myflo> & Signal,                    // входной одномерный массив сигнала
                                          _In_ const vector <myflo> & AdditionalData,      // дополнительные данные по импульсу
                                          _Out_ Plasma_proc_result & fdata);               // выходной класс с результатом обработки
extern "C" __declspec(dllexport) int Setka(_In_ vector <myflo> & Pila,                     // входной одномерный массив пилы
                                           _In_ vector <myflo> & Signal,                   // входной одномерный массив сигнала
                                           _In_ const vector <myflo> & AdditionalData,     // дополнительные данные по импульсу
                                           _Out_ Plasma_proc_result & fdata);              // выходной класс с результатом обработки
extern "C" __declspec(dllexport) int Cilinder(_In_ vector <myflo> & Pila,                  // входной одномерный массив пилы
                                              _In_ vector <myflo> & Signal,                // входной одномерный массив сигнала
                                              _In_ const vector <myflo> & AdditionalData,  // дополнительные данные по импульсу
                                              _Out_ Plasma_proc_result & fdata);           // выходной класс с результатом обработки

bool is_invalid(myflo val); // from SubFuncs.cpp
bool is_invalid(int val);   // from SubFuncs.cpp

bool is_signalpeakslookingdown(_In_ const vector <myflo> & v); // from SubFuncs.cpp

// Error codes
typedef enum {

    ERR_NoError = 0,                // No error
    ERR_BadInputVecs = -6201,       // Corrupted input vectors
    ERR_ZeroInputVals = -6202,      // Input data error, some values equals 0
    ERR_BadCutOffLeft = -6203,      // Сut-off points on the left value must be > 0.0 and < 0.5
    ERR_BadCutOffRight = -6204,     // Сut-off points on the right value must be > 0.0 and < 0.5
    ERR_BadFactorizing = -6205,     // Error after Pila|Signal factorizing
    ERR_BadNoise = -6206,           // Error after noise extracting
    ERR_BadSegInput = -6207,        // Input segment's values error while segmend approximating
    ERR_TooFewSegs = -6208,         // Less then 4 segments found, check input arrays  
    ERR_BadSegsLength = -6209,      // Error in finding segments length, check input params
    ERR_BadLinearPila = -6210,      // Error in pila linearizing, check cut-off params
    ERR_TooManyAttempts = -6211,    // More than 5 attempts to find signal, check if signal is noise or not
    ERR_BadStartEnd = -6212,        // Error in finding start|end of signal, check if signal is noise or not
    ERR_TooManySegs = -6213,        // Too many segments, check if signal is noise or not
    ERR_NoSegs = -6214,             // No segments found, check if signal is noise or not
    ERR_BadDiagNum = -6215,         // Diagnostics number must be > 0 and < 2
    ERR_IdxOutOfRange = -6216,      // Index is out of range. Programm continued running
    ERR_Exception = -6217,          // An exception been occured whule running, script did not stopt working
    ERR_BadStEndTime = -6218        // Error!: start time must be less then end time and total time, more than 0\n\
                                    end time must be less then total time, more then 0
};

extern "C" __declspec(dllexport) string ERR_GetErrorDescription(int);
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
             S;
extern double M_Ar,
              M_He;
////////////////////////////////////////////// GLOBAL VARIABLES //////////////////////////////////////////////