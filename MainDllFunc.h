#pragma once

#define WIN32_LEAN_AND_MEAN // Исключите редко используемые компоненты из заголовков Windows
#define KLUDGE

#define TOL 10E-30f /* smallest value allowed in cholesky_decomp() */
#define Pi 3.14159265f
#define ef_koef 10 /* коэффициент эффективности, если == 1 то метод обрабатывает все точки, 
                если == 2 каждую вторую (в два раза быстре, но точность меньше) и тд \
                ef_koef используется в **SubFuncs** **Lev-Marq** */

#define ERR(f) if (err = (f), err < 0) goto Error;

/* Файлы заголовков Windows */
#include <windows.h>    // for windows interface
#include <iostream>     // for cin/cout
#include <fstream>      // for file access
#include <sstream>      // for file access
#include <vector>       // for vector
#include <tuple>        // for tuple
#include <algorithm>    // for min_element|max_element
//#include <omp.h>        // for multithreading
#include <thread>       // for multithreading
#include <list>         // for list
#include <functional>   // for function
#include <CommCtrl.h>   // for progressbar

#include "ProgressBarWindow/resource.h"   // progress bar resource
#include "Processing_Result_class/ProcessingResultClass.h"

using namespace std;    // for std::

typedef float myflo;

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

extern HWND mywindow;
extern HINSTANCE hInstThisDll;
////////////////////////////////////////////// GLOBAL VARIABLES //////////////////////////////////////////////

////////////////////////////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////////
/* GAUSS */
myflo fx_GAUSS(myflo, const vector <myflo> &); // from MainDllFunc.cpp

/* STEP */
myflo fx_STEP(myflo, const vector <myflo> &); // from MainDllFunc.cpp

/* LINE */
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
    for (int i = 0; i < in.size(); ++i)
        in[i] *= scalar;
};

/* divide the given vector on given scalar */
template <typename v, typename s>
inline void vectordiv(_Inout_ vector <v> & in, _In_ s scalar)
{
    if (scalar == 0) return;

    for (int i = 0; i < in.size(); ++i)
        in[i] /= scalar;
};

/* find derivtive */
void diff(_In_ const vector <myflo> &,
          _Out_ vector <myflo> &,
          _In_opt_ int);  // если -1 - то переворачивает функцию, если 1 - оставляет

/* find peaks */
namespace PeakFinder {
    const myflo EPS = 2.2204E-16f;

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

/* Levenberg-Marqurdt curve fitting method */
void levmarq(_In_ const vector <myflo> &,					 // independent data
             _In_ const vector <myflo> &,					 // dependent data
             _Inout_ vector <myflo> &,					     // vector of parameters we are searching for
             _In_ vector <bool> &,						     // vector of param's fixed status 
             _In_ unsigned int,							     // number of max iterations this method does
             _In_ myflo(*)(myflo, const vector <myflo> &));  // the function that the data will be approximated by

/* my function which finds usefull signal out of noise and creates voltage approximation */
int find_signal_and_make_pila(_In_ const vector <myflo> &,
                              _In_ const vector <myflo> &,
                              _Out_ vector <myflo> &,
                              _Out_ vector <int> &);    // from SubFuncs.cpp

inline int match_pila_and_signal(_In_ const vector <myflo> &,
                                 _In_ const vector <myflo> &,
                                 _Inout_ vector <int> &);     // from SubFuncs.cpp

/* savitzky golay smoothing */
void sg_smooth(_In_ const vector <myflo> &,
               _Out_ vector <myflo> &,
               _In_ const int,
               _In_ int);  // from Savitsky-Golay.cpp

/* processes one segment out of the whole signal */
extern "C" __declspec(dllexport) int make_one_segment(_In_ int,					            // diagnostics type (zond::0|setka::1|cilind::2)
                                                      _In_ const vector <myflo> &,	        // X data
                                                      _In_ vector <myflo> &,			    // Y data
                                                      _Out_ vector <myflo> &,		        // vector to be filled with the result
                                                      _Out_ vector <myflo> &,		        // vector to be filled with the filtration
                                                      _Out_ vector <myflo> &,		        // vector to be filled with the differentiation
                                                      _Inout_ vector <myflo> &);		    // additional coeffs/results vector

/* processes the whole Lengmur probe signal */
extern "C" __declspec(dllexport) int Zond(_In_ vector <myflo> & Pila,                       // входной одномерный массив пилы
                                          _In_ vector <myflo> & Signal,                     // входной одномерный массив сигнала
                                          _In_ const vector <myflo> & AdditionalData,       // дополнительные данные по импульсу
                                          _Out_ Plasma_proc_result <myflo> & fdata);        // выходной класс с результатом обработки

/* processes the whole potential analizator signal */
extern "C" __declspec(dllexport) int Setka(_In_ vector <myflo> & Pila,                      // входной одномерный массив пилы
                                           _In_ vector <myflo> & Signal,                    // входной одномерный массив сигнала
                                           _In_ const vector <myflo> & AdditionalData,      // дополнительные данные по импульсу
                                           _Out_ Plasma_proc_result <myflo> & fdata);       // выходной класс с результатом обработки

/* processes the whole cilinder analizator signal */
extern "C" __declspec(dllexport) int Cilinder(_In_ vector <myflo> & Pila,                   // входной одномерный массив пилы
                                              _In_ vector <myflo> & Signal,                 // входной одномерный массив сигнала
                                              _In_ const vector <myflo> & AdditionalData,   // дополнительные данные по импульсу
                                              _Out_ Plasma_proc_result <myflo> & fdata);    // выходной класс с результатом обработки

/* get this dll code description */
extern "C" __declspec(dllexport) string ERR_GetErrorDescription(int);    // from GlobalVarsInit.cpp

/* check whether a value valid or not */
template <typename T>
bool is_invalid(T val)
{
    if (val == NAN ||
        val == INFINITY ||
        val == -INFINITY ||
        val != val)
        return true;

    return false;
};

/* summitate a vector */
template <typename T>
inline T vSum(const vector <T>& v)
{
    T sum = 0;
    for (const auto it : v)
        sum += it;
    return sum;
};

/* check whether signal looks like ^ or V */
template <typename T>
bool is_signalpeakslookingdown(_In_ const vector <T> & v)
{
    T min = *min_element(v.begin(), v.end()),
      max = *max_element(v.begin(), v.end()),
      vsuml,
      vsumr,
      vavg;

    vector <T> vholder;

    vholder.assign(v.begin(), v.begin() + v.size() * 0.1);
    vsuml = vSum(vholder);

    vholder.assign(v.begin() + v.size() * 0.9, v.end());
    vsumr = vSum(vholder);

    vavg = (vsuml + vsumr) / (v.size() * 0.2);

    if (abs(vavg - max) >= abs(vavg - min))
        return false;
    else
        return true;
};

/* funcs for progress bar */
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
struct ThreadArgs
{
    HINSTANCE hInstance = hInstThisDll;
    LPCWSTR lpTemplateName = MAKEINTRESOURCE(IDD_DIALOG1);
    HWND hWndParent = NULL;
    DLGPROC lpDialogFunc = (DLGPROC)WndProc;
    LPARAM dwInitParam = NULL;
};
DWORD WINAPI DialogBoxParamWrapper(ThreadArgs* lpParameter);

/* get win32 errors description */
string GetLastErrorAsString();

// Error codes
typedef enum 
{
    ERR_NoError = 0,                // No error
    ERR_BadInputVecs = -6201,       // Corrupted input vectors
    ERR_ZeroInputVals = -6202,      // Input data error, some values equals 0
    ERR_BadCutOff = -6203,          // Сut-off points must be less than 0.9 in total and more than 0 each
    ERR_BadLinFit = -6204,          // The length of linear fit must be less than 0.9
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
////////////////////////////////////////////// GLOBAL FUNCTIONS //////////////////////////////////////////////