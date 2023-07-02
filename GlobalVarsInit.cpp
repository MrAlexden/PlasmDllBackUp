#include "MainDllFunc.h"

int freqP = 50,
	Num_iter = 10,
	one_segment_width = 0,
	fuel = 0; /* 0 -> He, 1 -> Ar, 2 -> Ne*/

myflo leftP = 0.05f,
	  rightP = 0.05f,
	  linfitP = 0.5f,
	  filtS = 0.0f,
	  st_time_end_time[2] = { -1.0f, -1.0f },
	  //S = Pi * 0.0005f * 0.005f + Pi * 0.0005f * 0.0005f / 4, // вся боковая поверхность цилиндра
      S = 0.0005f * 0.005f, // только проекция для набегающего потока
      M_He = 6.6464731E-27 - 9.10938356E-31,
	  M_Ar = 6.6335209E-26 - 9.10938356E-31,
      M_Ne = 3.3509177E-26 - 9.10938356E-31;

HWND mywindow = NULL;
HINSTANCE hInstThisDll = NULL;

/* GAUSS */
myflo fx_GAUSS(myflo x, const vector <myflo> & vParams)
{
	myflo y0 = vParams[0], xc = vParams[1], w = vParams[2], A = vParams[3];
	return y0 + (A / (w * sqrt(Pi / 2))) * exp(-0.5 * pow(((x - xc) / w), 2));
}
int GAUSS_InitParams(_In_ const vector <myflo> & vx, 
                     _In_ const vector <myflo> & vy, 
                     _Out_ vector <myflo> & vParams)
{
    /* берем ветора с краев аппроксимируемых данных чтобы потом найти по ним y0 */
    vector <myflo> vL, vR;
    take_data_around_the_edges(vy.size() / 10, vy, vL, vR);

    /* среднее между средним с краев */
    vParams[0] = (vSum(vL) / vL.size() + vSum(vR) / vR.size()) * 0.5;							// y0 - offset
    /* значение X максимума */
    vParams[1] = vx[max_element(vy.begin(), vy.end()) - vy.begin()];						    // xc - center
    /* вся ширина пилы делить на 10 */
    vParams[2] = abs(vx.back() - vx[0]) / 10;													// w - width
    /* выражено из формулы yc = y0 + A / (w * sqrt(Pi / 2)) */
    vParams[3] = (*max_element(vy.begin(), vy.end()) - vParams[0]) * (vParams[2] * 1.25331414);	// A - area

    return 0;
}

/* STEP */
myflo fx_STEP(myflo x, const vector <myflo> & vParams)
{
	myflo A = vParams[0], B = vParams[1], C = vParams[2], D = vParams[3];
	return A + B * x + C * exp(x / D);
}

/* LINE */
myflo fx_LINE(myflo x, const vector <myflo> & vParams)
{
	myflo A = vParams[0], B = vParams[1];
	return A + B * x;
}

/* EXPONENTIAL */
myflo fx_EXP(myflo x, const vector <myflo> & vParams)
{
    myflo y0 = vParams[0], A = vParams[1], R0 = vParams[2];
    return y0 + A * exp(x / R0);
}

/* PEACEWISE LINEAR THREE */
myflo fx_PWL3(myflo x, const vector <myflo> & vParams)
{
    myflo a1 = vParams[0],
          k1 = vParams[1],
          xi1 = vParams[2],
          k2 = vParams[3],
          xi2 = vParams[4],
          k3 = vParams[5];
    myflo yi1 = a1 + k1 * xi1;
    myflo yi2 = yi1 + k2 * (xi2 - xi1);

    if (x < xi1)
        return a1 + k1 * x;
    else if (x < xi2)
        return yi1 + k2 * (x - xi1);
    else
        return yi2 + k3 * (x - xi2);

    return 0;
}
int PWL3_InitParams(_In_ const vector <myflo> & vx,
                    _In_ const vector <myflo> & vy,
                    _Out_ vector <myflo> & vParams)
{
    if (vx.size() != vy.size())
        return -1;

    myflo a1 = vParams[0],
          k1 = vParams[1],
          xi1 = vParams[2],
          k2 = vParams[3],
          xi2 = vParams[4],
          k3 = vParams[5];

    int n1, n2, ni1, ni2;

    auto it_minx = min_element(vx.begin(), vx.end()),
         it_maxx = max_element(vx.begin(), vx.end());
    myflo x1 = *it_minx, 
          x2 = *it_maxx, 
          y1, 
          y2, 
          yi1, 
          yi2;

    xi1 = x1 + (x2 - x1) / 3;
    xi2 = x1 + 2 * (x2 - x1) / 3;

    y1 = vy[it_minx - vx.begin()];
    y2 = vy[it_maxx - vx.begin()];

    vector <myflo> vd(vx.begin(), vx.end());

    for (int i = 0; i < vx.size(); ++i)
        vd[i] = abs(vx[i] - xi1);

    yi1 = vy[min_element(vd.begin(), vd.end()) - vd.begin()];

    for (int i = 0; i < vx.size(); ++i)
        vd[i] = abs(vx[i] - xi2);

    yi2 = vy[min_element(vd.begin(), vd.end()) - vd.begin()];

    k1 = (yi1 - y1) / (xi1 - x1);
    a1 = y1 - k1 * x1;
    k2 = (yi2 - yi1) / (xi2 - xi1);
    k3 = (y2 - yi2) / (x2 - xi2);

    vParams[0] = a1;
    vParams[1] = k1;
    vParams[2] = xi1;
    vParams[3] = k2;
    vParams[4] = xi2;
    vParams[5] = k3;

    return 0;
}

string ERR_GetErrorDescription(int err)
{
    switch (err)
    {
    case ERR_BadInputVecs:
        return "Corrupted input vectors";
    case ERR_ZeroInputVals:
        return "Input data error, some values equals 0";
    case ERR_BadCutOff:
        return "Сut-off points must be less than 0.9 in total and more than 0 each";
    case ERR_BadLinFit:
        return "The length of linear fit must be less than 0.9";
    case ERR_BadFactorizing:
        return "Error after Ramp|Signal factorizing";
    case ERR_BadNoise:
        return "Error after noise extracting";
    case ERR_BadSegInput:
        return "Input segment's values error while segmend approximating";
    case ERR_TooFewSegs:
        return "Less then 4 segments found, check input arrays";
    case ERR_BadSegsLength:
        return "Error in finding segments length, check input params";
    case ERR_BadLinearRamp:
        return "Error in Ramp linearizing, check cut-off params";
    case ERR_TooManyAttempts:
        return "More than 5 attempts to find signal, check if signal is noise or not";
    case ERR_BadStartEnd:
        return "Error in finding start|end of signal, check if signal is noise or not";
    case ERR_TooManySegs:
        return "Too many segments, check if signal is noise or not";
    case ERR_NoSegs:
        return "No segments found, check if signal is noise or not";
    case ERR_BadDiagNum:
        return "Diagnostics number must be > 0 and < 3";
    case ERR_IdxOutOfRange:
        return "Index is out of range. Programm continued running";
    case ERR_Exception:
        return "An exception been occured while running, script did not stoped working";
    case ERR_BadStEndTime:
        return "Error!: start time must be less then end time and total time, more than 0\n\
end time must be less then total time, more then 0";
    case ERR_BufferExtension:
        return "Buffer extension, bad impulse";
    default:
        return "No Error";
    }
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    mywindow = hWnd;
    
    switch (message)
    {
    case WM_CLOSE:
        DestroyWindow(hWnd);
        break;
    case WM_DESTROY:
        PostQuitMessage(NULL);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }

    return 0;
}

DWORD WINAPI DialogBoxParamWrapper(ThreadArgs* lpParameter)
{
    if (DialogBoxParam(lpParameter->hInstance, // СЮДА НУЖНО КЛАСТЬ ТОЛЬКО HINSTANCE ЭТОЙ DLL, ИНАЧЕ НЕ РАБОТАЕТ
                       lpParameter->lpTemplateName,
                       lpParameter->hWndParent,
                       lpParameter->lpDialogFunc,
                       lpParameter->dwInitParam) < 0)
        MessageBoxA(NULL, GetLastErrorAsString().c_str(), "Error!", MB_ICONINFORMATION | MB_OK);

    return 0;
}

//Returns the last Win32 error, in string format. Returns an empty string if there is no error.
string GetLastErrorAsString()
{
    //Get the error message ID, if any.
    DWORD errorMessageID = GetLastError();
    if (errorMessageID == 0)
        return string(); //No error message has been recorded

    LPSTR messageBuffer = nullptr;

    //Ask Win32 to give us the string version of that message ID.
    //The parameters we pass in, tell Win32 to create the buffer that holds the message for us (because we don't yet know how long the message string will be).
    size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                                 NULL,
                                 errorMessageID,
                                 MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                                 (LPSTR)&messageBuffer, 
                                 NULL,
                                 NULL);

    //Copy the error message into a std::string.
    string message(messageBuffer, size);

    //Free the Win32's string's buffer.
    LocalFree(messageBuffer);

    return message;
}