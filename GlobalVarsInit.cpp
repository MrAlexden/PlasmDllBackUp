#include "MainDllFunc.h"

int freqP = 50,
	Num_iter = 10,
	one_segment_width = 0,
	fuel = 0;

myflo leftP = 0.05f,
	  rightP = 0.05f,
	  linfitP = 0.5f,
	  filtS = 0.0f,
	  st_time_end_time[2] = { -1.0f, -1.0f },
	  S = 3.141592f * 0.0005f * 0.005f + Pi * 0.0005f * 0.0005f / 4,
      M_He = 6.6464731f * 10E-27f - 9.10938356f * 10E-31f,
	  M_Ar = 6.6335209f * 10E-26f - 9.10938356f * 10E-31f;

/*===================================== GAUSS =====================================*/
myflo fx_GAUSS(myflo x, const vector <myflo> & vParams)
{
	myflo y0 = vParams[0], xc = vParams[1], w = vParams[2], A = vParams[3];
	return y0 + (A / (w * sqrt(Pi / 2))) * exp(-2 * pow(((x - xc) / w), 2));
}

/*===================================== STEP =====================================*/
myflo fx_STEP(myflo x, const vector <myflo> & vParams)
{
	myflo A = vParams[0], B = vParams[1], C = vParams[2], D = vParams[3];
	return A + B * x + C * exp(x / D);
}

/*===================================== LINE =====================================*/
myflo fx_LINE(myflo x, const vector <myflo> & vParams)
{
	myflo A = vParams[0], B = vParams[1];
	return A + B * x;
}

string ERR_GetErrorDescription(int err)
{
    switch (err)
    {
    case ERR_BadInputVecs:
        return "Corrupted input vectors";
        break;
    case ERR_ZeroInputVals:
        return "Input data error, some values equals 0";
        break;
    case ERR_BadCutOff:
        return "Ñut-off points must be less than 0.9 in total and more than 0 each";
        break;
    case ERR_BadLinFit:
        return "The length of linear fit must be less than 0.9";
        break;
    case ERR_BadFactorizing:
        return "Error after Pila|Signal factorizing";
        break;
    case ERR_BadNoise:
        return "Error after noise extracting";
        break;
    case ERR_BadSegInput:
        return "Input segment's values error while segmend approximating";
        break;
    case ERR_TooFewSegs:
        return "Less then 4 segments found, check input arrays";
        break;
    case ERR_BadSegsLength:
        return "Error in finding segments length, check input params";
        break;
    case ERR_BadLinearPila:
        return "Error in pila linearizing, check cut-off params";
        break;
    case ERR_TooManyAttempts:
        return "More than 5 attempts to find signal, check if signal is noise or not";
        break;
    case ERR_BadStartEnd:
        return "Error in finding start|end of signal, check if signal is noise or not";
        break;
    case ERR_TooManySegs:
        return "Too many segments, check if signal is noise or not";
        break;
    case ERR_NoSegs:
        return "No segments found, check if signal is noise or not";
        break;
    case ERR_BadDiagNum:
        return "Diagnostics number must be > 0 and < 2";
        break;
    case ERR_IdxOutOfRange:
        return "Index is out of range. Programm continued running";
        break;
    case ERR_Exception:
        return "An exception been occured while running, script did not stoped working";
        break;
    case ERR_BadStEndTime:
        return "Error!: start time must be less then end time and total time, more than 0\n\
end time must be less then total time, more then 0";
        break;
    default:
        return "No Error";
        break;
    }
}