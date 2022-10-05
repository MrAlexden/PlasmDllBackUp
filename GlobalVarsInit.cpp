#include "MainDllFunc.h"

int freqP = 50,
	Num_iter = 10,
	one_segment_width = 0,
	fuel = 0;

myflo leftP = 0.05,
	  rightP = 0.05,
	  linfitP = 0.5,
	  filtS = 0.0,
	  st_time_end_time[2] = { -1, -1 },
	  S = 3.141592 * 0.0005 * 0.005 + 3.141592 * 0.0005 * 0.0005 / 4,
	  M_He = 6.6464731 * 10E-27 - 9.10938356 * 10E-31,
	  M_Ar = 6.6335209 * 10E-26 - 9.10938356 * 10E-31;

/*===================================== GAUSS =====================================*/
myflo fx_GAUSS(myflo x, vector <myflo>& vParams)
{
	myflo y0 = vParams[0], xc = vParams[1], w = vParams[2], A = vParams[3];
	return y0 + (A / (w * sqrt(Pi / 2))) * exp(-2 * pow(((x - xc) / w), 2));
}

/*===================================== STEP =====================================*/
myflo fx_STEP(myflo x, vector <myflo>& vParams)
{
	myflo A = vParams[0], B = vParams[1], C = vParams[2], D = vParams[3];
	return A + B * x + C * exp(x / D);
}

/*===================================== LINE =====================================*/
myflo fx_LINE(myflo x, vector <myflo>& vParams)
{
	myflo A = vParams[0], B = vParams[1];
	return A + B * x;
}