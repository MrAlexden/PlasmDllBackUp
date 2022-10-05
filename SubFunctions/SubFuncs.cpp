#include "SubFuncs.h"

int NumPointsOfOriginalPila = 0,
	i_nach = 0,
	i_konec = 0;

myflo AverageSignal = 0,
	  koeff_MaxMin = 1.5;

myflo vSum(vector <myflo> & v)
{
	myflo sum = 0.0;
	for (int i = 0; i < v.size(); ++i)
		sum += v[i];
	return sum;
}

/////////////////// ������� ������� �� ���� ///////////////////
/*
����� ������ ������� ������ �������� ���������, ���� ������� ��������� ����� ���� (mnL_mxL_mnR_mxR[5]), ���� ������� �������� ����� ���� (mnL_mxL_mnR_mxR[4])
-- ��� ����� ����� ��������� ��������� �� ����� ��� ����� ��������� ��������� �������(����� ���)  --
&&
	����� ������ ������� && ����������� ����� ������ �������� ����������, �������� �������� ���� ����� (mnL_mxL_mnR_mxR[0]) ������������ �� ����������� (koeff_MaxMin)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
	||
	����� ������ ������� && ����������� ����� ������ �������� ����������, ������� ��������� ���� ����� (mnL_mxL_mnR_mxR[1]) ������������ �� ����������� (koeff_MaxMin)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
	||
	����� ������ ������� && ����������� ����� ������ �������� ����������, �������� �������� ���� ������ (mnL_mxL_mnR_mxR[2]) ������������ �� ����������� (koeff_MaxMin)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
	||
	����� ������ ������� && ����������� ����� ������ �������� ����������, ������� ��������� ���� ������ (mnL_mxL_mnR_mxR[3]) ������������ �� ����������� (koeff_MaxMin)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
*/
/////////////////// ������� ������� �� ���� ///////////////////

int find_i_nach(vector <myflo>& vSignal, vector <myflo>& mnL_mxL_mnR_mxR, vector <int>& vnIndices)
{
	myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05) ? 0.05 : otsech;
	int nach_ignored_gap = vnIndices[0] + (1 - otsech) * one_segment_width,
		konec_ignored_gap = vnIndices[0] + (1 + otsech) * one_segment_width;

	/*
	����� �� �� 1, � �� 2, ��� ��������� � 2 ���� �������(��������), � �������� �� ��������
	����� �����, ��� ������� ���� ������� ������� �� �������� (vnIndices[0] + one_segment_width*otsech)
	*/
	for (int i = vnIndices[0] + one_segment_width * otsech, k = 0; i < vnIndices[vnIndices.size() - 1]; i += ef_koef)
	{
		if ((i < nach_ignored_gap) || (i > konec_ignored_gap))
		{
			if (i > konec_ignored_gap)
			{
				k++;
				nach_ignored_gap = vnIndices[k] + (1 - otsech) * one_segment_width;
				konec_ignored_gap = vnIndices[k] + (1 + otsech) * one_segment_width;
			}

			if ((vSignal[i] < mnL_mxL_mnR_mxR[4] || vSignal[i] > mnL_mxL_mnR_mxR[5])
				&& ((vSignal[i] < mnL_mxL_mnR_mxR[0] && vSignal[i + 1] < mnL_mxL_mnR_mxR[0])
					|| (vSignal[i] > mnL_mxL_mnR_mxR[1] && vSignal[i + 1] > mnL_mxL_mnR_mxR[1])
					|| (vSignal[i] < mnL_mxL_mnR_mxR[2] && vSignal[i + 1] < mnL_mxL_mnR_mxR[2])
					|| (vSignal[i] > mnL_mxL_mnR_mxR[3] && vSignal[i + 1] > mnL_mxL_mnR_mxR[3])))
			{
				return i;
			}
		}
	}
	return 0;
}

int find_i_konec(vector <myflo>& vSignal, vector <myflo>& mnL_mxL_mnR_mxR, vector <int>& vnIndices)
{
	myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05) ? 0.05 : otsech;
	int nach_ignored_gap = vnIndices[vnIndices.size() - 1] - (1 - otsech) * one_segment_width,
		konec_ignored_gap = vnIndices[vnIndices.size() - 1] - (1 + otsech) * one_segment_width;

	/*
	����� �� �� 1, � �� 2, ��� ��������� � 2 ���� �������(��������), � �������� �� ��������
	����� �����, ��� ������� ���� ������� ������� �� �������� (vnIndices[vnIndices.GetSize()-1] - one_segment_width*otsech)
	*/
	for (int i = vnIndices[vnIndices.size() - 1] - one_segment_width * otsech, counter = 0, k = 0; i > vnIndices[0]; i -= ef_koef)
	{
		if ((i > nach_ignored_gap) || (i < konec_ignored_gap))
		{
			if (i < konec_ignored_gap)
			{
				k++;
				nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - (1 - otsech) * one_segment_width;
				konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - (1 + otsech) * one_segment_width;
			}

			if ((vSignal[i] < mnL_mxL_mnR_mxR[4] || vSignal[i] > mnL_mxL_mnR_mxR[5])
				&& ((vSignal[i] < mnL_mxL_mnR_mxR[0] && vSignal[i - 1] < mnL_mxL_mnR_mxR[0])
					|| (vSignal[i] > mnL_mxL_mnR_mxR[1] && vSignal[i - 1] > mnL_mxL_mnR_mxR[1])
					|| (vSignal[i] < mnL_mxL_mnR_mxR[2] && vSignal[i - 1] < mnL_mxL_mnR_mxR[2])
					|| (vSignal[i] > mnL_mxL_mnR_mxR[3] && vSignal[i - 1] > mnL_mxL_mnR_mxR[3])))
			{
				return i;
			}
		}
	}
	return 0;
}

void noise_vecs(int k, vector <myflo>& vLnoise, vector <myflo>& vRnoise, vector <myflo>& vSignal, vector <int>& vnIndices)
{
	myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05) ? 0.05 : otsech;
	int i = 0, j = 0, n = 0, nach_ignored_gap = 0, konec_ignored_gap = 0;

	i = vnIndices[k] - one_segment_width * (1 - otsech);
	if (i < 0) i = ef_koef;

	nach_ignored_gap = vnIndices[k] - otsech * one_segment_width;
	konec_ignored_gap = vnIndices[k] + otsech * one_segment_width;

	for (j = 0; j < vLnoise.size(); i += ef_koef)
	{
		if ((i < nach_ignored_gap) || (i > konec_ignored_gap))
		{
			if (i > konec_ignored_gap)
			{
				k++;
				nach_ignored_gap = vnIndices[k] - otsech * one_segment_width;
				konec_ignored_gap = vnIndices[k] + otsech * one_segment_width;
			}

			for (n = 0; n < ef_koef && j < vLnoise.size(); ++n, ++j)
				vLnoise[j] = vSignal[i - ef_koef + n];
		}
	}

	i = vnIndices[vnIndices.size() - 1 - k] + one_segment_width * (1 - otsech);
	if (i > vSignal.size() - 1) i = vSignal.size() - 1 - ef_koef;

	nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] + otsech * one_segment_width;
	konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - otsech * one_segment_width;

	for (j = 0; j < vRnoise.size(); i -= ef_koef)
	{
		if ((i > nach_ignored_gap) || (i < konec_ignored_gap))
		{
			if (i < konec_ignored_gap)
			{
				k++;
				nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] + otsech * one_segment_width;
				konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - otsech * one_segment_width;
			}

			for (n = 0; n < ef_koef && j < vRnoise.size(); ++n, ++j)
				vRnoise[j] = vSignal[i + ef_koef - n];
		}
	}

	myflo avL = vSum(vLnoise), avR = vSum(vRnoise);
	AverageSignal = (avL + avR) / (vLnoise.size() + vRnoise.size());
}

void MinMaxNoise(vector <myflo>& vLnoise, vector <myflo>& vRnoise, vector <myflo>& mnL_mxL_mnR_mxR)
{
	myflo mnL = *min_element(vLnoise.begin(), vLnoise.end()),
		mxL = *max_element(vLnoise.begin(), vLnoise.end()),
		mnR = *min_element(vRnoise.begin(), vRnoise.end()),
		mxR = *max_element(vRnoise.begin(), vRnoise.end());

	if (AverageSignal != 0)
	{
		mnL_mxL_mnR_mxR[0] = /*mnL*/AverageSignal - abs(AverageSignal - mnL) * koeff_MaxMin;
		mnL_mxL_mnR_mxR[1] = /*mxL*/AverageSignal + abs(mxL - AverageSignal) * koeff_MaxMin;
		mnL_mxL_mnR_mxR[2] = /*mnR*/AverageSignal - abs(AverageSignal - mnR) * koeff_MaxMin;
		mnL_mxL_mnR_mxR[3] = /*mxR*/AverageSignal + abs(mxR - AverageSignal) * koeff_MaxMin;
	}
	else
	{
		if (mnL < 0)
			mnL_mxL_mnR_mxR[0] = mnL * koeff_MaxMin;
		else
			mnL_mxL_mnR_mxR[0] = mnL / koeff_MaxMin;
		if (mxL > 0)
			mnL_mxL_mnR_mxR[1] = mxL * koeff_MaxMin;
		else
			mnL_mxL_mnR_mxR[1] = mxL / koeff_MaxMin;
		if (mnR < 0)
			mnL_mxL_mnR_mxR[2] = mnR * koeff_MaxMin;
		else
			mnL_mxL_mnR_mxR[2] = mnR / koeff_MaxMin;
		if (mxR > 0)
			mnL_mxL_mnR_mxR[3] = mxR * koeff_MaxMin;
		else
			mnL_mxL_mnR_mxR[3] = mxR / koeff_MaxMin;
	}

	////////////////////////////////////// ���������� ������� ������� ���� //////////////////////////////////////
	if (mxL > mxR) mnL_mxL_mnR_mxR[5] = AverageSignal + abs(mxL - AverageSignal);
	else mnL_mxL_mnR_mxR[5] = AverageSignal + abs(mxR - AverageSignal);

	if (mnL < mnR) mnL_mxL_mnR_mxR[4] = AverageSignal - abs(AverageSignal - mnL);
	else mnL_mxL_mnR_mxR[4] = AverageSignal - abs(AverageSignal - mnR);
	////////////////////////////////////// ���������� ������� ������� ���� //////////////////////////////////////
}

void fit_linear_pila(vector <myflo>& vPila, int ind_st_of_pil, vector <myflo>& vSegPila)
{
	int i = 0, j = 0, ind_st, ind_end;
	ind_st = (leftP >= 0.2) ? ind_st_of_pil + leftP * NumPointsOfOriginalPila : ind_st_of_pil + 0.2 * NumPointsOfOriginalPila;
	ind_end = (rightP >= 0.2) ? ind_st_of_pil + (1 - rightP) * NumPointsOfOriginalPila : ind_st_of_pil + (1 - 0.2) * NumPointsOfOriginalPila;
	myflo stVP = vPila[ind_st], endVP = vPila[ind_end], delta;
	vSegPila.resize((one_segment_width) * (1 - leftP - rightP));

	/* ���������� ���������� � ����� ����� ���� �� ���������� */
	myflo segment_length = endVP - stVP, segment_points_amount = abs(ind_st - ind_end), transform_factor = (myflo)one_segment_width / NumPointsOfOriginalPila;
	delta = segment_length / (segment_points_amount * transform_factor);

	vSegPila[0] = (leftP >= 0.2) ? stVP : vPila[ind_st_of_pil + leftP * NumPointsOfOriginalPila];
	for (i = 1; i < vSegPila.size(); ++i)
		vSegPila[i] = vSegPila[i - 1] + delta;
}

void convert_time_to_pts(int v_tok_size)
{
	myflo tot_time = (myflo)v_tok_size / one_segment_width / freqP;

	myflo st_ot_tot = st_time_end_time[0] / tot_time;
	myflo end_ot_tot = st_time_end_time[1] / tot_time;

	if (st_time_end_time[0] <= tot_time && st_time_end_time[1] <= tot_time)
	{
		i_nach = v_tok_size * st_ot_tot;
		i_konec = v_tok_size * end_ot_tot;
	}
}

int find_signal_and_make_pila(_In_ vector <myflo> & vPila,
							  _In_ vector <myflo> & vSignal, 
							  _Out_ vector <myflo> & vSegPila,
							  _Out_ vector <int> & vStartSegIndxs)
{
	vector <myflo> vLnoise, vRnoise;

	int i = 0,
		j = 0,
		k = 0,
		segments_amount = 0;

	/* ����� ��������� ��� ���� ����� ������ �������� */
	//vector <myflo> vSignal_copy(vSignal);
	//sg_smooth(vSignal_copy, vSignal, vSignal_copy.size() / 1000, 4);

	/* ������� ��� ������� �� ���� ��� ������ ����� */
	vector <int> vnIndices;
	PeakFinder::findPeaks(vPila, vnIndices, false, 1);
	if (vnIndices.size() <= 3) // 3 ������ ��� ������ � ������� � ������� ������ ������� � �� ����� �������, ������ �� ������� �� �������� �������
	{
		MessageBoxA(NULL, "Less then 4 segments found", "Error!", MB_ICONWARNING | MB_OK);
		return -1;
	}

	/* ���������� ����� ������� ����������� ���� */
	int todel = (int)abs(vnIndices[2] - vnIndices[1]) % 100;
	if ((myflo)todel / 100 <= 0.5) 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) - todel;
	else 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) + 100 - todel;

	/* ���������� ����� ������� ������� */
	if (vPila.size() != vSignal.size())
		one_segment_width = NumPointsOfOriginalPila * ((myflo)vSignal.size() / vPila.size());
	else
		one_segment_width = NumPointsOfOriginalPila;

	/* �������� ������ */
	if (one_segment_width == 0 
		|| NumPointsOfOriginalPila == 0 
		|| is_invalid(one_segment_width) 
		|| is_invalid(NumPointsOfOriginalPila)
		|| one_segment_width > vSignal.size()
		|| NumPointsOfOriginalPila > vSignal.size())
	{
		MessageBoxA(NULL, "Error in finding segments length", "Error!", MB_ICONWARNING | MB_OK);
		return -1;
	}

	/* ������� �������������� ���� � ������ �� ������ � ������ leftP � rightP */
	fit_linear_pila(vPila, vnIndices[1] + 5, vSegPila); // +5 ������ ��� �� ����� ���� ����, � � ���� ������� ���� �������� �� �������� �����, �������������� ����� ���������� ���� ����
	if (vSegPila.size() <= 0)
	{
		MessageBoxA(NULL, "Error in pila linearizing", "Error!", MB_ICONWARNING | MB_OK);
		return -1;
	}

	/* �������� ������� ����� �������� � ����� ���� � ����� */
	if (vPila.size() != vSignal.size())
		vectordiv(vnIndices, (myflo)vPila.size() / vSignal.size());
	/* ������� ���� ������� � ����� (�� ������ ������) � ��������� ���������� �������� ������� */
	vnIndices.pop_back();
	if (vnIndices[vnIndices.size() - 1] > vSignal.size()) 
		vectormult(vnIndices, (myflo)vPila.size() / vSignal.size());

	/* �����!!!!!!!!!!!
	��� ����� ������������, ��� ���� � ����� �������, ���� � ������ ���� ������� ����(����� ����� ��� �� ������� � ������ ������ ������������ ����)
	*/
	vLnoise.resize(one_segment_width);
	vRnoise.resize(one_segment_width);

	/* ��������� ������� ������� */
	noise_vecs(k, vLnoise, vRnoise, vSignal, vnIndices);

	vector <myflo> mnL_mxL_mnR_mxR(6);
	/*
	������ mnL_mxL_mnR_mxR:
		-[0] ������� ���� �����
		-[1] �������� ���� ������
		-[2] ������� ���� ������
		-[3] �������� ���� ������
		-[4] ������� ����� ����
		-[5] �������� ����� ����
	*/

	/* ��������� ������ mnL_mxL_mnR_mxR */
	MinMaxNoise(vLnoise, vRnoise, mnL_mxL_mnR_mxR);

	/* ���� ������������ ��� ������ �������� ��������� -> �� ������ � ���� while */
	bool istr = false;
	if (st_time_end_time[0] >= 0 && st_time_end_time[1] >= 0 && st_time_end_time[0] < st_time_end_time[1])
		convert_time_to_pts(vSignal.size());
	else istr = true;

	/* �������� ���� ������ ������� (������� �� ����) */
	while (segments_amount == 0 && istr)
	{
		k++;
		if (k > 5)
		{
			MessageBoxA(NULL, "More than 5 attepts to find signal", "Error!", MB_ICONWARNING | MB_OK);
			return -1;
		}

		i_nach = find_i_nach(vSignal, mnL_mxL_mnR_mxR, vnIndices);
		i_konec = find_i_konec(vSignal, mnL_mxL_mnR_mxR, vnIndices);

		if (i_nach < i_konec 
			&& i_nach != -1
			&& i_konec != -1
			|| ( i_nach != 0
			&& i_konec != 0 ))
			segments_amount = abs(i_konec - i_nach) / one_segment_width;
		else
		{
			MessageBoxA(NULL, "Error in finding start|end of signal", "Error!", MB_ICONWARNING | MB_OK);
			return -1;
		}

		/* ���� ���������� ���� �������� ������ 0.1 �� ������������ �������������, �� �� ������� ���� -> ������������� (����� ��������� ����������) */
		if (segments_amount <= (int)(0.1 * vSignal.size() / one_segment_width))
		{
			/* ��������� ������������ ����� �������� ������ ������ �������� */
			koeff_MaxMin -= 0.3;

			/* ����� ��������� ������� ������� */
			noise_vecs(k, vLnoise, vRnoise, vSignal, vnIndices);

			/* ����� ��������� ������ mnL_mxL_mnR_mxR */
			MinMaxNoise(vLnoise, vRnoise, mnL_mxL_mnR_mxR);

			/* �������� ����� �������� � ������/����� ����� ����� ����� � ���� while */
			segments_amount = 0;
			i_nach = 0;
			i_konec = 0;
		}
		/* ���� ���������� ���� �������� ������ 0.8 �� ������������ �������������, �� �� ������� ����� -> ������������� (����� ��������� ����������) */
		if (segments_amount >= (int)(0.8 * vSignal.size() / one_segment_width))
		{
			/* ����������� ������������ ����� �������� ������ ������ �������� */
			koeff_MaxMin += 0.3;

			/* ����� ��������� ������� ������� */
			noise_vecs(k, vLnoise, vRnoise, vSignal, vnIndices);

			/* ����� ��������� ������ mnL_mxL_mnR_mxR */
			MinMaxNoise(vLnoise, vRnoise, mnL_mxL_mnR_mxR);

			/* �������� ����� �������� � ������/����� ����� ����� ����� � ���� while */
			segments_amount = 0;
			i_nach = 0;
			i_konec = 0;
		}
	}

	/* �������� �� ������ ���� ������� ������� ����� ��������, �������� ������ - ���, � �� �������� ��� */
	if (abs(i_konec - i_nach) / one_segment_width >= (int)(0.8 * vSignal.size() / (one_segment_width)))
	{
		MessageBoxA(NULL, "Too many segments. Probably the signal is noise", "Error!", MB_ICONWARNING | MB_OK);
		return -1;
	}

	/* �����! ������ ������� ������� ����� �������� */
	vStartSegIndxs.resize(abs(i_konec - i_nach) / one_segment_width);

	/* ���������� ������� ����� �������� � ������ ������ � ����� ������� */
	for (i = 0, j = 0; j < vStartSegIndxs.size() && i < vnIndices.size(); ++i)
	{
		if (vnIndices[i] > i_nach && vnIndices[i] < i_konec)
		{
			vStartSegIndxs[j] = vnIndices[i];
			j++;
		}
	}

	/* ������� ������ � ����� */
	for (i = vStartSegIndxs.size() - 1; i > 1; --i)
	{
		if (vStartSegIndxs[i] != 0)
		{
			vStartSegIndxs.resize(i + 1);
			break;
		}
	}

	if (vStartSegIndxs.size() <= 0) 
	{
		MessageBoxA(NULL, "No segments found", "Error!", MB_ICONWARNING | MB_OK);
		return -1;
	}

	return 0;
}

myflo metod_Newton(myflo x0, vector <myflo>& vParams, myflo(*fx)(myflo, vector <myflo>&))
{
	myflo x1 = 0.0, x_cashe = 2 * x0;
	int n = 0;
	while (abs(x1 - x0) > TOL)
	{
		x1 = x0 - fx(x0, vParams) / partial_derivative(x0, vParams, fx, -1);
		x0 = x1 - fx(x1, vParams) / partial_derivative(x1, vParams, fx, -1);
		n++;
		if (x1 > x_cashe || x0 > x_cashe || n > 500)
		{
			x1 = x_cashe;
			break;
		}
	}
	return x1;
}

myflo dens(vector <myflo> & vPila, vector <myflo> & vParams)
{
	myflo x = metod_Newton(vPila[vPila.size() - 1], vParams, fx_STEP), ans = 0.0, M = 0.0;

	ans = abs(vParams[0] + vParams[1] * x);

	switch (fuel)
	{
	case 0:
		M = M_He;
		break;
	case 1:
		M = M_Ar;
		break;
	}

	return ans * (10e-6) / (0.52026 * S * (1.602 * 10e-19) * sqrt((1.380649 * 10e-23) * vParams[3] * 11604.51812 / M)); //��� �������
	//return ans * (10e-6) / (S * (1.602 * 10e-19) * sqrt(2 * vParams[3] * (1.602 * 10e-19) / M)); //������� ������
}

int make_one_segment(_In_ int diagnostics,					 // diagnostics type (zond::0|setka::1|cilind::2)
					 _In_  vector <myflo> & vPila,			 // X data
					 _In_  vector <myflo> & vSignal,		 // Y data
					 _Out_ vector <myflo> & vres,			 // vector to be filled with the result
					 _Out_ vector <myflo> & vfilt,			 // vector to be filled with the filtration
					 _Out_ vector <myflo> & vcoeffs)		 // additional coeffs/results vector
{
	myflo filtS = vcoeffs[0],
		linfitP = vcoeffs[1];

	switch (diagnostics)
	{
		case 0: // Zond
		{
			if (vPila.size() != vSignal.size())
			{
				MessageBoxA(NULL, "Input segment's values error", "Error!", MB_ICONWARNING | MB_OK);
				return -1;
			}

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

			if (linfitP != 0)
			{
				vector <myflo> vX, vY;

				vX.assign(vPila.begin(), vPila.begin() + vPila.size() * linfitP);
				vY.assign(vSignal.begin(), vSignal.begin() + vPila.size() * linfitP);

				vector <myflo> vAB = linear_fit(vX, vY);

				vParams[0] = vAB[0];
				vParams[1] = vAB[1];

				vFixed[0] = vFixed[1] = true;
			}
			else
			{
				vParams[0] = vPila[0];
				vParams[1] = 1;
			}

			vParams[2] = -vParams[0];
			vParams[3] = 10;

			if (linfitP != 0)
			{
				vector <myflo> vX, vY;

				vX.assign(vPila.begin() + vPila.size() * linfitP, vPila.end());
				vY.assign(vSignal.begin() + vPila.size() * linfitP, vSignal.end());

				levmarq(vX, vY, vParams, vFixed, Num_iter, fx_STEP);
			}
			else
				levmarq(vPila, vSignal, vParams, vFixed, Num_iter, fx_STEP);

			if (vParams[3] < 0 || vParams[3] > 100)
			{
				vParams[0] = -1;
				vParams[1] = 1;
				vParams[2] = 1;
				vParams[3] = 10;
				vFixed[0] = vFixed[1] = vFixed[2] = vFixed[3] = false;
				levmarq(vPila, vSignal, vParams, vFixed, Num_iter * 2, fx_STEP);
			}
			/*if (vParams[2] <= 0)
				vParams[2] = -vParams[0] / 4;*/

			vcoeffs = vParams;

			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_STEP(vPila[i], vParams);

			vcoeffs.push_back(dens(vPila, vParams));

			break;
		}
		case 1: // Setka
		{
			if (vPila.size() != vSignal.size())
			{
				MessageBoxA(NULL, "Input segment's values error", "Error!", MB_ICONWARNING | MB_OK);
				return -1;
			}

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());

			if (filtS == 0)
				filtS = 0.4;
			/* ��������� ������ � ���������� � vfilt */
			sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
			/* �������������� */
			diff(vfilt, vSignal, -1);
			/* ��������� ����������� */
			sg_smooth(vSignal, vres, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));

			vParams[0] = 0;													// y0 - offset
			vParams[1] = vPila[vPila.size() * 0.5];							// xc - center
			vParams[2] = abs(vPila[vPila.size() - 1] - vPila[0]) / 10;	// w - width
			vParams[3] = 1;													// A - area

			/* �������������� ����������� ������� */
			levmarq(vPila, vres, vParams, vFixed, Num_iter * 2, fx_GAUSS);

			vcoeffs.resize(4);
			vcoeffs[0] = vParams[0] + vParams[3] / (vParams[2] * sqrt(Pi / 2));	// Max Value
			vcoeffs[1] = vParams[2];												// Temp
			vcoeffs[2] = vParams[1];												// Peak Voltage
			vcoeffs[3] = sqrt(log(4)) * vParams[2];							// Energy

			/* ���������� � vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);

			break;
		}
		case 2: // Cilinder|Magnit
		{
			if (vPila.size() != vSignal.size())
			{
				MessageBoxA(NULL, "Input segment's values error", "Error!", MB_ICONWARNING | MB_OK);
				return -1;
			}

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

			vParams[0] = 0;													// y0 - offset
			vParams[1] = vPila[vPila.size() * 0.5];							// xc - center
			vParams[2] = abs(vPila[vPila.size() - 1] - vPila[0]) / 10;	// w - width
			vParams[3] = 1;													// A - area

			/* �������������� ����������� ������� */
			levmarq(vPila, vSignal, vParams, vFixed, Num_iter * 2, fx_GAUSS);

			vcoeffs.resize(4);
			vcoeffs[0] = vParams[0] + vParams[3] / (vParams[2] * sqrt(Pi / 2));	// Max Value
			vcoeffs[1] = vParams[2];												// Temp
			vcoeffs[2] = vParams[1];												// Peak Voltage
			vcoeffs[3] = sqrt(log(4)) * vParams[2];							// Energy

			/* ���������� � vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);

			break;
		}
	}

	return 0;
}

bool is_invalid(myflo val)
{
	if (val == NAN ||
		val == INFINITY ||
		val == -INFINITY ||
		val != val)
		return true;

	return false;
}
bool is_invalid(int val)
{
	if (val == NAN ||
		val == INFINITY ||
		val == -INFINITY ||
		val != val)
		return true;

	return false;
}

bool is_signalpeakslookingdown(vector <myflo>& v)
{
	myflo min = *min_element(v.begin(), v.end()), 
		  max = *max_element(v.begin(), v.end()),
		  vsum = vSum(v), 
		  vavg;
	vavg = vsum / v.size();

	if (abs(vavg - max) >= abs(vavg - min))
		return false;
	else
		return true;
}