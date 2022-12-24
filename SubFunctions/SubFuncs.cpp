#include "SubFuncs.h"

int NumPointsOfOriginalPila = 0,
	i_nach = 0,
	i_konec = 0;

myflo AverageSignal = 0,
	  koeff_MaxMin = 1.5;

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

static inline int find_i_nach(_In_ const vector <myflo> & vSignal, 
							  _In_ const vector <myflo> & mnL_mxL_mnR_mxR,
							  _In_ const vector <int> & vnIndices)
{
	/*myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05f) ? 0.05f : otsech;*/
	int nach_ignored_gap = vnIndices[0] + (1 - rightP) * one_segment_width,
		konec_ignored_gap = vnIndices[0] + (1 + leftP) * one_segment_width;

	for (int i = vnIndices[0] + one_segment_width * leftP, k = 0; i < vnIndices.back(); i += ef_koef)
	{
		if ((i < nach_ignored_gap) || (i > konec_ignored_gap))
		{
			if (i > konec_ignored_gap)
			{
				++k;
				nach_ignored_gap = vnIndices[k] + (1 - rightP) * one_segment_width;
				konec_ignored_gap = vnIndices[k] + (1 + leftP) * one_segment_width;
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

static inline int find_i_konec(_In_ const vector <myflo> & vSignal,
							   _In_ const vector <myflo> & mnL_mxL_mnR_mxR,
							   _In_ const vector <int> & vnIndices)
{
	/*myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05f) ? 0.05f : otsech;*/
	int nach_ignored_gap = vnIndices.back() - (1 - rightP) * one_segment_width,
		konec_ignored_gap = vnIndices.back() - (1 + leftP) * one_segment_width;

	for (int i = vnIndices.back() - one_segment_width * rightP, counter = 0, k = 0; i > vnIndices[0]; i -= ef_koef)
	{
		if ((i > nach_ignored_gap) || (i < konec_ignored_gap))
		{
			if (i < konec_ignored_gap)
			{
				++k;
				nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - (1 - rightP) * one_segment_width;
				konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - (1 + leftP) * one_segment_width;
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

static inline void noise_vecs(_In_ int k,
							  _Out_ vector <myflo> & vLnoise, 
							  _Out_ vector <myflo> & vRnoise, 
							  _In_ const vector <myflo> & vSignal,
							  _In_ const vector <int> & vnIndices)
{
	/*myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05f) ? 0.05f : otsech;*/
	int i = 0, j = 0, n = 0, nach_ignored_gap = 0, konec_ignored_gap = 0;
	
	i = vnIndices[k] - one_segment_width * (1 - rightP);
	if (i < 0 || i < ef_koef) i = ef_koef;

	nach_ignored_gap = vnIndices[k] - rightP * one_segment_width;
	konec_ignored_gap = vnIndices[k] + leftP * one_segment_width;

	for (j = 0; j < vLnoise.size() && i < vSignal.size() - ef_koef; i += ef_koef)
	{
		if (i < nach_ignored_gap || i > konec_ignored_gap)
		{
			if (i > konec_ignored_gap)
			{
				++k;
				nach_ignored_gap = vnIndices[k] - rightP * one_segment_width;
				konec_ignored_gap = vnIndices[k] + leftP * one_segment_width;
			}

			for (n = 0; n < ef_koef && j < vLnoise.size(); ++n, ++j)
				vLnoise[j] = vSignal[i - ef_koef + n];
		}
	}

	i = vnIndices[vnIndices.size() - 1 - k] + one_segment_width * (1 - leftP);
	if (i > vSignal.size() - 1 || i > vSignal.size() - 1 - ef_koef) i = vSignal.size() - 1 - ef_koef;

	nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] + rightP * one_segment_width;
	konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - leftP * one_segment_width;

	for (j = 0; j < vRnoise.size() && i > ef_koef; i -= ef_koef)
	{
		if (i > nach_ignored_gap || i < konec_ignored_gap)
		{
			if (i < konec_ignored_gap)
			{
				++k;
				nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] + rightP * one_segment_width;
				konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - leftP * one_segment_width;
			}

			for (n = 0; n < ef_koef && j < vRnoise.size(); ++n, ++j)
				vRnoise[j] = vSignal[i + ef_koef - n];
		}
	}

	myflo avL = vSum(vLnoise), avR = vSum(vRnoise);
	AverageSignal = (avL + avR) / (vLnoise.size() + vRnoise.size());
}

static inline void MinMaxNoise(_In_ const vector <myflo> & vLnoise,
							   _In_ const vector <myflo> & vRnoise,
							   _Out_ vector <myflo> & mnL_mxL_mnR_mxR)
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

inline void fit_linear_pila(_In_ const vector <myflo> & vPila, 
							_In_ int ind_st_of_pil, 
							_Out_ vector <myflo> & vSegPila)
{
	int i = 0, j = 0, ind_st, ind_end;
	ind_st = (leftP >= 0.2) ? ind_st_of_pil + leftP * NumPointsOfOriginalPila : ind_st_of_pil + 0.2 * NumPointsOfOriginalPila;
	ind_end = (rightP >= 0.2) ? ind_st_of_pil + (1 - rightP) * NumPointsOfOriginalPila : ind_st_of_pil + (1 - 0.2) * NumPointsOfOriginalPila;
	myflo stVP = vPila[ind_st], endVP = vPila[ind_end], delta;
	vSegPila.resize((one_segment_width) * (1 - leftP - rightP));

	/* ���������� ���������� � ����� ����� ���� �� ���������� */
	myflo segment_length = endVP - stVP, segment_points_amount = abs(ind_st - ind_end), transform_factor = (myflo)one_segment_width / NumPointsOfOriginalPila;
	delta = segment_length / (segment_points_amount * transform_factor);

#define iterator vector <myflo>::iterator
	function<iterator(iterator, iterator)> minormax;

	if (stVP > endVP)
	{
		minormax = [](iterator begin, iterator end)
		{
			return max_element(begin, end);
		};
	}

	if (stVP < endVP)
	{
		minormax = [](iterator begin, iterator end)
		{
			return min_element(begin, end);
		};
	}
#undef iterator
	
	/* ����� ��������� ����� ������ ������ �������� �� ����������� �������� */
	vector <myflo> vPilaHolder(vPila.begin() + ind_st_of_pil, vPila.begin() + ind_st_of_pil + NumPointsOfOriginalPila);
	int index_of_minormax = minormax(vPilaHolder.begin(), vPilaHolder.end()) - vPilaHolder.begin();

	vSegPila[0] = *(vPila.begin() + ind_st_of_pil + index_of_minormax + leftP * NumPointsOfOriginalPila);

	for (i = 1; i < vSegPila.size(); ++i)
		vSegPila[i] = vSegPila[i - 1] + delta;

	if (vSegPila.back() < vSegPila.front())
		vectormult(vSegPila, -1);
}

static inline int convert_time_to_pts(_In_ int v_tok_size)
{
	if (st_time_end_time[0] == st_time_end_time[1])
		return 1;

	myflo tot_time = (myflo)v_tok_size / one_segment_width / freqP;

	if (st_time_end_time[0] < 0
		|| st_time_end_time[0] > tot_time
		|| st_time_end_time[0] > st_time_end_time[1]
		|| st_time_end_time[1] < 0
		|| st_time_end_time[1] > tot_time)
		return ERR_BadStEndTime;

	myflo st_ot_tot = st_time_end_time[0] / tot_time;
	myflo end_ot_tot = st_time_end_time[1] / tot_time;

	if (st_time_end_time[0] <= tot_time && st_time_end_time[1] <= tot_time)
	{
		i_nach = v_tok_size * st_ot_tot;
		i_konec = v_tok_size * end_ot_tot;
	}

	return 0;
}

inline int easy_find_peaks(_In_ const vector <myflo> & v,	// input data in which peaks are serched
						   _In_ const unsigned int step,	// the approximate length by which the peaks are distant 
						   _In_ const int min_max,			// 1 - max, -1 - min
						   _Out_ vector <int> & vnIndices)	// output vector with indexes
{
	if (min_max < -1 || min_max > 1 || min_max == 0)
		return -1;

#define iterator vector <myflo>::iterator
	function<iterator(iterator, iterator)> minormax;

	if (min_max == 1)
	{
		minormax = [](iterator begin, iterator end)
		{
			return max_element(begin, end);
		};
	}

	if (min_max == -1)
	{
		minormax = [](iterator begin, iterator end)
		{
			return min_element(begin, end);
		};
	}
#undef iterator

	vector <myflo> container;
	list <int> l;
	int k = 0;

	for (auto i = v.begin(), j = v.begin(); i - v.begin() < v.size() - 2 * step; i += step + (j - container.begin()))
	{
		container.assign(i, i + step);
		j = minormax(container.begin(), container.end());
		l.push_back((i - v.begin()) + (j - container.begin()));
	}

	vnIndices.assign(l.begin(), l.end());

	return 0;
}

inline int match_pila_and_signal(_In_ const vector <myflo> & vPila,
								 _In_ const vector <myflo> & vSignal,
								 _Inout_ vector <int> & vnPilaPeaks)
{
	/* ����� ��������� ����� ������� ������� int � ��������� ������, ��� ������� ���������� ���� */
	vector <myflo> IndHandler(vnPilaPeaks.begin(), vnPilaPeaks.end()); 
	if (vPila.size() != vSignal.size())
		(vPila.size() > vSignal.size()) ? vectordiv(IndHandler, (myflo)vPila.size() / vSignal.size())
										: vectormult(IndHandler, (myflo)vSignal.size() / vPila.size());

	/* ���������� ������ � �������� ������ ���� int */
	vnPilaPeaks.assign(IndHandler.begin(), IndHandler.end());

	/*��������� ���������� �������� �������*/
	while (vnPilaPeaks.back() > vSignal.size())
		vnPilaPeaks.pop_back();

	/* ������������� ����������, �� ��������� ��������� ����������� */
	//vector <int> vnSignalPeaks;
	//easy_find_peaks(vSignal, one_segment_width * 0.9, -1, vnSignalPeaks);

	return 0;
}

int find_signal_and_make_pila(_In_ const vector <myflo> & vPila,
							  _In_ const vector <myflo> & vSignal,
							  _Out_ vector <myflo> & vSegPila,
							  _Out_ vector <int> & vStartSegIndxs)
{
	vector <myflo> vLnoise, vRnoise;

	int i = 0,
		j = 0,
		k = 0,
		segments_amount = 0,
		err = 0;

	/* ����� ��������� ��� ���� ����� ������ �������� */
	//vector <myflo> vSignal_copy(vSignal);
	//sg_smooth(vSignal_copy, vSignal, vSignal_copy.size() / 1000, 4);

	/* ������� ��� ������� �� ���� ��� ������ ����� */
	vector <int> vnIndices;
	err = PeakFinder::findPeaks(vPila, vnIndices, false, 1);
	if (vnIndices.size() <= 3) // 3 ������ ��� ������ � ������� � ������� ������ ������� � �� ����� �������, ������ �� ������� �� �������� �������
		return ERR_TooFewSegs;

	/* ���������� ����� ������� ����������� ���� */
	int todel = (int)abs(vnIndices[2] - vnIndices[1]) % 100;
	if ((myflo)todel / 100 <= 0.5) 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) - todel;
	else 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) + 100 - todel;

	/* ���������� ����� ������� ������� */
	one_segment_width = (vPila.size() != vSignal.size()) ? NumPointsOfOriginalPila * ((myflo)vSignal.size() / vPila.size()) :
		NumPointsOfOriginalPila;

	/* �������� ������ */
	if (one_segment_width == 0 
		|| NumPointsOfOriginalPila == 0 
		|| is_invalid(one_segment_width) 
		|| is_invalid(NumPointsOfOriginalPila)
		|| one_segment_width > vSignal.size()
		|| NumPointsOfOriginalPila > vSignal.size())
		return ERR_BadSegsLength;

	/* ������� �������������� ���� � ������ �� ������ � ������ leftP � rightP */
	fit_linear_pila(vPila, vnIndices[1] + 5, vSegPila); // +5 ������ ��� �� ����� ���� ����, � � ���� ������� ���� �������� �� �������� �����, �������������� ����� ���������� ���� ����
	if (vSegPila.size() <= 0)
		return ERR_BadLinearPila;

	/* �������� ������� ����� �������� � ����� ���� � ����� */
	match_pila_and_signal(vPila, vSignal, vnIndices);

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
	int temp = convert_time_to_pts(vSignal.size());
	bool dowhile = true;
	if (temp < 0) err = temp;
	if (temp == 0) dowhile = false;

	/* �������� ���� ������ ������� (������� �� ����) */
	while (segments_amount == 0 && dowhile)
	{
		++k;
		if (k > 5)
			return ERR_TooManyAttempts;

		i_nach = find_i_nach(vSignal, mnL_mxL_mnR_mxR, vnIndices);
		i_konec = find_i_konec(vSignal, mnL_mxL_mnR_mxR, vnIndices);

		if (i_nach < i_konec 
			&& i_nach != -1
			&& i_konec != -1
			|| ( i_nach != 0
			&& i_konec != 0 ))
			segments_amount = abs(i_konec - i_nach) / one_segment_width;
		else
			return ERR_BadStartEnd;

		/* ���� ���������� ���� �������� ������ 0.1 �� ������������ �������������, �� �� ������� ���� -> ������������� (����� ��������� ����������) */
		if (segments_amount <= (int)(0.1f * vSignal.size() / one_segment_width))
		{
			/* ��������� ������������ ����� �������� ������ ������ �������� */
			koeff_MaxMin -= 0.3f;

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
		if (segments_amount >= (int)(0.8f * vSignal.size() / one_segment_width))
		{
			/* ����������� ������������ ����� �������� ������ ������ �������� */
			koeff_MaxMin += 0.3f;

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
	if (dowhile && ( abs(i_konec - i_nach) / one_segment_width >= (int)(0.9f * vSignal.size() / (one_segment_width)) ))
		return ERR_TooManySegs;

	/* ���������� ������� ����� �������� � ������ ������ � ����� ������� */
	for (i = 0; i < vnIndices.size(); ++i)
	{
		if (vnIndices[i] > i_nach)
		{
			vStartSegIndxs.assign(vnIndices.begin() + i, vnIndices.begin() + i + abs(i_konec - i_nach) / one_segment_width);
			break;
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
		return ERR_NoSegs;

	return err;
}

inline myflo metod_hord(_In_ myflo x0,
						_In_ myflo x1,
						_In_ const vector <myflo> & vParams,
						_In_ myflo(*fx)(myflo, const vector <myflo> &))
{
	int n = 0;
	while (abs(x1 - x0) > 10E-10 && n < 100)
	{
		x0 = x1 - (x1 - x0) * fx(x1, vParams) / (fx(x1, vParams) - fx(x0, vParams));
		x1 = x0 - (x0 - x1) * fx(x0, vParams) / (fx(x0, vParams) - fx(x1, vParams));
		++n;
	}
	return x1;
}

inline myflo metod_Newton(_In_ myflo x0,
						  _In_ const vector <myflo> & vParams,
						  _In_ myflo(*fx)(myflo, const vector <myflo> &))
{
	myflo x1 = 0.0;
	int n = 0;
	while (abs(x1 - x0) > 10E-10 && n < 100)
	{
		x1 = x0 - fx(x0, vParams) / partial_derivative(x0, vParams, fx, -1);
		x0 = x1 - fx(x1, vParams) / partial_derivative(x1, vParams, fx, -1);
		++n;
	}
	return x1;
}

inline myflo find_floating_potential(_In_ const vector <myflo> & vPila,
									 _In_ const vector <myflo> & vParams)
{
	myflo x = metod_Newton(vPila.back(), vParams, fx_STEP);

	if (is_invalid(x) || !x)
	{
		x = metod_hord(vPila[0], vPila.back(), vParams, fx_STEP);
		if (is_invalid(x) || !x)
			x = vPila.back();
	}

	return x;
}

inline myflo find_density(_In_ const myflo ion_current,
						  _In_ const myflo temperature)
{
	myflo M = 0.0f;
	
	switch (fuel)
	{
	case 0:
		M = M_He;
		break;
	case 1:
		M = M_Ar;
		break;
	default:
		__assume(0);
	}

	return abs(ion_current) * 1E-6f / (0.52026f * S * e * sqrt(kb * temperature * eVtoK / M)); //��� �������
	//return abs(ion_current) * 1E-6f / (S * e * sqrt(2 * temperature * e / M)); //������� ������
}

inline myflo find_ion_current(_In_ const myflo A,
							  _In_ const myflo B,
							  _In_ const myflo x)
{
	myflo ic = abs(A + B * x);

	if (is_invalid(ic) || !ic)
		ic = A;

	return ic;
}

inline myflo find_width_at_half_maximum(_In_ const vector <myflo> vx, _In_ const vector <myflo> vy)
{
	auto max = max_element(vy.begin(), vy.end()),
		min = min_element(vy.begin(), vy.end()),
		itl = max,
		itr = max;

	for (bool lflag = NULL, rflag = NULL; itl != vy.begin() && itr != vy.end();)
	{
		if (*itl > *min + (*max - *min) / 2)
			itl--;
		else 
			lflag = true;

		if (*itr > *min + (*max - *min) / 2)
			itr++;
		else
			rflag = true;

		if (lflag && rflag)
			break;
	}

	return vx[itr - vy.begin()] - vx[itl - vy.begin()];
}

int make_one_segment(_In_ const int diagnostics,			 // diagnostics type (zond::0|setka::1|cilind::2)
					 _In_ const vector <myflo> & vPila,		 // X data
					 _In_ vector <myflo> & vSignal,			 // Y data
					 _Out_ vector <myflo> & vres,			 // vector to be filled with the result
					 _Out_ vector <myflo> & vfilt,			 // vector to be filled with the filtration
					 _Out_ vector <myflo> & vdiff,			 // vector to be filled with the differentiation
					 _Inout_ vector <myflo> & vcoeffs)		 // additional coeffs/results vector
{
	S = vcoeffs[0];
	myflo linfitP = vcoeffs[1],
		  filtS = vcoeffs[2];
	fuel = (int)vcoeffs[3];
	Num_iter = (int)vcoeffs[4];

	switch (diagnostics)
	{
		case 0: // Zond
		{
			if (vPila.size() != vSignal.size() 
				|| vPila.empty() 
				|| vSignal.empty())
				return ERR_BadSegInput;

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());
			vdiff.resize(vPila.size());

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
				vParams[0] = vSignal[0];
				vParams[1] = 1;
			}

			vParams[2] = -vParams[0];
			vParams[3] = 15;

			if (linfitP != 0)
			{
				vector <myflo> vX, vY;

				vX.assign(vPila.begin() + vPila.size() * linfitP, vPila.end());
				vY.assign(vSignal.begin() + vPila.size() * linfitP, vSignal.end());

				levmarq(vX, vY, vParams, vFixed, Num_iter, fx_STEP);
			}
			else
				levmarq(vPila, vSignal, vParams, vFixed, Num_iter, fx_STEP);

			if (vParams[3] < 0.0 
				|| vParams[3] > 100.0
				|| is_invalid(vParams[3]))
			{
				vector <myflo> vX, vY;

				vX.assign(vPila.begin(), vPila.begin() + vPila.size() * 0.5);
				vY.assign(vSignal.begin(), vSignal.begin() + vPila.size() * 0.5);

				vector <myflo> vAB = linear_fit(vX, vY);

				vParams[0] = vAB[0];
				vParams[1] = vAB[1];
				vParams[2] = -vParams[0];
				vParams[3] = 15;

				vFixed[0] = vFixed[1] = true;

				vX.assign(vPila.begin() + vPila.size() * (1 - 0.5), vPila.end());
				vY.assign(vSignal.begin() + vPila.size() * (1 - 0.5), vSignal.end());

				levmarq(vX, vY, vParams, vFixed, max(Num_iter, 400), fx_STEP);

				if (vParams[3] < 0.0
					|| vParams[3] > 100.0
					|| is_invalid(vParams[3]))
				{
					vParams[2] = - vAB[0] / 1E5;
					vParams[3] = 15;
				}
			}

			/* ���������� � vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_STEP(vPila[i], vParams);

			/* ���������� � vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = find_floating_potential(vPila, vParams);							// Floating potential
			vcoeffs[1] = find_ion_current(vParams[0], vParams[1], vcoeffs[0]);				// Saturate ion current
			vcoeffs[2] = vParams[3];														// Temperature
			vcoeffs[3] = find_density(vcoeffs[1], vcoeffs[2]);								// Density ne

			break;
		}
		case 1: // Setka
		{
			if (vPila.size() != vSignal.size()
				|| vPila.empty()
				|| vSignal.empty())
				return ERR_BadSegInput;

			int i = 0;

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());
			vdiff.resize(vPila.size());

			if (filtS == 0)
				filtS = 0.4f;
			/* ��������� ������ � ���������� � vfilt */
			sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
			/* �������������� (������ �����������) */
			diff(vfilt, vdiff, -1); // -1 -> ����� ����� �����������
			/* ��������� ������ ����������� */
			sg_smooth(vdiff, vres, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
			vdiff = vres;

#ifdef GAUSSFIT
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			/* �������������� ���������� ��������� �-�� ������ */
			GAUSS_InitParams(vPila, vres, vParams);

			/* �������������� ����������� ������� */
			levmarq(vPila, vres, vParams, vFixed, Num_iter/* * 2*/, fx_GAUSS);

			if (vParams[2] < 0.0 
				|| vParams[2] > abs(vPila.back() - vPila[0]) / 2
				|| is_invalid(vParams[2]))
			{
				filtS = max(filtS, 0.8f);

				/* ��������� ������ � ���������� � vfilt */
				sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				/* �������������� */
				diff(vfilt, vdiff, -1);
				/* ��������� ����������� */
				sg_smooth(vdiff, vres, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vdiff = vres;

				/* �������������� ���������� ��������� �-�� ������ */
				GAUSS_InitParams(vPila, vres, vParams);

				levmarq(vPila, vres, vParams, vFixed, min(Num_iter, 2), fx_GAUSS);

				if (vParams[2] < 0.0
					|| vParams[2] > 100.0
					|| is_invalid(vParams[2]))
					vParams[2] = 15;
			}

			/* ���������� � vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);
#else
			vector <bool> vFixed = { false, false, false };
			vector <myflo> vParams(3);
			vector <myflo> vhandler(vPila.size());

			int middlepoint_of_curvature = max_element(vdiff.begin(), vdiff.end()) - vdiff.begin();

			vhandler.assign(vdiff.begin() + middlepoint_of_curvature, vdiff.end());

			/* �������������� (������ �����������) */
			diff(vhandler, vres, -1); // -1 -> ����� ����� �����������
			/* ��������� ������ ����������� */
			sg_smooth(vres, vhandler, vPila.size()* (filtS / 10), (5/*OriginLab poly order*/ - 1));
			
			int endpoint_of_curvature = middlepoint_of_curvature + 4 * abs(max_element(vhandler.begin(), vhandler.end()) - vhandler.begin());

			if (endpoint_of_curvature < vPila.size() && middlepoint_of_curvature < endpoint_of_curvature)
			{
				vhandler.assign(vPila.begin() + middlepoint_of_curvature, vPila.begin() + endpoint_of_curvature);
				vres.assign(vfilt.begin() + middlepoint_of_curvature, vfilt.begin() + endpoint_of_curvature);
			}
			else
			{
				vhandler.assign(vPila.begin() + middlepoint_of_curvature, vPila.end());
				vres.assign(vfilt.begin() + middlepoint_of_curvature, vfilt.end());
			}

			vParams[0] = *min_element(vres.begin(), vres.end());
			vParams[1] = abs(vParams[0]);
			vParams[2] = -15;

			/* �������������� ��������� ������� ����������� */
			levmarq(vhandler, vres, vParams, vFixed, max(Num_iter, 100), fx_EXP);

			if (vParams[2] > 0.0
				|| vParams[2] < -100.0
				|| is_invalid(vParams[2]))
			{
				vParams[2] = -15;
			}

			vres.resize(vPila.size());
			/* ���������� � vres */
			for (i = 0; i < vPila.size(); ++i)
			{
				if (i >= middlepoint_of_curvature &&
					i <= middlepoint_of_curvature + vhandler.size())
					vres[i] = fx_EXP(vPila[i], vParams);
				else
					vres[i] = NULL;
			}
#endif

			/* ���������� � vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = *max_element(vSignal.begin(), vSignal.end());						// Max Value
			vcoeffs[1] = abs(vParams[2]);													// Temp
			vcoeffs[2] = vPila[max_element(vdiff.begin(), vdiff.end()) - vdiff.begin()];	// Peak Voltage
			vcoeffs[3] = find_width_at_half_maximum(vPila, vdiff);							// Energy

			break;
		}
		case 2: // Cilinder|Magnit
		{
			if (vPila.size() != vSignal.size()
				|| vPila.empty()
				|| vSignal.empty())
				return ERR_BadSegInput;

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());
			vdiff.resize(vPila.size());

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

			/* �������������� ���������� ��������� �-�� ������ */
			GAUSS_InitParams(vPila, vSignal, vParams);

			/* �������������� ����������� ������� */
			levmarq(vPila, vSignal, vParams, vFixed, Num_iter/* * 2*/, fx_GAUSS);

			if (vParams[2] < 0.0 
				|| vParams[2] > abs(vPila.back() - vPila[0]) / 2
				|| is_invalid(vParams[2]))
			{
				filtS = max(filtS, 0.8f);

				/* ��������� ������ � ���������� � vfilt */
				sg_smooth(vSignal, vres, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));

				/* �������������� ���������� ��������� �-�� ������ */
				GAUSS_InitParams(vPila, vres, vParams);

				levmarq(vPila, vres, vParams, vFixed, min(Num_iter, 2), fx_GAUSS);

				if (vParams[2] < 0.0
					|| vParams[2] > 100.0
					|| is_invalid(vParams[2]))
					vParams[2] = 15;
			}

			/* ���������� � vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);

			/* ���������� � vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = *max_element(vSignal.begin(), vSignal.end());		// Max Value
			vcoeffs[1] = vParams[2];										// Temp
			vcoeffs[2] = abs(vParams[1]);									// Peak Voltage
			vcoeffs[3] = find_width_at_half_maximum(vPila, vSignal);		// Energy

			break;
		}
		case 3:	// double langmuir probe
		{
			if (vPila.size() != vSignal.size()
				|| vPila.empty()
				|| vSignal.empty())
				return ERR_BadSegInput;

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());
			vdiff.resize(vPila.size());

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

			PWL3_InitParams(vPila, vSignal, vParams);

			levmarq(vPila, vSignal, vParams, vFixed, Num_iter, fx_PWL3);

			if (vParams[3] < 0.0
				|| vParams[3] > 100.0
				|| is_invalid(vParams[3]))
			{
				vector <myflo> vX, vY;

				vX.assign(vPila.begin(), vPila.begin() + vPila.size() * 0.5);
				vY.assign(vSignal.begin(), vSignal.begin() + vPila.size() * 0.5);

				vector <myflo> vAB = linear_fit(vX, vY);

				vParams[0] = vAB[0];
				vParams[1] = vAB[1];
				vParams[2] = -vParams[0];
				vParams[3] = 15;

				vFixed[0] = vFixed[1] = true;

				vX.assign(vPila.begin() + vPila.size() * (1 - 0.5), vPila.end());
				vY.assign(vSignal.begin() + vPila.size() * (1 - 0.5), vSignal.end());

				levmarq(vX, vY, vParams, vFixed, max(Num_iter, 400), fx_STEP);

				if (vParams[3] < 0.0
					|| vParams[3] > 100.0
					|| is_invalid(vParams[3]))
				{
					vParams[2] = -vAB[0] / 1E5;
					vParams[3] = 15;
				}
			}

			/* ���������� � vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_STEP(vPila[i], vParams);

			/* ���������� � vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = find_floating_potential(vPila, vParams);							// Floating potential
			vcoeffs[1] = find_ion_current(vParams[0], vParams[1], vcoeffs[0]);				// Saturate ion current
			vcoeffs[2] = vParams[3];														// Temperature
			vcoeffs[3] = find_density(vcoeffs[1], vcoeffs[2]);								// Density ne

			break;
		}
		default:
			__assume(0);
	}

	return 0;
}