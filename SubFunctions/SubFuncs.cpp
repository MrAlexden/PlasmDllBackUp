#include "SubFuncs.h"

int NumPointsOfOriginalRamp = 0,
	i_nach = 0,
	i_konec = 0;

myflo AverageSignal = 0,
	  koeff_MaxMin = 1.5;

/////////////////// УСЛОВИЕ ОЧИСТКИ ОТ ШУМА ///////////////////
/*
Точка начала сигнала должна обладать значением, либо большим максимума всего шума (mnL_mxL_mnR_mxR[5]), либо меньшим минимума всего шума (mnL_mxL_mnR_mxR[4])
-- это нужно чтобы исключить вхождение по массе при очень маленькой амплитуде сигнала(почти шум)  --
&&
	Точка начала сигнала && последующая точка должны обладать значениями, меньшими минимума шума слева (mnL_mxL_mnR_mxR[0]) домноженного на коэффициент (koeff_MaxMin)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
	||
	Точка начала сигнала && последующая точка должны обладать значениями, большим максимума шума слева (mnL_mxL_mnR_mxR[1]) домноженного на коэффициент (koeff_MaxMin)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
	||
	Точка начала сигнала && последующая точка должны обладать значениями, меньшими минимума шума справа (mnL_mxL_mnR_mxR[2]) домноженного на коэффициент (koeff_MaxMin)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
	||
	Точка начала сигнала && последующая точка должны обладать значениями, большим максимума шума справа (mnL_mxL_mnR_mxR[3]) домноженного на коэффициент (koeff_MaxMin)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
*/
/////////////////// УСЛОВИЕ ОЧИСТКИ ОТ ШУМА ///////////////////

static inline int find_i_nach(_In_ const vector <myflo> & vSignal, 
							  _In_ const vector <myflo> & mnL_mxL_mnR_mxR,
							  _In_ const vector <int> & vnIndices)
{
	/*myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05f) ? 0.05f : otsech;*/
	int nach_ignored_gap = vnIndices[0] + (1 - rightP) * one_segment_width,
		konec_ignored_gap = vnIndices[0] + (1 + leftP) * one_segment_width;

	for (int i = vnIndices[0] + one_segment_width * leftP, k = 0; i < vnIndices.back(); i += boost)
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

	for (int i = vnIndices.back() - one_segment_width * rightP, counter = 0, k = 0; i > vnIndices[0]; i -= boost)
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
	if (i < 0 || i < boost) i = boost;

	nach_ignored_gap = vnIndices[k] - rightP * one_segment_width;
	konec_ignored_gap = vnIndices[k] + leftP * one_segment_width;

	for (j = 0; j < vLnoise.size() && i < vSignal.size() - boost; i += boost)
	{
		if (i < nach_ignored_gap || i > konec_ignored_gap)
		{
			if (i > konec_ignored_gap)
			{
				++k;
				nach_ignored_gap = vnIndices[k] - rightP * one_segment_width;
				konec_ignored_gap = vnIndices[k] + leftP * one_segment_width;
			}

			for (n = 0; n < boost && j < vLnoise.size(); ++n, ++j)
				vLnoise[j] = vSignal[i - boost + n];
		}
	}

	i = vnIndices[vnIndices.size() - 1 - k] + one_segment_width * (1 - leftP);
	if (i > vSignal.size() - 1 || i > vSignal.size() - 1 - boost) i = vSignal.size() - 1 - boost;

	nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] + rightP * one_segment_width;
	konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - leftP * one_segment_width;

	for (j = 0; j < vRnoise.size() && i > boost; i -= boost)
	{
		if (i > nach_ignored_gap || i < konec_ignored_gap)
		{
			if (i < konec_ignored_gap)
			{
				++k;
				nach_ignored_gap = vnIndices[vnIndices.size() - 1 - k] + rightP * one_segment_width;
				konec_ignored_gap = vnIndices[vnIndices.size() - 1 - k] - leftP * one_segment_width;
			}

			for (n = 0; n < boost && j < vRnoise.size(); ++n, ++j)
				vRnoise[j] = vSignal[i + boost - n];
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

	////////////////////////////////////// Определяем средний уровень шума //////////////////////////////////////
	if (mxL > mxR) mnL_mxL_mnR_mxR[5] = AverageSignal + abs(mxL - AverageSignal);
	else mnL_mxL_mnR_mxR[5] = AverageSignal + abs(mxR - AverageSignal);

	if (mnL < mnR) mnL_mxL_mnR_mxR[4] = AverageSignal - abs(AverageSignal - mnL);
	else mnL_mxL_mnR_mxR[4] = AverageSignal - abs(AverageSignal - mnR);
	////////////////////////////////////// Определяем средний уровень шума //////////////////////////////////////
}

inline void fit_linear_Ramp(_In_ const vector <myflo> & vRamp, 
							_In_ int ind_st_of_pil, 
							_Out_ vector <myflo> & vSegRamp)
{
	int i = 0, j = 0, ind_st, ind_end;
	ind_st = (leftP >= 0.2) ? ind_st_of_pil + leftP * NumPointsOfOriginalRamp : ind_st_of_pil + 0.2 * NumPointsOfOriginalRamp;
	ind_end = (rightP >= 0.2) ? ind_st_of_pil + (1 - rightP) * NumPointsOfOriginalRamp : ind_st_of_pil + (1 - 0.2) * NumPointsOfOriginalRamp;
	myflo stVP = vRamp[ind_st], endVP = vRamp[ind_end], delta;
	vSegRamp.resize((one_segment_width) * (1 - leftP - rightP));

	/* определяем приращение в одной точке пилы от предыдущей */
	myflo segment_length = endVP - stVP, segment_points_amount = abs(ind_st - ind_end), transform_factor = (myflo)one_segment_width / NumPointsOfOriginalRamp;
	delta = segment_length / (segment_points_amount * transform_factor);

#define iterator vector <myflo>::iterator
	function<iterator(iterator, iterator)> minormax;

	if (stVP > endVP) /* пила выглядит как '\' */
	{
		minormax = [](iterator begin, iterator end)
		{
			return max_element(begin, end);
		};
	}

	if (stVP < endVP) /* пила выглядит как '/' */
	{
		minormax = [](iterator begin, iterator end)
		{
			return min_element(begin, end);
		};
	}
#undef iterator
	
	/* Нужен держатель иначе дельта функция ругается на константный итератор */
	vector <myflo> vRampHolder(vRamp.begin() + ind_st_of_pil, vRamp.begin() + ind_st_of_pil + NumPointsOfOriginalRamp);
	int index_of_minormax = minormax(vRampHolder.begin(), vRampHolder.end()) - vRampHolder.begin();

	vSegRamp[0] = *(vRamp.begin() + ind_st_of_pil + index_of_minormax + leftP * NumPointsOfOriginalRamp);

	/* после этого пила всегда будет выглядеть как '/' */
	for (i = 1; i < vSegRamp.size(); ++i)
		vSegRamp[i] = vSegRamp[i - 1] + delta;

	/*if (vSegRamp.back() < vSegRamp.front())
		vectormult(vSegRamp, -1);*/
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

inline int match_Ramp_and_signal(_In_ const vector <myflo> & vRamp,
								 _In_ const vector <myflo> & vSignal,
								 _Inout_ vector <int> & vnRampPeaks)
{
	/* нужен держатель чтобы сделать индексы int с плавающей точкой, для точного приведения фазы */
	vector <myflo> IndHandler(vnRampPeaks.begin(), vnRampPeaks.end()); 
	if (vRamp.size() != vSignal.size())
		(vRamp.size() > vSignal.size()) ? vectordiv(IndHandler, (myflo)vRamp.size() / vSignal.size())
										: vectormult(IndHandler, (myflo)vSignal.size() / vRamp.size());

	/* возвращаем данные в основной вектор типа int */
	vnRampPeaks.assign(IndHandler.begin(), IndHandler.end());

	/*проверяем превышение крайнего индекса*/
	while (vnRampPeaks.back() > vSignal.size())
		vnRampPeaks.pop_back();

	/* окончательное приведение, во избежании небольшой погрешности */
	//vector <int> vnSignalPeaks;
	//easy_find_peaks(vSignal, one_segment_width * 0.9, -1, vnSignalPeaks);

	return 0;
}

int find_signal_and_make_Ramp(_In_ const vector <myflo> & vRamp,
							  _In_ const vector <myflo> & vSignal,
							  _Out_ vector <myflo> & vSegRamp,
							  _Out_ vector <int> & vStartSegIndxs)
{
	vector <myflo> vLnoise, vRnoise;

	int i = 0,
		j = 0,
		k = 0,
		segments_amount = 0,
		err = 0;

	/* слабо фильтруем для того чтобы убрать всплески */
	//vector <myflo> vSignal_copy(vSignal);
	//sg_smooth(vSignal_copy, vSignal, vSignal_copy.size() / 1000, 4);

	/* находим все участки на пиле при помощи пиков */
	vector <int> vnIndices;
	err = PeakFinder::findPeaks(vRamp, vnIndices, false, 1);
	if (vnIndices.size() <= 3) // 3 потому что дальше в функцию я передаю начало второго и он нужен целиком, тоесть от второго до третьего индекса
		return ERR_TooFewSegs;

	/* определяем длину участка изначальной пилы */
	int todel = (int)abs(vnIndices[2] - vnIndices[1]) % 100;
	if ((myflo)todel / 100 <= 0.5) 
		NumPointsOfOriginalRamp = abs(vnIndices[2] - vnIndices[1]) - todel;
	else 
		NumPointsOfOriginalRamp = abs(vnIndices[2] - vnIndices[1]) + 100 - todel;

	/* определяем длину участка сигнала */
	one_segment_width = (vRamp.size() != vSignal.size()) ? NumPointsOfOriginalRamp * ((myflo)vSignal.size() / vRamp.size()) :
		NumPointsOfOriginalRamp;

	/* проверка ошибок */
	if (one_segment_width == 0 
		|| NumPointsOfOriginalRamp == 0 
		|| is_invalid(one_segment_width) 
		|| is_invalid(NumPointsOfOriginalRamp)
		|| one_segment_width > vSignal.size()
		|| NumPointsOfOriginalRamp > vSignal.size())
		return ERR_BadSegsLength;

	/* линейно аппроксимируем пилу и задаем ей размер с учетом leftP и rightP */
	fit_linear_Ramp(vRamp, vnIndices[1] + 5, vSegRamp); // +5 потому что мы нашли пики пилы, а в этой функции пила строится от заданной точки, соответственно нужно спуститься вниз пилы
	if (vSegRamp.size() <= 0)
		return ERR_BadLinearRamp;

	/* приводим индексы начал отрезков к одной фазе с током */
	match_Ramp_and_signal(vRamp, vSignal, vnIndices);

	/* ВАЖНО!!!!!!!!!!!
	Мой метод предполагает, что либо в конце сигнала, либо в начале есть немного шума(чтобы взять его за образец и искать сигнал относительно него)
	*/
	vLnoise.resize(one_segment_width);
	vRnoise.resize(one_segment_width);

	/* заполняем шумовые вектора */
	noise_vecs(k, vLnoise, vRnoise, vSignal, vnIndices);

	vector <myflo> mnL_mxL_mnR_mxR(6);
	/*
	вектор mnL_mxL_mnR_mxR:
		-[0] минимум шума слева
		-[1] максимум шума справа
		-[2] минимум шума справа
		-[3] максимум шума справа
		-[4] минимум всего тока
		-[5] максимум всего тока
	*/

	/* заполняем вектор mnL_mxL_mnR_mxR */
	MinMaxNoise(vLnoise, vRnoise, mnL_mxL_mnR_mxR);

	/* если пользователь сам выбрал диапазон обработки -> не зайдем в цикл while */
	int temp = convert_time_to_pts(vSignal.size());
	bool dowhile = true;
	if (temp < 0) err = temp;
	if (temp == 0) dowhile = false;

	/* ОСНОВНОЙ ЦИКЛ ПОИСКА СИГНАЛА (ОЧИСТКИ ОТ ШУМА) */
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

		/* если количество всех отрезков меньше 0.1 от потенциально максимального, то их слишком мало -> пересчитываем (нужно увеличить количество) */
		if (segments_amount <= (int)(0.1f * vSignal.size() / one_segment_width))
		{
			/* уменьшаем коэффициенты чтобы понизить пороги отбора отрезков */
			koeff_MaxMin -= 0.3f;

			/* снова заполняем шумовые вектора */
			noise_vecs(k, vLnoise, vRnoise, vSignal, vnIndices);

			/* снова заполняем вектор mnL_mxL_mnR_mxR */
			MinMaxNoise(vLnoise, vRnoise, mnL_mxL_mnR_mxR);

			/* зануляем колво отрезков и начала/конца чтобы снова зайти в цикл while */
			segments_amount = 0;
			i_nach = 0;
			i_konec = 0;
		}
		/* если количество всех отрезков больше 0.8 от потенциально максимального, то их слишком много -> пересчитываем (нужно уменьшить количество) */
		if (segments_amount >= (int)(0.8f * vSignal.size() / one_segment_width))
		{
			/* увеличиваем коэффициенты чтобы повысить пороги отбора отрезков */
			koeff_MaxMin += 0.3f;

			/* снова заполняем шумовые вектора */
			noise_vecs(k, vLnoise, vRnoise, vSignal, vnIndices);

			/* снова заполняем вектор mnL_mxL_mnR_mxR */
			MinMaxNoise(vLnoise, vRnoise, mnL_mxL_mnR_mxR);

			/* зануляем колво отрезков и начала/конца чтобы снова зайти в цикл while */
			segments_amount = 0;
			i_nach = 0;
			i_konec = 0;
		}
	}

	/* проверка на случай если найдено слишком много отрезков, например сигнал - шум, а не полезный ток */
	if (dowhile && ( abs(i_konec - i_nach) / one_segment_width >= (int)(0.9f * vSignal.size() / (one_segment_width)) ))
		return ERR_TooManySegs;

	/* записываем индексы начал отрезков с учетом начала и конца сигнала */
	for (i = 0; i < vnIndices.size(); ++i)
	{
		if (vnIndices[i] > i_nach)
		{
			vStartSegIndxs.assign(vnIndices.begin() + i, vnIndices.begin() + i + abs(i_konec - i_nach) / one_segment_width);
			break;
		}
	}

	/* убираем лишнее с конца */
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

inline myflo find_floating_potential(_In_ const vector <myflo> & vRamp,
									 _In_ const vector <myflo> & vParams)
{
	myflo x = metod_Newton(vRamp.back(), vParams, fx_STEP);

	if (is_invalid(x) || !x)
	{
		x = metod_hord(vRamp[0], vRamp.back(), vParams, fx_STEP);
		if (is_invalid(x) || !x)
			x = vRamp.back();
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
	case 2:
		M = M_Ne;
		break;
	default:
		__assume(0);
	}

	return abs(ion_current) * 1E-6f / (0.52026f * S * ze * sqrt(kb * temperature * eVtoK / M)); //моя формула
	//return abs(ion_current) * 1E-6f / (S * ze * sqrt(2 * temperature * ze / M)); //формула сергея
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

	for (bool lflag = NULL, rflag = NULL; itl != vy.begin() && itr != vy.end() - 1;)
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

inline int make_one_segment(_In_ const int diagnostics,			 // diagnostics type (zond::0|setka::1|cilind::2)
							 _In_ const vector <myflo> & vRamp,		 // X data
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
			if (vRamp.size() != vSignal.size()
				|| vRamp.empty()
				|| vSignal.empty())
			{
				vcoeffs.resize(4);
				vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
				return ERR_BadSegInput;
			}

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vRamp.size());
			vfilt.resize(vRamp.size());
			vdiff.resize(vRamp.size());

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

			if (linfitP != 0)
			{
				vector <myflo> vX, vY;

				vX.assign(vRamp.begin(), vRamp.begin() + vRamp.size() * linfitP);
				vY.assign(vSignal.begin(), vSignal.begin() + vRamp.size() * linfitP);

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

			vParams[3] = 15;
			vParams[2] = abs((vSignal.back() - vParams[0] - vParams[1] * vRamp.back()) / (exp(vRamp.back() / vParams[3])));
			//vParams[2] = -vParams[0];

			if (linfitP != 0)
			{
				vector <myflo> vX, vY;

				vX.assign(vRamp.begin() + vRamp.size() * linfitP, vRamp.end());
				vY.assign(vSignal.begin() + vRamp.size() * linfitP, vSignal.end());

				levmarq(vX, vY, vParams, vFixed, Num_iter, fx_STEP);
				//LMfit(vX, vY, vParams, vFixed, Num_iter, fx_STEP, 1/*boost*/);
			}
			else
				levmarq(vRamp, vSignal, vParams, vFixed, Num_iter, fx_STEP);
				//LMfit(vRamp, vSignal, vParams, vFixed, Num_iter, fx_STEP, 1/*boost*/);

			if (vParams[3] < 0.0 
				|| vParams[3] > 100.0
				|| is_invalid(vParams[3]))
			{
				vector <myflo> vX, vY;

				vX.assign(vRamp.begin(), vRamp.begin() + vRamp.size() * 0.5);
				vY.assign(vSignal.begin(), vSignal.begin() + vRamp.size() * 0.5);

				vector <myflo> vAB = linear_fit(vX, vY);

				vParams[0] = vAB[0];
				vParams[1] = vAB[1];
				vParams[3] = 15;
				vParams[2] = abs((vSignal.back() - vParams[0] - vParams[1] * vRamp.back()) / (exp(vRamp.back() / vParams[3])));
				//vParams[2] = -vParams[0];
				
				vFixed[0] = vFixed[1] = true;

				vX.assign(vRamp.begin() + vRamp.size() * (1 - 0.5), vRamp.end());
				vY.assign(vSignal.begin() + vRamp.size() * (1 - 0.5), vSignal.end());

				levmarq(vX, vY, vParams, vFixed, max(Num_iter, 400), fx_STEP);
				//LMfit(vX, vY, vParams, vFixed, max(Num_iter, 400), fx_STEP, 1/*boost*/);

				if (vParams[3] < 0.0
					|| vParams[3] > 100.0
					|| is_invalid(vParams[3]))
				{
					vParams[3] = 15;
					vParams[2] = abs((vSignal.back() - vParams[0] - vParams[1] * vRamp.back()) / (exp(vRamp.back() / vParams[3])));
					//vParams[2] = - vAB[0] / 1E5;
				}
			}

			/* записываем в vres */
			for (i = 0; i < vRamp.size(); ++i)
				vres[i] = fx_STEP(vRamp[i], vParams);

			/* записываем в vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = find_floating_potential(vRamp, vParams);							// Floating potential
			vcoeffs[1] = find_ion_current(vParams[0], vParams[1], vcoeffs[0]);				// Saturate ion current
			vcoeffs[2] = vParams[3];														// Temperature
			vcoeffs[3] = find_density(vcoeffs[1], vcoeffs[2]);								// Density ne

			break;
		}
		case 1: // Setka
		{
			if (vRamp.size() != vSignal.size()
				|| vRamp.empty()
				|| vSignal.empty())
			{
				vcoeffs.resize(4);
				vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
				return ERR_BadSegInput;
			}

			int i = 0;

			vres.resize(vRamp.size());
			vfilt.resize(vRamp.size());
			vdiff.resize(vRamp.size());

			if (filtS == 0)
				filtS = 0.4f;

			/* фильтруем сигнал и записываем в vfilt */
			sg_smooth(vSignal, vfilt, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
			/* дифференцируем (первая производная) */
			diff(vfilt, vdiff, -1); // -1 -> чтобы сразу перевернуть
			/* фильтруем первую производную */
			sg_smooth(vdiff, vres, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
			vdiff = vres;

#ifdef GAUSSFITS
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			/* инициализируем прицельные параметры ф-ии Гаусса */
			GAUSS_InitParams(vRamp, vres, vParams);

			/* аппроксимируем производную Гауссом */
			levmarq(vRamp, vres, vParams, vFixed, Num_iter/* * 2*/, fx_GAUSS);

			if (vParams[2] < 0.0 
				|| vParams[2] > abs(vRamp.back() - vRamp[0]) / 2
				|| is_invalid(vParams[2]))
			{
				filtS = max(filtS, 0.8f);

				/* фильтруем сигнал и записываем в vfilt */
				sg_smooth(vSignal, vfilt, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				/* дифференцируем */
				diff(vfilt, vdiff, -1);
				/* фильтруем производную */
				sg_smooth(vdiff, vres, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vdiff = vres;

				/* инициализируем прицельные параметры ф-ии Гаусса */
				GAUSS_InitParams(vRamp, vres, vParams);

				levmarq(vRamp, vres, vParams, vFixed, min(Num_iter, 2), fx_GAUSS);

				if (vParams[2] < 0.0
					|| vParams[2] > 100.0
					|| is_invalid(vParams[2]))
					vParams[2] = 15;
			}

			/* записываем в vres */
			for (i = 0; i < vRamp.size(); ++i)
				vres[i] = fx_GAUSS(vRamp[i], vParams);

			/* записываем в vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = *max_element(vSignal.begin(), vSignal.end());						// Max Value
			vcoeffs[1] = abs(vParams[2]) * 0.5;												// Temp
			vcoeffs[2] = vRamp[max_element(vdiff.begin(), vdiff.end()) - vdiff.begin()];	// Peak Voltage
			vcoeffs[3] = find_width_at_half_maximum(vRamp, vdiff);							// Energy
#else
			vector <bool> vFixed = { false, false, false };
			vector <myflo> vParams(3);
			vector <myflo> vhandler(vRamp.size());

			int middlepoint_of_curvature = max_element(vdiff.begin(), vdiff.end()) - vdiff.begin();
			if (middlepoint_of_curvature >= vRamp.size() - 1)
			{
				vcoeffs.resize(4);
				vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
				return -1;
			}

			vhandler.assign(vdiff.begin() + middlepoint_of_curvature, vdiff.end());

			/* дифференцируем (вторая производная) */
			diff(vhandler, vres, -1); // -1 -> чтобы сразу перевернуть
			/* фильтруем вторую производную */
			sg_smooth(vres, vhandler, vRamp.size() * (filtS / 10), (5/*OriginLab poly order*/ - 1));
			
			/* отрезок аппроксимации будет длинной 5 полуширин изгиба */
			int endpoint_of_curvature = middlepoint_of_curvature + 5 * abs(max_element(vhandler.begin(), vhandler.end()) - vhandler.begin());

			if (endpoint_of_curvature < vRamp.size() && middlepoint_of_curvature < endpoint_of_curvature)
			{
				vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.begin() + endpoint_of_curvature);
				vres.assign(vfilt.begin() + middlepoint_of_curvature, vfilt.begin() + endpoint_of_curvature);
			}
			else
			{
				vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.end());
				vres.assign(vfilt.begin() + middlepoint_of_curvature, vfilt.end());
			}

			vParams[2] = -15;
			vParams[0] = *min_element(vres.begin(), vres.end());
			vParams[1] = /*log(*/abs((vres[0] - vParams[0]) / (exp(vhandler[0] / vParams[2])))/*) - vParams[2] * vhandler[0]*/;

			/* аппроксимируем выбранный участок экспонентой */
			levmarq(vhandler, vres, vParams, vFixed, max(Num_iter, 100), fx_EXP);

			vParams[2] = (vParams[2] > 0.0
						|| vParams[2] < -100.0
						|| is_invalid(vParams[2])) ? - 15 : vParams[2];
			vParams[0] = (is_invalid(vParams[0])) ? *min_element(vres.begin(), vres.end()) : vParams[0];
			vParams[1] = (vParams[1] < 0
						|| is_invalid(vParams[1])) ? log(abs((vres[0] - vParams[0]) / (exp(vhandler[0] / vParams[2])))) - vParams[2] * vhandler[0] : vParams[1];

			/* записываем в vres */
			vres.resize(vRamp.size());
			for (i = 0; i < vRamp.size(); ++i)
			{
				if (i >= middlepoint_of_curvature &&
					i <= middlepoint_of_curvature + vhandler.size())
					vres[i] = fx_EXP(vRamp[i], vParams);
				else
					vres[i] = NULL;
			}

			/* записываем в vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = *max_element(vSignal.begin(), vSignal.end());						// Max Value
			vcoeffs[1] = abs(vParams[2]);													// Temp
			vcoeffs[2] = vRamp[max_element(vdiff.begin(), vdiff.end()) - vdiff.begin()];	// Peak Voltage
			vcoeffs[3] = find_width_at_half_maximum(vRamp, vdiff);							// Energy
#endif
			break;
		}
		case 2: // Cilinder|Magnit
		{
			if (vRamp.size() != vSignal.size()
				|| vRamp.empty()
				|| vSignal.empty())
			{
				vcoeffs.resize(4);
				vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
				return ERR_BadSegInput;
			}

			int i = 0;

			vres.resize(vRamp.size());
			vfilt.resize(vRamp.size());
			vdiff.resize(vRamp.size());

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

#ifdef GAUSSFITC
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			/* инициализируем прицельные параметры ф-ии Гаусса */
			GAUSS_InitParams(vRamp, vSignal, vParams);

			/* аппроксимируем производную Гауссом */
			levmarq(vRamp, vSignal, vParams, vFixed, Num_iter/* * 2*/, fx_GAUSS);

			if (vParams[2] < 0.0 
				|| vParams[2] > abs(vRamp.back() - vRamp[0]) / 2
				|| is_invalid(vParams[2]))
			{
				filtS = max(filtS, 0.8f);

				/* фильтруем сигнал и записываем в vfilt */
				sg_smooth(vSignal, vres, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));

				/* инициализируем прицельные параметры ф-ии Гаусса */
				GAUSS_InitParams(vRamp, vres, vParams);

				levmarq(vRamp, vres, vParams, vFixed, min(Num_iter, 2), fx_GAUSS);

				if (vParams[2] < 0.0
					|| vParams[2] > 100.0
					|| is_invalid(vParams[2]))
					vParams[2] = 15;
			}

			/* записываем в vres */
			for (i = 0; i < vRamp.size(); ++i)
				vres[i] = fx_GAUSS(vRamp[i], vParams);

			/* записываем в vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = *max_element(vSignal.begin(), vSignal.end());					// Max Value
			vcoeffs[1] = abs(vParams[2]) * 0.5;											// Temp
			vcoeffs[2] = vRamp[max_element(vres.begin(), vres.end()) - vres.begin()];	// Peak Voltage
			vcoeffs[3] = find_width_at_half_maximum(vRamp, vres);						// Energy
#else
#ifndef debug
			vector <bool> vFixed = { false, false, false };
			vector <myflo> vParams(3);
			vector <myflo> vhandler(vRamp.size());

			/* дифференцируем (первая производная) */
			diff(vSignal, vhandler, -1); // -1 -> чтобы сразу перевернуть
			/* фильтруем первую производную */
			sg_smooth(vhandler, vres, vRamp.size() * ((filtS != 0) ? (filtS / 2) : (0.4 / 2)), (5/*OriginLab poly order*/ - 1));

			int middlepoint_of_curvature = max_element(vres.begin(), vres.end()) - vres.begin();
			if (middlepoint_of_curvature >= vRamp.size() - 1)
			{
				vcoeffs.resize(4);
				vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
				return -1;
			}

			/*vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.end());
			vres.assign(vSignal.begin() + middlepoint_of_curvature, vSignal.end());*/

			vhandler.assign(vres.begin() + middlepoint_of_curvature, vres.end());

			/* дифференцируем (вторая производная) */
			diff(vhandler, vres, -1); // -1 -> чтобы сразу перевернуть
			/* фильтруем вторую производную */
			sg_smooth(vres, vhandler, vRamp.size() * (filtS / 10), (5/*OriginLab poly order*/ - 1));

			int endpoint_of_curvature = middlepoint_of_curvature + 5 * abs(max_element(vhandler.begin(), vhandler.end()) - vhandler.begin());

			if (endpoint_of_curvature < vRamp.size() && middlepoint_of_curvature < endpoint_of_curvature)
			{
				vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.begin() + endpoint_of_curvature);
				vres.assign(vSignal.begin() + middlepoint_of_curvature, vSignal.begin() + endpoint_of_curvature);
			}
			else
			{
				vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.end());
				vres.assign(vSignal.begin() + middlepoint_of_curvature, vSignal.end());
			}

			vParams[2] = -abs((vhandler[0] - vhandler.back()) / log(vres[0] / vres.back())); /* из OroginLab */
			vParams[2] = (vParams[2] < -100 || is_invalid(vParams[2])) ? -100 : vParams[2];

			vParams[0] = *min_element(vres.begin(), vres.end());

			vParams[1] = log(abs((vres[0] - vParams[0]) / (exp(vhandler[0] / vParams[2])))) - vParams[2] * vhandler[0];
			vParams[1] = (vParams[1] < 0 || is_invalid(vParams[1])) ?
						log(abs(exp((vres.back() - vres[0]) / (exp((vhandler.back() - vhandler[0]) / vParams[2]) - 1)))) /* из OroginLab */
						: vParams[1];

			//MessageBoxA(NULL, ("vParams[0]: " + to_string(vParams[0])
			//			   + "\nvParams[1]: " + to_string(vParams[1])
			//			   + "\nvParams[2]: " + to_string(vParams[2])).c_str(), "Error!", MB_ICONINFORMATION | MB_OK);

			/* аппроксимируем выбранный участок экспонентой */
			levmarq(vhandler, vres, vParams, vFixed, max(Num_iter, 400), fx_EXP);

			//vParams[2] = (vParams[2] > 0.0
			//			|| vParams[2] < -100.0
			//			|| is_invalid(vParams[2])) ? -15 : vParams[2];
			vParams[0] = (is_invalid(vParams[0])) ? *min_element(vres.begin(), vres.end()) : vParams[0];
			vParams[1] = (vParams[1] < 0
						|| is_invalid(vParams[1])) ? log(abs((vres[0] - vParams[0]) / (exp(vhandler[0] / vParams[2])))) - vParams[2] * vhandler[0]
						: vParams[1];

			/* записываем в vres */
			vres.resize(vRamp.size());
			for (i = 0; i < vRamp.size(); ++i)
			{
				if (i >= middlepoint_of_curvature &&
					i <= middlepoint_of_curvature + vhandler.size())
					vres[i] = fx_EXP(vRamp[i], vParams);
				else
					vres[i] = NULL;
			}

			/* записываем в vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = *max_element(vSignal.begin(), vSignal.end());							// Max Value
			vcoeffs[1] = abs(vParams[2]);														// Temp
			vcoeffs[2] = vRamp[max_element(vSignal.begin(), vSignal.end()) - vSignal.begin()];	// Peak Voltage
			vcoeffs[3] = find_width_at_half_maximum(vRamp, vSignal);							// Energy
#else
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);
			vector <myflo> vhandler(vRamp.size());

			/* дифференцируем (первая производная) */
			diff(vSignal, vhandler, -1); // -1 -> чтобы сразу перевернуть
			/* фильтруем первую производную */
			sg_smooth(vhandler, vres, vRamp.size()* ((filtS != 0) ? (filtS / 2) : (0.4 / 2)), (5/*OriginLab poly order*/ - 1));

			int middlepoint_of_curvature = max_element(vres.begin(), vres.end()) - vres.begin();
			if (middlepoint_of_curvature >= vRamp.size() - 1)
			{
				vcoeffs.resize(4);
				vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
				return -1;
			}

			vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.end());
			vres.assign(vSignal.begin() + middlepoint_of_curvature, vSignal.end());

			//vhandler.assign(vres.begin() + middlepoint_of_curvature, vres.end());

			///* дифференцируем (вторая производная) */
			//diff(vhandler, vres, -1); // -1 -> чтобы сразу перевернуть
			///* фильтруем вторую производную */
			//sg_smooth(vres, vhandler, vRamp.size() * (filtS / 10), (5/*OriginLab poly order*/ - 1));

			//int endpoint_of_curvature = middlepoint_of_curvature + 5 * abs(max_element(vhandler.begin(), vhandler.end()) - vhandler.begin());

			//if (endpoint_of_curvature < vRamp.size() && middlepoint_of_curvature < endpoint_of_curvature)
			//{
			//	vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.begin() + endpoint_of_curvature);
			//	vres.assign(vSignal.begin() + middlepoint_of_curvature, vSignal.begin() + endpoint_of_curvature);
			//}
			//else
			//{
			//	vhandler.assign(vRamp.begin() + middlepoint_of_curvature, vRamp.end());
			//	vres.assign(vSignal.begin() + middlepoint_of_curvature, vSignal.end());
			//}

			{
				vector <myflo> vX, vY;

				vX.assign(vRamp.begin() + vRamp.size() * linfitP, vRamp.end());
				vY.assign(vSignal.begin() + vRamp.size() * linfitP, vSignal.end());

				vector <myflo> vAB = linear_fit(vX, vY);

				vParams[0] = vAB[0];
				vParams[1] = vAB[1];

				vFixed[0] = vFixed[1] = true;
			}

			vParams[3] = -15;
			vParams[2] = abs((vSignal.front() - vParams[0] - vParams[1] * vRamp.front()) / (exp(vRamp.front() / vParams[3])));

			/* аппроксимируем выбранный участок экспонентой */
			levmarq(vhandler, vres, vParams, vFixed, max(Num_iter, 400), fx_STEP);

			/* записываем в vres */
			vres.resize(vRamp.size());
			for (i = 0; i < vRamp.size(); ++i)
			{
				if (i >= middlepoint_of_curvature &&
					i <= middlepoint_of_curvature + vhandler.size())
					vres[i] = fx_STEP(vRamp[i], vParams);
				else
					vres[i] = NULL;
			}

			/* записываем в vcoeffs */
			vcoeffs.resize(4);
			vcoeffs[0] = *max_element(vSignal.begin(), vSignal.end());							// Max Value
			vcoeffs[1] = abs(vParams[3]);														// Temp
			vcoeffs[2] = vRamp[max_element(vSignal.begin(), vSignal.end()) - vSignal.begin()];	// Peak Voltage
			vcoeffs[3] = find_width_at_half_maximum(vRamp, vSignal);							// Energy
#endif
#endif
			break;
		}
		case 3:	// double langmuir probe
		{
			if (vRamp.size() != vSignal.size()
				|| vRamp.empty()
				|| vSignal.empty())
			{
				vcoeffs.resize(4);
				vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
				return ERR_BadSegInput;
			}

			int i = 0;
			vector <bool> vFixed = { false, false, false, false, false, false };
			vector <myflo> vParams(6);

			vres.resize(vRamp.size());
			vfilt.resize(vRamp.size());
			vdiff.resize(vRamp.size());

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vRamp.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

			PWL3_InitParams(vRamp, vSignal, vParams);

			levmarq(vRamp, vSignal, vParams, vFixed, Num_iter, fx_PWL3);

			/* записываем в vres */
			for (i = 0; i < vRamp.size(); ++i)
				vres[i] = fx_PWL3(vRamp[i], vParams);

			/* записываем в vcoeffs */
			myflo y1 = vParams[0] + vParams[1] * vParams[2],				// central PWL3 line offset 
				  Inegative;												// Saturate negative ion current
			vcoeffs.resize(4);
			vcoeffs[0] = vParams[2] - y1 / vParams[3];						// x of Y == 0 (OX crossing point)					
			Inegative = vParams[0] + vParams[1] * vcoeffs[0];				
			vcoeffs[1] = y1 + vParams[3] * (vParams[4] - vParams[2]) +
				vParams[5] * (vcoeffs[0] - vParams[4]);						// Saturate positive ion current
			vcoeffs[2] = abs((vParams[4] - vParams[2]) /
					(2 * (vcoeffs[1] - Inegative) / Inegative));			// Temperature
			vcoeffs[3] = find_density(vcoeffs[1], vcoeffs[2]);				// Density ne

			break;
		}
		default:
			__assume(0);
	}

	return 0;
}