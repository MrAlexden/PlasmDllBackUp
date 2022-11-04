#include "SubFuncs.h"

int NumPointsOfOriginalPila = 0,
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

	for (int i = vnIndices[0] + one_segment_width * leftP, k = 0; i < vnIndices[vnIndices.size() - 1]; i += ef_koef)
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
	int nach_ignored_gap = vnIndices[vnIndices.size() - 1] - (1 - rightP) * one_segment_width,
		konec_ignored_gap = vnIndices[vnIndices.size() - 1] - (1 + leftP) * one_segment_width;

	for (int i = vnIndices[vnIndices.size() - 1] - one_segment_width * rightP, counter = 0, k = 0; i > vnIndices[0]; i -= ef_koef)
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
	if (i < 0) i = ef_koef;

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
	if (i > vSignal.size() - 1) i = vSignal.size() - 1 - ef_koef;

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

	////////////////////////////////////// Определяем средний уровень шума //////////////////////////////////////
	if (mxL > mxR) mnL_mxL_mnR_mxR[5] = AverageSignal + abs(mxL - AverageSignal);
	else mnL_mxL_mnR_mxR[5] = AverageSignal + abs(mxR - AverageSignal);

	if (mnL < mnR) mnL_mxL_mnR_mxR[4] = AverageSignal - abs(AverageSignal - mnL);
	else mnL_mxL_mnR_mxR[4] = AverageSignal - abs(AverageSignal - mnR);
	////////////////////////////////////// Определяем средний уровень шума //////////////////////////////////////
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

	/* определяем приращение в одной точке пилы от предыдущей */
	myflo segment_length = endVP - stVP, segment_points_amount = abs(ind_st - ind_end), transform_factor = (myflo)one_segment_width / NumPointsOfOriginalPila;
	delta = segment_length / (segment_points_amount * transform_factor);

	vSegPila[0] = (leftP >= 0.2) ? stVP : vPila[ind_st_of_pil + leftP * NumPointsOfOriginalPila];
	for (i = 1; i < vSegPila.size(); ++i)
		vSegPila[i] = vSegPila[i - 1] + delta;
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

static inline int match_pila_and_signal(_In_ const vector <myflo> & vPila,
										_In_ const vector <myflo> & vSignal,
										_Inout_ vector <int> & vnPilaPeaks)
{
	/* нужен держатель чтобы сделать индексы int с плавающей точкой, для точного приведения фазы */
	vector <myflo> IndHandler(vnPilaPeaks.begin(), vnPilaPeaks.end()); 
	if (vPila.size() != vSignal.size())
		(vPila.size() > vSignal.size()) ? vectordiv(IndHandler, (vPila.size() / vSignal.size())) : vectormult(IndHandler, (vSignal.size() / vPila.size()));

	/*проверяем превышение крайнего индекса*/
	if (IndHandler[IndHandler.size() - 1] > vSignal.size())
		IndHandler.pop_back();

	/* возвращаем данные в основной вектор типа int */
	vnPilaPeaks.assign(IndHandler.begin(), IndHandler.end());

	/* окончательное приведение, во избежании небольшой погрешности */
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

	/* слабо фильтруем для того чтобы убрать всплески */
	//vector <myflo> vSignal_copy(vSignal);
	//sg_smooth(vSignal_copy, vSignal, vSignal_copy.size() / 1000, 4);

	/* находим все участки на пиле при помощи пиков */
	vector <int> vnIndices;
	err = PeakFinder::findPeaks(vPila, vnIndices, false, 1);
	if (vnIndices.size() <= 3) // 3 потому что дальше в функцию я передаю начало второго и он нужен целиком, тоесть от второго до третьего индекса
		return ERR_TooFewSegs;

	/* определяем длину участка изначальной пилы */
	int todel = (int)abs(vnIndices[2] - vnIndices[1]) % 100;
	if ((myflo)todel / 100 <= 0.5) 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) - todel;
	else 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) + 100 - todel;

	/* определяем длину участка сигнала */
	one_segment_width = (vPila.size() != vSignal.size()) ? NumPointsOfOriginalPila * ((myflo)vSignal.size() / vPila.size()) :
		NumPointsOfOriginalPila;

	/* проверка ошибок */
	if (one_segment_width == 0 
		|| NumPointsOfOriginalPila == 0 
		|| is_invalid(one_segment_width) 
		|| is_invalid(NumPointsOfOriginalPila)
		|| one_segment_width > vSignal.size()
		|| NumPointsOfOriginalPila > vSignal.size())
		return ERR_BadSegsLength;

	/* линейно аппроксимируем пилу и задаем ей размер с учетом leftP и rightP */
	fit_linear_pila(vPila, vnIndices[1] + 5, vSegPila); // +5 потому что мы нашли пики пилы, а в этой функции пила строится от заданной точки, соответственно нужно спуститься вниз пилы
	if (vSegPila.size() <= 0)
		return ERR_BadLinearPila;

	/* приводим индексы начал отрезков к одной фазе с током */
	match_pila_and_signal(vPila, vSignal, vnIndices);

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

inline myflo find_floating_potential(_In_ const vector <myflo> & vPila,
									 _In_ const vector <myflo> & vParams)
{
	myflo x = metod_Newton(vPila[vPila.size() - 1], vParams, fx_STEP);

	if (is_invalid(x) || !x)
	{
		x = metod_hord(vPila[0], vPila[vPila.size() - 1], vParams, fx_STEP);
		if (is_invalid(x) || !x)
			x = vPila[vPila.size() - 1];
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

	return ion_current * (10E-6f) / (0.52026f * S * (1.602f * 10E-19f) * sqrt((1.380649f * 10E-23f) * temperature * 11604.51812f / M)); //моя формула
	//return ion_current * (10E-6f) / (S * (1.602f * 10E-19f) * sqrt(2 * temperature * (1.602f * 10E-19f) / M)); //формула сергея
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

template <typename T>
inline int take_data_around_the_edges(_In_ const unsigned int size,
									  _In_ const vector <T> & vIN,
									  _Out_ vector <T> & vL,
									  _Out_ vector <T> & vR)
{
	vL.assign(vIN.begin(), vIN.begin() + size);
	vR.assign(vIN.end() - size, vIN.end());

	return 0;
}

int make_one_segment(_In_ const int diagnostics,					 // diagnostics type (zond::0|setka::1|cilind::2)
					 _In_ const vector <myflo> & vPila,		 // X data
					 _In_ vector <myflo> & vSignal,			 // Y data
					 _Out_ vector <myflo> & vres,			 // vector to be filled with the result
					 _Out_ vector <myflo> & vfilt,			 // vector to be filled with the filtration
					 _Inout_ vector <myflo> & vcoeffs)		 // additional coeffs/results vector
{
	myflo filtS = vcoeffs[0],
		linfitP = vcoeffs[1];

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

			if (vParams[3] < 0.0 || vParams[3] > 100.0)
			{
				vector <myflo> vX, vY;

				vX.assign(vPila.begin(), vPila.begin() + vPila.size() * 0.5);
				vY.assign(vSignal.begin(), vSignal.begin() + vPila.size() * 0.5);

				vector <myflo> vAB = linear_fit(vX, vY);

				vParams[0] = vAB[0];
				vParams[1] = vAB[1];
				vParams[2] = -vParams[0];
				vParams[3] = 10;

				//vFixed[0] = vFixed[1] = true;

				vX.assign(vPila.begin() + vPila.size() * (1 - 0.5), vPila.end());
				vY.assign(vSignal.begin() + vPila.size() * (1 - 0.5), vSignal.end());

				levmarq(vX, vY, vParams, vFixed, min(Num_iter, 2), fx_STEP);

				if (vParams[3] < 0.0 || vParams[3] > 100.0)
					vParams[3] = 10;
			}

			/* записываем в vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_STEP(vPila[i], vParams);

			vcoeffs.resize(4);
			vcoeffs[0] = find_floating_potential(vPila, vParams);					// Floating potential
			vcoeffs[1] = find_ion_current(vParams[0], vParams[1], vcoeffs[0]);		// Saturate ion current
			vcoeffs[2] = vParams[3];												// Temperature
			vcoeffs[3] = find_density(vcoeffs[1], vcoeffs[2]);						// Density ne

			break;
		}
		case 1: // Setka
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

			if (filtS == 0)
				filtS = 0.4f;
			/* фильтруем сигнал и записываем в vfilt */
			sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
			/* дифференцируем */
			diff(vfilt, vSignal, -1);
			/* фильтруем производную */
			sg_smooth(vSignal, vres, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));

			/* берем ветора с краев аппроксимируемых данных чтобы потом найти по ним y0 */
			vector <myflo> vL, vR;
			take_data_around_the_edges(vres.size() / 10, vres, vL, vR);

																											/* среднее между средним с краев */
			vParams[0] = (vSum(vL) / vL.size() + vSum(vR) / vR.size()) * 0.5;								// y0 - offset
																											/* значение X максимума */
			vParams[1] = vPila[max_element(vres.begin(), vres.end()) - vres.begin()];						// xc - center
																											/* вся ширина пилы делить на 10 */
			vParams[2] = abs(vPila[vPila.size() - 1] - vPila[0]) / 10;										// w - width
																											/* выражено из формулы yc = y0 + A / (w * sqrt(Pi / 2)) */
			vParams[3] = (*max_element(vres.begin(), vres.end()) - vParams[0]) * (vParams[2] * 1.25331414);	// A - area

			/* аппроксимируем производную Гауссом */
			levmarq(vPila, vres, vParams, vFixed, Num_iter/* * 2*/, fx_GAUSS);

			/* записываем в vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);

			vcoeffs.resize(4);
			vcoeffs[0] = vParams[0] + vParams[3] / (vParams[2] * 1.25331414/*sqrt(Pi / 2)*/);	// Max Value
			vcoeffs[1] = vParams[2];															// Temp
			vcoeffs[2] = vParams[1];															// Peak Voltage
			vcoeffs[3] = /*sqrt(log(4))*/1.17741001 * vParams[2];								// Energy

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

			if (filtS != 0)
			{
				sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
				vSignal = vfilt;
			}

			/* берем ветора с краев аппроксимируемых данных чтобы потом найти по ним y0 */
			vector <myflo> vL, vR;
			take_data_around_the_edges(vSignal.size() / 10, vSignal, vL, vR);

																													/* среднее между средним с краев */
			vParams[0] = (vSum(vL) / vL.size() + vSum(vR) / vR.size()) * 0.5;										// y0 - offset
																													/* значение X максимума */
			vParams[1] = vPila[max_element(vSignal.begin(), vSignal.end()) - vSignal.begin()];						// xc - center
																													/* вся ширина пилы делить на 10 */
			vParams[2] = abs(vPila[vPila.size() - 1] - vPila[0]) / 10;												// w - width
																													/* выражено из формулы yc = y0 + A / (w * sqrt(Pi / 2)) */
			vParams[3] = (*max_element(vSignal.begin(), vSignal.end()) - vParams[0]) * (vParams[2] * 1.25331414);	// A - area

			/* аппроксимируем производную Гауссом */
			levmarq(vPila, vSignal, vParams, vFixed, Num_iter/* * 2*/, fx_GAUSS);

			/* записываем в vres */
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);

			vcoeffs.resize(4);
			vcoeffs[0] = vParams[0] + vParams[3] / (vParams[2] * 1.25331414/*sqrt(Pi / 2)*/);	// Max Value
			vcoeffs[1] = vParams[2];															// Temp
			vcoeffs[2] = vParams[1];															// Peak Voltage
			vcoeffs[3] = /*sqrt(log(4))*/1.17741001 * vParams[2];								// Energy

			break;
		}
		default:
			__assume(0);
	}

	return 0;
}