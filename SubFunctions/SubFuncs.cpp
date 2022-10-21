#include "SubFuncs.h"

int NumPointsOfOriginalPila = 0,
	i_nach = 0,
	i_konec = 0;

myflo AverageSignal = 0,
	  koeff_MaxMin = 1.5;

inline myflo vSum(const vector <myflo> & v)
{
	myflo sum = 0.0;
	for (int i = 0; i < v.size(); ++i)
		sum += v[i];
	return sum;
}

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

inline int find_i_nach(const vector <myflo> & vSignal, const vector <myflo> & mnL_mxL_mnR_mxR, const vector <int> & vnIndices)
{
	myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05) ? 0.05 : otsech;
	int nach_ignored_gap = vnIndices[0] + (1 - otsech) * one_segment_width,
		konec_ignored_gap = vnIndices[0] + (1 + otsech) * one_segment_width;

	/*
	Шагаю не по 1, а по 2, так обработка в 2 раза быстрее(очевидно), а точность не теряется
	ТАКЖЕ ВАЖНО, ЧТО НАЧИНАЮ ИДТИ НЕМНОГО ОТСТУПЯ ОТ ВСПЛЕСКА (vnIndices[0] + one_segment_width*otsech)
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

inline int find_i_konec(const vector <myflo> & vSignal, const vector <myflo> & mnL_mxL_mnR_mxR, const vector <int> & vnIndices)
{
	myflo otsech = (leftP >= rightP) ? leftP : rightP;
	otsech = (otsech < 0.05) ? 0.05 : otsech;
	int nach_ignored_gap = vnIndices[vnIndices.size() - 1] - (1 - otsech) * one_segment_width,
		konec_ignored_gap = vnIndices[vnIndices.size() - 1] - (1 + otsech) * one_segment_width;

	/*
	Шагаю не по 1, а по 2, так обработка в 2 раза быстрее(очевидно), а точность не теряется
	ТАКЖЕ ВАЖНО, ЧТО НАЧИНАЮ ИДТИ НЕМНОГО ОТСТУПЯ ОТ ВСПЛЕСКА (vnIndices[vnIndices.GetSize()-1] - one_segment_width*otsech)
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

inline void noise_vecs(int k, vector <myflo> & vLnoise, vector <myflo> & vRnoise, const vector <myflo> & vSignal, const vector <int> & vnIndices)
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

inline void MinMaxNoise(const vector <myflo> & vLnoise, const vector <myflo> & vRnoise, vector <myflo> & mnL_mxL_mnR_mxR)
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

inline int convert_time_to_pts(int v_tok_size)
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
	{
		//MessageBoxA(NULL, "Less then 4 segments found", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_TooFewSegs;
	}

	/* определяем длину участка изначальной пилы */
	int todel = (int)abs(vnIndices[2] - vnIndices[1]) % 100;
	if ((myflo)todel / 100 <= 0.5) 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) - todel;
	else 
		NumPointsOfOriginalPila = abs(vnIndices[2] - vnIndices[1]) + 100 - todel;

	/* определяем длину участка сигнала */
	if (vPila.size() != vSignal.size())
		one_segment_width = NumPointsOfOriginalPila * ((myflo)vSignal.size() / vPila.size());
	else
		one_segment_width = NumPointsOfOriginalPila;

	/* проверка ошибок */
	if (one_segment_width == 0 
		|| NumPointsOfOriginalPila == 0 
		|| is_invalid(one_segment_width) 
		|| is_invalid(NumPointsOfOriginalPila)
		|| one_segment_width > vSignal.size()
		|| NumPointsOfOriginalPila > vSignal.size())
	{
		//MessageBoxA(NULL, "Error in finding segments length", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_BadSegsLength;
	}

	/* линейно аппроксимируем пилу и задаем ей размер с учетом leftP и rightP */
	fit_linear_pila(vPila, vnIndices[1] + 5, vSegPila); // +5 потому что мы нашли пики пилы, а в этой функции пила строится от заданной точки, соответственно нужно спуститься вниз пилы
	if (vSegPila.size() <= 0)
	{
		//MessageBoxA(NULL, "Error in pila linearizing", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_BadLinearPila;
	}

	/* приводим индексы начал отрезков к одной фазе с током */
	if (vPila.size() != vSignal.size())
		(vPila.size() > vSignal.size()) ? vectordiv<int, int>(vnIndices, (int)(vPila.size() / vSignal.size())) : vectormult<int, int>(vnIndices, (int)(vSignal.size() / vPila.size()));
	/* удаляем один элемент с конца (на всякий случай) и проверяем превышение крайнего индекса */
	vnIndices.pop_back();
	if (vnIndices[vnIndices.size() - 1] > vSignal.size())
		vectormult<int, int>(vnIndices, (int)(vPila.size() / vSignal.size()));

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
		k++;
		if (k > 5)
		{
			//MessageBoxA(NULL, "More than 5 attempts to find signal", "Error!", MB_ICONWARNING | MB_OK);
			return ERR_TooManyAttempts;
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
			//MessageBoxA(NULL, "Error in finding start|end of signal", "Error!", MB_ICONWARNING | MB_OK);
			return ERR_BadStartEnd;
		}

		/* если количество всех отрезков меньше 0.1 от потенциально максимального, то их слишком мало -> пересчитываем (нужно увеличить количество) */
		if (segments_amount <= (int)(0.1 * vSignal.size() / one_segment_width))
		{
			/* уменьшаем коэффициенты чтобы понизить пороги отбора отрезков */
			koeff_MaxMin -= 0.3;

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
		if (segments_amount >= (int)(0.8 * vSignal.size() / one_segment_width))
		{
			/* увеличиваем коэффициенты чтобы повысить пороги отбора отрезков */
			koeff_MaxMin += 0.3;

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
	if (dowhile && ( abs(i_konec - i_nach) / one_segment_width >= (int)(0.8 * vSignal.size() / (one_segment_width)) ))
	{
		//MessageBoxA(NULL, "Too many segments. Probably the signal is noise", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_TooManySegs;
	}

	/* ВАЖНО! ЗАДАЕМ РАЗМЕРА ВЕКТОРА НАЧАЛ ОТРЕЗКОВ */
	vStartSegIndxs.resize(abs(i_konec - i_nach) / one_segment_width);

	/* записываем индексы начал отрезков с учетом начала и конца сигнала */
	for (i = 0, j = 0; j < vStartSegIndxs.size() && i < vnIndices.size(); ++i)
	{
		if (vnIndices[i] > i_nach && vnIndices[i] < i_konec)
		{
			vStartSegIndxs[j] = vnIndices[i];
			j++;
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
	{
		//MessageBoxA(NULL, "No segments found", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_NoSegs;
	}

	return err;
}

inline myflo metod_hord(myflo x0, myflo x1, const vector <myflo> & vParams, myflo(*fx)(myflo, const vector <myflo> &))
{
	myflo x_cashe = 2 * x1;
	int n = 0;
	while (abs(x1 - x0) > 10E-10)
	{
		x0 = x1 - (x1 - x0) * fx(x1, vParams) / (fx(x1, vParams) - fx(x0, vParams));
		x1 = x0 - (x0 - x1) * fx(x0, vParams) / (fx(x0, vParams) - fx(x1, vParams));
		n++;
		if (x1 > x_cashe || x0 > x_cashe || n > 500)
		{
			x1 = x_cashe;
			break;
		}
	}
	return x1;
}

inline myflo metod_Newton(myflo x0, const vector <myflo> & vParams, myflo(*fx)(myflo, const vector <myflo> &))
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

inline myflo dens(const vector <myflo> & vPila, const vector <myflo> & vParams)
{
	myflo x = metod_Newton(vPila[vPila.size() - 1], vParams, fx_STEP), ans = 0.0, M = 0.0;
	
	if (is_invalid(x) || !x)
	{
		x = metod_hord(vPila[0], vPila[vPila.size() - 1], vParams, fx_STEP);
		if (is_invalid(x) || !x)
			x = vPila[vPila.size() - 1];
	}
	ans = abs(vParams[0] + vParams[1] * x);
	if (is_invalid(ans) || !ans)
		ans = vParams[0];

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

	return (myflo)ans * (10E-6) / (0.52026 * S * (1.602 * 10E-19) * sqrt((1.380649 * 10E-23) * (myflo)vParams[3] * 11604.51812 / M)); //моя формула
	//return (myflo)ans * (10E-6) / (S * (1.602 * 10E-19) * sqrt(2 * (myflo)vParams[3] * (1.602 * 10E-19) / M)); //формула сергея
}

int make_one_segment(_In_ int diagnostics,					 // diagnostics type (zond::0|setka::1|cilind::2)
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
			{
				//MessageBoxA(NULL, "Input segment's values error", "Error!", MB_ICONWARNING | MB_OK);
				return ERR_BadSegInput;
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
				vParams[0] = vSignal[0];
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

			vcoeffs = vParams;

			/* записываем в vres */
#pragma omp parallel for schedule(static, 1) 
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_STEP(vPila[i], vParams);

			vcoeffs.push_back(dens(vPila, vParams));

			break;
		}
		case 1: // Setka
		{
			if (vPila.size() != vSignal.size()
				|| vPila.empty()
				|| vSignal.empty())
			{
				//MessageBoxA(NULL, "Input segment's values error", "Error!", MB_ICONWARNING | MB_OK);
				return ERR_BadSegInput;
			}

			int i = 0;
			vector <bool> vFixed = { false, false, false, false };
			vector <myflo> vParams(4);

			vres.resize(vPila.size());
			vfilt.resize(vPila.size());

			if (filtS == 0)
				filtS = 0.4;
			/* фильтруем сигнал и записываем в vfilt */
			sg_smooth(vSignal, vfilt, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
			/* дифференцируем */
			diff(vfilt, vSignal, -1);
			/* фильтруем производную */
			sg_smooth(vSignal, vres, vPila.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));

			vParams[0] = 0;													// y0 - offset
			vParams[1] = vPila[vPila.size() * 0.5];							// xc - center
			vParams[2] = abs(vPila[vPila.size() - 1] - vPila[0]) / 10;		// w - width
			vParams[3] = 1;													// A - area

			/* аппроксимируем производную Гауссом */
			levmarq(vPila, vres, vParams, vFixed, Num_iter * 2, fx_GAUSS);

			vcoeffs.resize(4);
			vcoeffs[0] = vParams[0] + vParams[3] / (vParams[2] * sqrt(Pi / 2));		// Max Value
			vcoeffs[1] = vParams[2];												// Temp
			vcoeffs[2] = vParams[1];												// Peak Voltage
			vcoeffs[3] = sqrt(log(4)) * vParams[2];									// Energy

			/* записываем в vres */
#pragma omp parallel for schedule(static, 1) 
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);

			break;
		}
		case 2: // Cilinder|Magnit
		{
			if (vPila.size() != vSignal.size()
				|| vPila.empty()
				|| vSignal.empty())
			{
				//MessageBoxA(NULL, "Input segment's values error", "Error!", MB_ICONWARNING | MB_OK);
				return ERR_BadSegInput;
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
			vParams[2] = abs(vPila[vPila.size() - 1] - vPila[0]) / 10;		// w - width
			vParams[3] = 1;													// A - area

			/* аппроксимируем производную Гауссом */
			levmarq(vPila, vSignal, vParams, vFixed, Num_iter * 2, fx_GAUSS);

			vcoeffs.resize(4);
			vcoeffs[0] = vParams[0] + vParams[3] / (vParams[2] * sqrt(Pi / 2));		// Max Value
			vcoeffs[1] = vParams[2];												// Temp
			vcoeffs[2] = vParams[1];												// Peak Voltage
			vcoeffs[3] = sqrt(log(4)) * vParams[2];									// Energy

			/* записываем в vres */
#pragma omp parallel for schedule(static, 1) 
			for (i = 0; i < vPila.size(); ++i)
				vres[i] = fx_GAUSS(vPila[i], vParams);

			break;
		}
		default:
			__assume(0);
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

bool is_signalpeakslookingdown(_In_ const vector <myflo> & v)
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