#include "Cilinder.h"

extern "C" __declspec(dllexport) Plasma_proc_result * Cilinder(vector <myflo> vPila, vector <myflo> vSignal, vector <myflo> AdditionalData)
{
	if (vPila.size() == 0
		|| vPila.empty()
		|| vSignal.size() == 0
		|| vSignal.empty()
		|| AdditionalData.size() == 0
		|| AdditionalData.empty())
	{
		//MessageBoxA(NULL, "Corrupted input vectors", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;//ERR_BadInputVecs
	}
	if ((int)AdditionalData[7] == 0 || is_invalid(AdditionalData[7]) ||
		(int)AdditionalData[8] == 0 || is_invalid(AdditionalData[8]) ||
		(int)AdditionalData[9] == 0 || is_invalid(AdditionalData[9]) ||
		(int)AdditionalData[11] == 0 || is_invalid(AdditionalData[11]))
	{
		//MessageBoxA(NULL, "Input data error, some values equals 0", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;//ZeroInputVals
	}
	if (AdditionalData[3] >= 0.5 || AdditionalData[3] < 0.0)
	{
		//MessageBoxA(NULL, "Сut-off points on the left value must be > 0.0 and < 0.5", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;//BadCutOffLeft
	}
	if (AdditionalData[4] >= 0.5 || AdditionalData[4] < 0.0)
	{
		//MessageBoxA(NULL, "Сut-off points on the right value must be > 0.0 and < 0.5", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;//BadCutOffRight
	}

	vector <myflo> vSegPila;
	vector <int> vStartSegIndxs;

	int	numSegments = 0,	// количество отрезков в импульсе
		resistance = 0,
		coefPila = 0,
		dimension = 5;		// количество столбиков parameters

	S = 0.0;										// площадь поверхности зонда
	st_time_end_time[0] = AdditionalData[1];		// время начала обработки (если этот параметр не выбран: -1)
	st_time_end_time[1] = AdditionalData[2];		// время конца обработки (если этот параметр не выбран: -1)
	leftP = AdditionalData[3];					    // часть точек отсечки слева
	rightP = AdditionalData[4];					    // часть точек отсечки справа
	linfitP = 0.0;									// часть точек линейно аппроксимации
	filtS = AdditionalData[6];					    // часть точек фильтрации сигнала
	freqP = (int)AdditionalData[7];					// частота пилы
	resistance = (int)AdditionalData[8];			// сопротивление на цилиндре|магните
	coefPila = (int)AdditionalData[9];				// коэффициент усиления пилы
	fuel = (int)AdditionalData[10];					// рабочее вещество (He::0|Ar::1)
	Num_iter = (int)AdditionalData[11];				// количество итераций аппроксимации(сильно влияет на скорость работы программы)

	/* домножаем пилу на коэффициент усиления */
	vectormult(vPila, coefPila);
	/* переворачиваем ток чтобы смотрел вверх(если нужно), и делим на сопротивление */
	if (is_signalpeakslookingdown(vSignal)) 
		vectordiv(vSignal, -resistance);
	else
		vectordiv(vSignal, resistance);

	if (is_invalid(vPila[0])
		|| is_invalid(vPila[vPila.size() - 1])
		|| is_invalid(vSignal[0])
		|| is_invalid(vSignal[vSignal.size() - 1]))
	{
		//MessageBoxA(NULL, "Error after Pila|Signal factorizing", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;//BadFactorizing
	}

	if (find_signal_and_make_pila(vPila, vSignal, vSegPila, vStartSegIndxs) == -1) return nullptr;
	numSegments = vStartSegIndxs.size();

	if (vStartSegIndxs.size() == 0
		|| vStartSegIndxs.empty()
		|| is_invalid(vSegPila[0])
		|| is_invalid(vSegPila[vSegPila.size() - 1])
		|| is_invalid(vStartSegIndxs[0])
		|| is_invalid(vStartSegIndxs[vStartSegIndxs.size() - 1]))
	{
		//MessageBoxA(NULL, "Error after noise extracting", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;//BadNoise
	}

	Plasma_proc_result* fdata = new Plasma_proc_result(numSegments, vSegPila.size(), dimension);

	fdata->SetPila(vSegPila);

#pragma omp parallel for schedule(static, 1) 
	for (int segnum = 0; segnum < numSegments; ++segnum)
	{
		vector <myflo> vY, vres, vfilt, vcoeffs = { filtS , linfitP };

		vY.assign(vSignal.begin() + vStartSegIndxs[segnum] + one_segment_width * leftP,
			vSignal.begin() + vStartSegIndxs[segnum] + one_segment_width * leftP + vSegPila.size());

		fdata->SetOriginSegment(vY, segnum);

		if (make_one_segment(2, vSegPila, vY, vres, vfilt, vcoeffs) == -1)
			continue;

		vcoeffs.insert(vcoeffs.begin(), vStartSegIndxs[segnum] * (1.0 / (one_segment_width * freqP)));

		fdata->SetFiltedSegment(vfilt, segnum);
		fdata->SetApproxSegment(vres, segnum);
		fdata->SetParamsSegment(vcoeffs, segnum);
	}

	return fdata;
}