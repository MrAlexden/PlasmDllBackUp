#include "Zond.h"

extern "C" __declspec(dllexport) Plasma_proc_result * Zond(vector <myflo> vPila, vector <myflo> vSignal, vector <myflo> AdditionalData)
{
	if (vPila.size() == 0
		|| vPila.empty()
		|| vSignal.size() == 0
		|| vSignal.empty()
		|| AdditionalData.size() == 0
		|| AdditionalData.empty())
	{
		MessageBoxA(NULL, "Corrupted input vectors", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;
	}
	if (AdditionalData[0] == 0 || is_invalid(AdditionalData[0]) ||
		(int)AdditionalData[7] == 0 || is_invalid(AdditionalData[7]) ||
		(int)AdditionalData[8] == 0 || is_invalid(AdditionalData[8]) ||
		(int)AdditionalData[9] == 0 || is_invalid(AdditionalData[9]) ||
		(int)AdditionalData[11] == 0 || is_invalid(AdditionalData[11]))
	{
		MessageBoxA(NULL, "Input data error, some values equals 0", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;
	}
	if (AdditionalData[3] >= 0.5 || AdditionalData[3] < 0.0)
	{
		MessageBoxA(NULL, "Сut-off points on the left value must be > 0.0 and < 0.5", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;
	}
	if (AdditionalData[4] >= 0.5 || AdditionalData[4] < 0.0)
	{
		MessageBoxA(NULL, "Сut-off points on the right value must be > 0.0 and < 0.5", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;
	}

	vector <myflo> vSegPila;
	vector <int> vStartSegIndxs;

	int	numSegments = 0,	// количество отрезков в импульсе
		resistance = 0,
		coefPila = 0,
		dimension = 6;		// количество столбиков parameters

	S = AdditionalData[0];						    // площадь поверхности зонда
	st_time_end_time[0] = AdditionalData[1];		// время начала обработки (если этот параметр не выбран: -1)
	st_time_end_time[1] = AdditionalData[2];		// время конца обработки (если этот параметр не выбран: -1)
	leftP = AdditionalData[3];					    // часть точек отсечки слева
	rightP = AdditionalData[4];					    // часть точек отсечки справа
	linfitP = AdditionalData[5];					// часть точек линейно аппроксимации
	filtS = AdditionalData[6];					    // часть точек фильтрации сигнала
	freqP = (int)AdditionalData[7];					// частота пилы
	resistance = (int)AdditionalData[8];			// сопротивление на зонде
	coefPila = (int)AdditionalData[9];				// коэффициент усиления пилы
	fuel = (int)AdditionalData[10];					// рабочее вещество (He::0|Ar::1)
	Num_iter = (int)AdditionalData[11];				// количество итераций аппроксимации(сильно влияет на скорость работы программы)

	/* домножаем пилу на коэффициент усиления */
	vectormult(vPila, coefPila);
	/* переворачиваем ток чтобы смотрел вверх(если нужно), и делим на сопротивление */
	vectordiv(vSignal, -resistance);

	if (is_invalid(vPila[0])
		|| is_invalid(vPila[vPila.size() - 1])
		|| is_invalid(vSignal[0])
		|| is_invalid(vSignal[vSignal.size() - 1]))
	{
		MessageBoxA(NULL, "Error after Pila|Signal factorizing", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;
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
		MessageBoxA(NULL, "Error after noise extracting", "Error!", MB_ICONWARNING | MB_OK);
		return nullptr;
	}

	myflo *** fdata = new myflo ** [5];
	for (int i = 0; i < 3; ++i)
	{
		fdata[i] = new myflo * [numSegments + 1];
		for (int j = 0; j < numSegments + 1; ++j)
		{
			fdata[i][j] = new myflo[vSegPila.size()];
		}
	}
	fdata[3] = new myflo * [numSegments];
	for (int j = 0; j < numSegments; ++j)
	{
		fdata[3][j] = new myflo[dimension];
	}
	fdata[4] = new myflo * [1];
	fdata[4][0] = new myflo[3];

	/*=============================================================================================
	-- fdata[0] - оригинальные данные
		--[i][j]
			i::столбики(количество отрезков + 1, т.к. первая пила)
			j::строки(размер отрезка)
	-- fdata[1] - обработанные данные(результат)
		--[i][j]
			i::столбики(количество отрезков + 1, т.к. первая пила)
			j::строки(размер отрезка)
	-- fdata[2] - фильтрованный сигнал
		--[i][j]
			i::столбики(количество отрезков + 1, т.к. первая пила)
			j::строки(размер отрезка)
	-- fdata[3] - массив с параметрами
		--[i][j]
			i::строки(количество отрезков)
			j::столбики, их всегда 6 для зонда
				--1. время
				--2. A
				--3. B
				--4. C
				--5. D(температура)
				--6. Плотность
	-- fdata[4] - массив из одного элемента из трех элементов
		--[0][i] 
			i::3 размера данных
			--1. столбики для fdata[0],fdata[1],fdata[2](количество отрезков + 1, т.к. первая пила)
				и количество строк для fdata[3](нужно отнять 1, т.к. нет пилы)
			--2. размер отрезка для fdata[0],fdata[1],fdata[2]
			--3. столбики для fdata[3], их всегда 6 для зонда
	=============================================================================================*/

	fdata[4][0][0] = (myflo)numSegments;
	fdata[4][0][1] = (myflo)vSegPila.size();
	fdata[4][0][2] = (myflo)dimension;

	memcpy(fdata[0][0], vSegPila.data(), sizeof myflo * vSegPila.size());
	memcpy(fdata[1][0], vSegPila.data(), sizeof myflo * vSegPila.size());
	memcpy(fdata[2][0], vSegPila.data(), sizeof myflo * vSegPila.size());

#pragma omp parallel for schedule(static, 1) 
	for (int segnum = 1; segnum < numSegments + 1; ++segnum)
	{
		vector <myflo> vY, vres, vfilt, vcoeffs = { filtS , linfitP };

		vY.assign(vSignal.begin() + vStartSegIndxs[segnum - 1] + one_segment_width * leftP,
			vSignal.begin() + vStartSegIndxs[segnum - 1] + one_segment_width * leftP + vSegPila.size());

		memcpy(fdata[0][segnum], vY.data(), sizeof myflo * vSegPila.size());

		if (make_one_segment(0, vSegPila, vY, vres, vfilt, vcoeffs) == -1)
			continue;

		memcpy(fdata[1][segnum], vres.data(), sizeof myflo * vSegPila.size());
		memcpy(fdata[2][segnum], vfilt.data(), sizeof myflo * vSegPila.size());

		fdata[3][segnum - 1][0] = vStartSegIndxs[segnum - 1] * (1.0 / (one_segment_width * freqP));
		for (int i = 1, j = 0; i < dimension; ++i, ++j)
			fdata[3][segnum - 1][i] = vcoeffs[j];
	}

	return fdata;
}