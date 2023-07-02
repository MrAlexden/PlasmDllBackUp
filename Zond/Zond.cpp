#include "Zond.h"

int Zond(_In_ vector <myflo> & vRamp,
		 _In_ vector <myflo> & vSignal, 
		 _In_ const vector <myflo> & AdditionalData,
		 _Out_ Plasma_proc_result <myflo> & fdata)
{
	if (vRamp.size() == 0
		|| vRamp.empty()
		|| vSignal.size() == 0
		|| vSignal.empty()
		|| AdditionalData.size() == 0
		|| AdditionalData.empty())
		return ERR_BadInputVecs;
	if (AdditionalData[0] == 0 || is_invalid(AdditionalData[0]) ||
		(int)AdditionalData[7] == 0 || is_invalid(AdditionalData[7]) ||
		(int)AdditionalData[8] == 0 || is_invalid(AdditionalData[8]) ||
		(int)AdditionalData[9] == 0 || is_invalid(AdditionalData[9]) ||
		(int)AdditionalData[11] == 0 || is_invalid(AdditionalData[11]))
		return ERR_ZeroInputVals;
	if (AdditionalData[3] + AdditionalData[4] > 0.9
		|| AdditionalData[3] < 0.0
		|| AdditionalData[4] < 0.0)
		return ERR_BadCutOff;
	if (AdditionalData[5] < 0 || AdditionalData[5] > 0.9)
		return ERR_BadLinFit;

	vector <myflo> vSegRamp;
	vector <int> vStartSegIndxs;

	int	numSegments = 0,	// количество отрезков в импульсе
		resistance = 0,
		coefRamp = 0,
		dimension = 5,		// количество столбиков parameters
		err = 0;

	ThreadArgs args;
	HANDLE hThread;

	S = AdditionalData[0];						    // площадь поверхности зонда
	st_time_end_time[0] = AdditionalData[1];		// врем€ начала обработки (если этот параметр не выбран: -1)
	st_time_end_time[1] = AdditionalData[2];		// врем€ конца обработки (если этот параметр не выбран: -1)
	leftP = AdditionalData[3];					    // часть точек отсечки слева
	rightP = AdditionalData[4];					    // часть точек отсечки справа
	linfitP = AdditionalData[5];					// часть точек линейно аппроксимации
	filtS = AdditionalData[6];					    // часть точек фильтрации сигнала
	freqP = (int)AdditionalData[7];					// частота пилы
	resistance = (int)AdditionalData[8];			// сопротивление на зонде
	coefRamp = (int)AdditionalData[9];				// коэффициент усилени€ пилы
	fuel = (int)AdditionalData[10];					// рабочее вещество (He::0|Ar::1|Ne::2)
	Num_iter = (int)AdditionalData[11];				// количество итераций аппроксимации(сильно вли€ет на скорость работы программы)

	/* домножаем пилу на коэффициент усилени€ */
	thread T1(vectormult<myflo, int>, ref(vRamp), coefRamp);
	//vectormult(vRamp, coefRamp);
	/* переворачиваем ток чтобы смотрел вверх(если нужно), и делим на сопротивление */
	thread T2(vectordiv<myflo, int>, ref(vSignal), -resistance);
	//vectordiv(vSignal, resistance);

	T1.join();
	T2.join();

	if (is_invalid(vRamp.at(0))
		|| is_invalid(vRamp.at(vRamp.size() - 1))
		|| is_invalid(vSignal.at(0))
		|| is_invalid(vSignal.at(vSignal.size() - 1)))
		return ERR_BadFactorizing;

	ERR(find_signal_and_make_Ramp(vRamp, vSignal, vSegRamp, vStartSegIndxs));
	numSegments = vStartSegIndxs.size();

#ifdef KLUDGE
	{/* ¬–≈ћ≈ЌЌџ…  ќ—“џЋ№ */
		vector <myflo> IndHandler(vStartSegIndxs.begin(), vStartSegIndxs.end());
		vectormult(IndHandler, 1.001f);
		vStartSegIndxs.assign(IndHandler.begin(), IndHandler.end());
		/*провер€ем превышение крайнего индекса*/
		while (vStartSegIndxs.back() > vSignal.size())
			vStartSegIndxs.pop_back();
	}/* ¬–≈ћ≈ЌЌџ…  ќ—“џЋ№ */
#endif

	if (vStartSegIndxs.size() == 0
		|| vStartSegIndxs.empty()
		|| is_invalid(vSegRamp.at(0))
		|| is_invalid(vSegRamp.at(vSegRamp.size() - 1))
		|| is_invalid(vStartSegIndxs.at(0))
		|| is_invalid(vStartSegIndxs.at(vStartSegIndxs.size() - 1)))
		return ERR_BadNoise;

	fdata.SetSegmentsNumber(numSegments);
	fdata.SetSegmentsSize(vSegRamp.size());
	fdata.SetParamsNumber(dimension);
	fdata.SetRamp(vSegRamp);

	hThread = CreateThread(NULL, NULL, (LPTHREAD_START_ROUTINE)&DialogBoxParamWrapper, &args, NULL, NULL);
	if (hThread)
	{
		WaitForSingleObject(hThread, 10);
		SendMessage(GetDlgItem(mywindow, IDC_PROGRESS1), PBM_SETRANGE, 0, MAKELPARAM(0, numSegments - 1));
	}

	for (int segnum = 0; segnum < numSegments; ++segnum)
	{	
		if (segnum % 4 == 0) // отправл€ю мэсседж раз в 4 итерации чтобы меньше нагружать
			SendMessage(GetDlgItem(mywindow, IDC_PROGRESS1), PBM_SETPOS, segnum, NULL);
		wstring wstr = L"In Progress... " + to_wstring(int(((myflo)segnum / numSegments) * 100)) + L" %";
		SetDlgItemText(mywindow, IDC_STATIC, wstr.c_str());

		vector <myflo> vY, 
					   vres, 
					   vfilt, 
					   vdiff, 
					   vcoeffs = { S, linfitP, filtS , (myflo)fuel, (myflo)Num_iter };

		vY.assign(vSignal.begin() + vStartSegIndxs.at(segnum) + one_segment_width * leftP,
			vSignal.begin() + vStartSegIndxs.at(segnum) + one_segment_width * leftP + vSegRamp.size());

		fdata.SetOriginSegment(vY, segnum);

		/*if (make_one_segment(0, vSegRamp, vY, vres, vfilt, vdiff, vcoeffs) < 0)
			continue;*/
		make_one_segment(0, vSegRamp, vY, vres, vfilt, vdiff, vcoeffs);

		vcoeffs.insert(vcoeffs.begin(), vStartSegIndxs.at(segnum) * (1.0 / (one_segment_width * freqP)));
	
		fdata.SetApproxSegment(vres, segnum);
		fdata.SetFiltedSegment(vfilt, segnum);
		fdata.SetDiffedSegment(vdiff, segnum);
		fdata.SetParamsSegment(vcoeffs, segnum);
	}	
	
	if (hThread)
	{
		TerminateThread(hThread, -1);
		CloseHandle(hThread);
	}

Error:
	return err;
}