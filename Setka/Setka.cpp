#include "Setka.h"

int Setka(_In_ vector <myflo> & vPila,
		  _In_ vector <myflo> & vSignal,
		  _In_ const vector <myflo> & AdditionalData,
		  _Out_ Plasma_proc_result <myflo> & fdata)
{
	if (vPila.size() == 0
		|| vPila.empty()
		|| vSignal.size() == 0
		|| vSignal.empty()
		|| AdditionalData.size() == 0
		|| AdditionalData.empty())
		return ERR_BadInputVecs;
	if ((int)AdditionalData[7] == 0 || is_invalid(AdditionalData[7]) ||
		(int)AdditionalData[8] == 0 || is_invalid(AdditionalData[8]) ||
		(int)AdditionalData[9] == 0 || is_invalid(AdditionalData[9]) ||
		(int)AdditionalData[11] == 0 || is_invalid(AdditionalData[11]))
		return ERR_ZeroInputVals;
	if (AdditionalData[3] + AdditionalData[4] > 0.9
		|| AdditionalData[3] < 0.0
		|| AdditionalData[4] < 0.0)
		return ERR_BadCutOff;

	vector <myflo> vSegPila;
	vector <int> vStartSegIndxs;

	int	numSegments = 0,	// ���������� �������� � ��������
		resistance = 0,
		coefPila = 0,
		dimension = 5,		// ���������� ��������� parameters
		err = 0;

	ThreadArgs args;
	HANDLE hThread;

	S = 0.0;										// ������� ����������� �����
	st_time_end_time[0] = AdditionalData[1];		// ����� ������ ��������� (���� ���� �������� �� ������: -1)
	st_time_end_time[1] = AdditionalData[2];		// ����� ����� ��������� (���� ���� �������� �� ������: -1)
	leftP = AdditionalData[3];					    // ����� ����� ������� �����
	rightP = AdditionalData[4];					    // ����� ����� ������� ������
	linfitP = 0.0;									// ����� ����� ������� �������������
	filtS = AdditionalData[6];					    // ����� ����� ���������� �������
	freqP = (int)AdditionalData[7];					// ������� ����
	resistance = (int)AdditionalData[8];			// ������������� �� �����
	coefPila = (int)AdditionalData[9];				// ����������� �������� ����
	fuel = (int)AdditionalData[10];					// ������� �������� (He::0|Ar::1|Ne::2)
	Num_iter = (int)AdditionalData[11];				// ���������� �������� �������������(������ ������ �� �������� ������ ���������)

	/* ��������� ���� �� ����������� �������� */
	thread T1(vectormult<myflo, int>, ref(vPila), coefPila);
	//vectormult(vPila, coefPila);
	/* �������������� ��� ����� ������� �����(���� �����), � ����� �� ������������� */
	thread T2(vectordiv<myflo, int>, ref(vSignal), resistance);
	//vectordiv(vSignal, resistance);

	T1.join();
	T2.join();

	if (is_invalid(vPila.at(0))
		|| is_invalid(vPila.at(vPila.size() - 1))
		|| is_invalid(vSignal.at(0))
		|| is_invalid(vSignal.at(vSignal.size() - 1)))
		return ERR_BadFactorizing;

	ERR(find_signal_and_make_pila(vPila, vSignal, vSegPila, vStartSegIndxs));
	numSegments = vStartSegIndxs.size();

#ifdef KLUDGE
	{/* ��������� ������� */
		vector <myflo> IndHandler(vStartSegIndxs.begin(), vStartSegIndxs.end());
		vectormult(IndHandler, 1.001f);
		vStartSegIndxs.assign(IndHandler.begin(), IndHandler.end());
		/*��������� ���������� �������� �������*/
		while (vStartSegIndxs.back() > vSignal.size())
			vStartSegIndxs.pop_back();
	}/* ��������� ������� */
#endif

	if (vStartSegIndxs.size() == 0
		|| vStartSegIndxs.empty()
		|| is_invalid(vSegPila.at(0))
		|| is_invalid(vSegPila.at(vSegPila.size() - 1))
		|| is_invalid(vStartSegIndxs.at(0))
		|| is_invalid(vStartSegIndxs.at(vStartSegIndxs.size() - 1)))
		return ERR_BadNoise;

	fdata.SetSegmentsNumber(numSegments);
	fdata.SetSegmentsSize(vSegPila.size());
	fdata.SetParamsNumber(dimension);
	fdata.SetPila(vSegPila);

	hThread = CreateThread(NULL, NULL, (LPTHREAD_START_ROUTINE)&DialogBoxParamWrapper, &args, NULL, NULL);
	if (hThread)
	{
		WaitForSingleObject(hThread, 10);
		SendMessage(GetDlgItem(mywindow, IDC_PROGRESS1), PBM_SETRANGE, 0, MAKELPARAM(0, numSegments - 1));
	}

	for (int segnum = 0; segnum < numSegments; ++segnum)
	{
		if (segnum % 4 == 0) // ��������� ������� ��� � 4 �������� ����� ������ ���������
			SendMessage(GetDlgItem(mywindow, IDC_PROGRESS1), PBM_SETPOS, segnum, NULL);
		wstring wstr = L"In Progress... " + to_wstring(int(((myflo)segnum / numSegments) * 100)) + L" %";
		SetDlgItemText(mywindow, IDC_STATIC, wstr.c_str());

		vector <myflo> vY,
					   vres, 
					   vfilt,
					   vdiff,
					   vcoeffs = { S, linfitP, filtS , (myflo)fuel, (myflo)Num_iter };

		vY.assign(vSignal.begin() + vStartSegIndxs.at(segnum) + one_segment_width * leftP,
			vSignal.begin() + vStartSegIndxs.at(segnum) + one_segment_width * leftP + vSegPila.size());

		fdata.SetOriginSegment(vY, segnum);

		/*if (make_one_segment(1, vSegPila, vY, vres, vfilt, vdiff, vcoeffs) < 0)
			continue;*/
		make_one_segment(1, vSegPila, vY, vres, vfilt, vdiff, vcoeffs);

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