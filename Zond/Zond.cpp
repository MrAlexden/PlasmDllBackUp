#include "Zond.h"

int Zond(_In_ vector <myflo> vPila, 
	     _In_ vector <myflo> vSignal, 
	     _In_ vector <myflo> AdditionalData,
	     _Out_ Plasma_proc_result & fdata)
{
	if (vPila.size() == 0
		|| vPila.empty()
		|| vSignal.size() == 0
		|| vSignal.empty()
		|| AdditionalData.size() == 0
		|| AdditionalData.empty())
	{
		//MessageBoxA(NULL, "Corrupted input vectors", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_BadInputVecs;
	}
	if (AdditionalData[0] == 0 || is_invalid(AdditionalData[0]) ||
		(int)AdditionalData[7] == 0 || is_invalid(AdditionalData[7]) ||
		(int)AdditionalData[8] == 0 || is_invalid(AdditionalData[8]) ||
		(int)AdditionalData[9] == 0 || is_invalid(AdditionalData[9]) ||
		(int)AdditionalData[11] == 0 || is_invalid(AdditionalData[11]))
	{
		//MessageBoxA(NULL, "Input data error, some values equals 0", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_ZeroInputVals;
	}
	if (AdditionalData[3] >= 0.5 || AdditionalData[3] < 0.0)
	{
		//MessageBoxA(NULL, "�ut-off points on the left value must be > 0.0 and < 0.5", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_BadCutOffLeft;
	}
	if (AdditionalData[4] >= 0.5 || AdditionalData[4] < 0.0)
	{
		//MessageBoxA(NULL, "�ut-off points on the right value must be > 0.0 and < 0.5", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_BadCutOffRight;
	}

	vector <myflo> vSegPila;
	vector <int> vStartSegIndxs;

	int	numSegments = 0,	// ���������� �������� � ��������
		resistance = 0,
		coefPila = 0,
		dimension = 6,		// ���������� ��������� parameters
		err = 0;

	S = AdditionalData[0];						    // ������� ����������� �����
	st_time_end_time[0] = AdditionalData[1];		// ����� ������ ��������� (���� ���� �������� �� ������: -1)
	st_time_end_time[1] = AdditionalData[2];		// ����� ����� ��������� (���� ���� �������� �� ������: -1)
	leftP = AdditionalData[3];					    // ����� ����� ������� �����
	rightP = AdditionalData[4];					    // ����� ����� ������� ������
	linfitP = AdditionalData[5];					// ����� ����� ������� �������������
	filtS = AdditionalData[6];					    // ����� ����� ���������� �������
	freqP = (int)AdditionalData[7];					// ������� ����
	resistance = (int)AdditionalData[8];			// ������������� �� �����
	coefPila = (int)AdditionalData[9];				// ����������� �������� ����
	fuel = (int)AdditionalData[10];					// ������� �������� (He::0|Ar::1)
	Num_iter = (int)AdditionalData[11];				// ���������� �������� �������������(������ ������ �� �������� ������ ���������)

	/* ��������� ���� �� ����������� �������� */
	vectormult(vPila, coefPila);
	/* �������������� ��� ����� ������� �����(���� �����), � ����� �� ������������� */
	vectordiv(vSignal, -resistance);

	if (is_invalid(vPila[0])
		|| is_invalid(vPila[vPila.size() - 1])
		|| is_invalid(vSignal[0])
		|| is_invalid(vSignal[vSignal.size() - 1]))
	{
		//MessageBoxA(NULL, "Error after Pila|Signal factorizing", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_BadFactorizing;
	}

	ERR(find_signal_and_make_pila(vPila, vSignal, vSegPila, vStartSegIndxs));
	numSegments = vStartSegIndxs.size();

	if (vStartSegIndxs.size() == 0
		|| vStartSegIndxs.empty()
		|| is_invalid(vSegPila[0])
		|| is_invalid(vSegPila[vSegPila.size() - 1])
		|| is_invalid(vStartSegIndxs[0])
		|| is_invalid(vStartSegIndxs[vStartSegIndxs.size() - 1]))
	{
		//MessageBoxA(NULL, "Error after noise extracting", "Error!", MB_ICONWARNING | MB_OK);
		return ERR_BadNoise;
	}

	fdata.SetSegmentsNumber(numSegments);
	fdata.SetSegmentsSize(vSegPila.size());
	fdata.SetParamsNumber(dimension);

	fdata.SetPila(vSegPila);
	
#pragma omp parallel for schedule(static, 1) 
	for (int segnum = 0; segnum < numSegments; ++segnum)
	{
		vector <myflo> vY, vres, vfilt, vcoeffs = { filtS , linfitP };

		vY.assign(vSignal.begin() + vStartSegIndxs[segnum] + one_segment_width * leftP,
			vSignal.begin() + vStartSegIndxs[segnum] + one_segment_width * leftP + vSegPila.size());

		fdata.SetOriginSegment(vY, segnum);

		if (make_one_segment(0, vSegPila, vY, vres, vfilt, vcoeffs) < 0)
			continue;

		vcoeffs.insert(vcoeffs.begin(), vStartSegIndxs[segnum] * (1.0 / (one_segment_width * freqP)));

		fdata.SetFiltedSegment(vfilt, segnum);
		fdata.SetApproxSegment(vres, segnum);
		fdata.SetParamsSegment(vcoeffs, segnum);
	}

Error:
	return err;
}