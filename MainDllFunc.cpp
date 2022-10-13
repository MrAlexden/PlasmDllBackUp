#include "MainDllFunc.h"

extern "C" __declspec(dllexport) int MainWrapper(_In_ int diagnostics,               // diagnostics type (zond::0|setka::1|cilind::2)
                                                 _In_ myflo * arrPila,               // 1д массив с данными пилы
                                                 _In_ myflo * arrSignal,             // 1д массив с данными сигнала
                                                 _In_ myflo * AdditionalData,        // 1д дополнительная информация об импульсе (!размер 14!)
                                                 _Out_ myflo * OriginalData,         // возвращаемый 1д массив с разделенными на отрезки оригинальными данными
                                                 _Out_ myflo * ResultData,           // возвращяемый 1д массив с результатом обработки
                                                 _Out_ myflo * FiltedData,           // возвращяемый 1д массив с отфильтрованными данными (если эта функция выбрана)
                                                 _Out_ myflo * Parameters)           // возвращаемый 1д массив с параметрами резултата обработки (температура, плотность и тд)
{
    if (arrPila == nullptr || 
        arrSignal == nullptr || 
        AdditionalData == nullptr ||
        OriginalData == nullptr ||
        ResultData == nullptr ||
        FiltedData == nullptr ||
        Parameters == nullptr)
    {
        //MessageBoxA(NULL, "Input data pointers equals nullptr", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_BadInputVecs;
    }

    if (AdditionalData[0] == 0 ||
        AdditionalData[7] == 0 ||
        AdditionalData[8] == 0 ||
        AdditionalData[9] == 0 ||
        AdditionalData[11] == 0 ||
        AdditionalData[12] == 0 ||
        AdditionalData[13] == 0)
    {
        //MessageBoxA(NULL, "Input data error, some values equals 0", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_ZeroInputVals;
    }

    if (diagnostics < 0 ||
        diagnostics > 2)
    {
        //MessageBoxA(NULL, "Diagnostics number < 0 or > 2", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_BadDiagNum;
    }

    int arrPsize = (int)AdditionalData[12];				// размер входного массива пилы
    int arrSsize = (int)AdditionalData[13];				// размер входного массива сигнала
    int d1, d2, d3, err = 0;

    vector <myflo> vPila(arrPila, arrPila + arrPsize);
    vector <myflo> vSignal(arrSignal, arrSignal + arrSsize);
    vector <myflo> vAdd(AdditionalData, AdditionalData + 14);

    Plasma_proc_result * fdata = new Plasma_proc_result;
    switch (diagnostics)
    {
    case 0: // Zond
        ERR(Zond(vPila, vSignal, vAdd, *fdata)); // вызов обработчика зонда
        break;
    case 1: // Setka
        ERR(Setka(vPila, vSignal, vAdd, *fdata)); // вызов обработчика сетки
        break;
    case 2: // Cilinder|Magnit
        ERR(Cilinder(vPila, vSignal, vAdd, *fdata)); // вызов обработчика цилиндра|магнита
        break;
    }

    d1 = fdata->Get_NumberOfSegments(), d2 = fdata->Get_SizeOfSegment(), d3 = fdata->Get_NumberOfParameters();
    try
    {
        for (int i = 0, j = 0, k = 0; i < d1; ++i, j += d2, k += d3)
        {
            memcpy(&OriginalData[j], &fdata->Get_mOriginalData()[i][0], sizeof myflo * d2);
            memcpy(&FiltedData[j], &fdata->Get_mFeltrationData()[i][0], sizeof myflo * d2);
            memcpy(&ResultData[j], &fdata->Get_mApproximatedData()[i][0], sizeof myflo * d2);
            memcpy(&Parameters[k], &fdata->Get_mParametersData()[i][0], sizeof myflo * d3);
        }
    }
    catch (...)
    {
        //MessageBoxA(NULL, "Index is out of range. Programm continued running", "Error!", MB_ICONWARNING | MB_OK);
    }

Error:
    return err;
}

extern "C" __declspec(dllexport) int FindSignalWrapper(_In_ myflo * arrPila,           // 1д массив с данными пилы
                                                       _In_ myflo * arrSignal,         // 1д массив с данными сигнала
                                                       _In_ myflo * AdditionalData,    // 1д дополнительная информация об импульсе (!размер 14!)
                                                       _Out_ int & DIM1,               // выходное значение (количество строк в матрице (количество отрезков))
                                                       _Out_ int & DIM2,               // выходное значение (количество столбиков в матрице (количество точек на отрезк))
                                                       _Out_ myflo * vResP)            // возвращяемый 1д массив с одним сегментом пилы (X для построения графика)
{
    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        AdditionalData == nullptr)
    {
        //MessageBoxA(NULL, "Input data pointers equals nullptr", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_BadInputVecs;
    }

    if (AdditionalData[0] == 0 ||
        AdditionalData[7] == 0 ||
        AdditionalData[8] == 0 ||
        AdditionalData[9] == 0 ||
        AdditionalData[11] == 0 ||
        AdditionalData[12] == 0 ||
        AdditionalData[13] == 0)
    {
        //MessageBoxA(NULL, "Input data error, some values equals 0", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_ZeroInputVals;
    }

    st_time_end_time[0] = AdditionalData[1];		    // время начала обработки (если этот параметр не выбран: -1)
    st_time_end_time[1] = AdditionalData[2];		    // время конца обработки (если этот параметр не выбран: -1)
    leftP = AdditionalData[3];					        // часть точек отсечки слева
    rightP = AdditionalData[4];					        // часть точек отсечки справа
    freqP = (int)AdditionalData[7];					    // частота пилы
    int resistance = (int)AdditionalData[8];			// сопротивление на зонде
    int coefPila = (int)AdditionalData[9];				// коэффициент усиления пилы
    int arrPsize = (int)AdditionalData[12];				// размер входного массива пилы
    int arrSsize = (int)AdditionalData[13];				// размер входного массива сигнала

    int numSegments = 0,
        err = 0;

    vector <myflo> vPila(arrPila, arrPila + arrPsize);
    vector <myflo> vSignal(arrSignal, arrSignal + arrSsize);
    vector <myflo> vSegPila;
    vector <int> vStartSegIndxs;

    /* домножаем пилу на коэффициент усиления */
    vectormult(vPila, coefPila);
    /* переворачиваем ток чтобы смотрел вверх(если нужно), и делим на сопротивление */
    vectordiv(vSignal, resistance);

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

    DIM1 = numSegments;
    DIM2 = vSegPila.size();
    memcpy(vResP, vSegPila.data(), sizeof myflo * vSegPila.size());

Error:
    if (err < 0)
        MessageBoxA(NULL, ERR_GetErrorDescription(err).c_str(), "Error!", MB_ICONWARNING | MB_OK);
    return err;
}

extern "C" __declspec(dllexport) int SetUpPila(_In_ myflo * arrPila,           // 1д массив с данными пилы
                                               _In_ myflo * arrSignal,         // 1д массив с данными сигнала
                                               _In_ int arrPsize,              // размер пилы
                                               _In_ int arrSsize,              // размер сигнала
                                               _Out_ myflo * vResP)            // возвращяемый 1д массив с пилой для построения графика
{
    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        vResP == nullptr)
    {
        //MessageBoxA(NULL, "Input data pointers equals nullptr", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_BadInputVecs;
    }

    if (arrPsize == 0 ||
        arrSsize == 0)
    {
        //MessageBoxA(NULL, "Input data error, some values equals 0", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_ZeroInputVals;
    }

    st_time_end_time[0] = -1;		    // время начала обработки (если этот параметр не выбран: -1)
    st_time_end_time[1] = -1;		    // время конца обработки (если этот параметр не выбран: -1)
    leftP = 0.0;					    // часть точек отсечки слева
    rightP = 0.0;					    // часть точек отсечки справа

    int err = 0;

    vector <myflo> vPila(arrPila, arrPila + arrPsize);
    vector <myflo> vSignal(arrSignal, arrSignal + arrSsize);
    vector <myflo> vSegPila;
    vector <int> vStartSegIndxs;
    vector <int> vnIndices;

    ERR(find_signal_and_make_pila(vPila, vSignal, vSegPila, vStartSegIndxs));

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

    /* находим все участки на пиле при помощи пиков */
    PeakFinder::findPeaks(vPila, vnIndices, false);
    if (vnIndices.size() <= 3) // 3 потому что дальше в функцию я передаю начало второго и он нужен целиком, тоесть от второго до третьего индекса
    {
        //MessageBoxA(NULL, "Less then 4 segments found", "Error!", MB_ICONWARNING | MB_OK);
        return ERR_TooFewSegs;
    }
    if (vPila.size() != vSignal.size())
        vectordiv(vnIndices, (myflo)vPila.size() / vSignal.size());

    try
    {
        for (int i = 0; i < vnIndices.size() - 1; ++i)
            memcpy(&vResP[vnIndices[i]], vSegPila.data(), sizeof myflo * vSegPila.size());
    }
    catch (...)
    {
        //MessageBoxA(NULL, "Index is out of range. Programm continued running", "Error!", MB_ICONWARNING | MB_OK);
    }

    try
    {
        if (vSegPila.size() >= vnIndices[0])
            memcpy(&vResP[0], &vSegPila[vSegPila.size() - vnIndices[0]], sizeof myflo * vnIndices[0]);
        else
        {
            memcpy(&vResP[vnIndices[0] - vSegPila.size()], vSegPila.data(), sizeof myflo * vSegPila.size());
            memcpy(&vResP[0], &vSegPila[vSegPila.size() - (vnIndices[0] - vSegPila.size())], sizeof myflo * (vnIndices[0] - vSegPila.size()));
        }
    }
    catch (...)
    {
        //MessageBoxA(NULL, "Index is out of range. Programm continued running", "Error!", MB_ICONWARNING | MB_OK);
    }

    try
    {
        if (vSegPila.size() >= (arrSsize - vnIndices.back()))
            memcpy(&vResP[vnIndices.back()], vSegPila.data(), sizeof myflo * (arrSsize - vnIndices.back()));
        else
        {
            memcpy(&vResP[vnIndices.back()], vSegPila.data(), sizeof myflo * vSegPila.size());
            memcpy(&vResP[vnIndices.back() + vSegPila.size()], vSegPila.data(), sizeof myflo * (arrSsize - vnIndices.back() + vSegPila.size()));
        }
    }
    catch (...)
    {
        //MessageBoxA(NULL, "Index is out of range. Programm continued running", "Error!", MB_ICONWARNING | MB_OK);
    }

    try
    {
        for (int i = 1; i < arrSsize - 1; ++i)
            if (vResP[i] == 0) 
                vResP[i] = vResP[i - 1];
    }
    catch (...)
    {
        //MessageBoxA(NULL, "Index is out of range. Programm continued running", "Error!", MB_ICONWARNING | MB_OK);
    }

Error:
    return err;
}