#include "MainDllFunc.h"

BOOL WINAPI DllMain(HINSTANCE hinstDLL,  // handle to DLL module
                    DWORD fdwReason,     // reason for calling function
                    LPVOID lpvReserved)  // reserved
{
    // Perform actions based on the reason for calling.
    switch (fdwReason)
    {
    case DLL_PROCESS_ATTACH:
        // Initialize once for each new process.
        // Return FALSE to fail DLL load.
        hInstThisDll = hinstDLL;
        break;

    case DLL_THREAD_ATTACH:
        // Do thread-specific initialization.
        break;

    case DLL_THREAD_DETACH:
        // Do thread-specific cleanup.
        break;

    case DLL_PROCESS_DETACH:

        if (lpvReserved != nullptr)
        {
            break; // do not do cleanup if process termination scenario
        }

        // Perform any necessary cleanup.
        break;
    }
    return TRUE;  // Successful DLL_PROCESS_ATTACH.
}

extern "C" __declspec(dllexport) int MainWrapper(_In_ int diagnostics,               // diagnostics type (zond::0|setka::1|cilind::2)
                                                 _In_ myflo * arrPila,               // 1д массив с данными пилы
                                                 _In_ myflo * arrSignal,             // 1д массив с данными сигнала
                                                 _In_ myflo * AdditionalData,        // 1д дополнительная информация об импульсе (!размер 14!)
                                                 _Out_ myflo * OriginalData,         // возвращаемый 1д массив с разделенными на отрезки оригинальными данными
                                                 _Out_ myflo * ResultData,           // возвращяемый 1д массив с результатом обработки
                                                 _Out_ myflo * FiltedData,           // возвращяемый 1д массив с отфильтрованными данными (если эта функция выбрана)
                                                 _Out_ myflo * DiffedData,			 // возвращяемый 1д массив с данными производной (если эта функция выбрана)
                                                 _Out_ myflo * Parameters)           // возвращаемый 1д массив с параметрами резултата обработки (температура, плотность и тд)
{
    int err = 0;

    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        AdditionalData == nullptr ||
        OriginalData == nullptr ||
        ResultData == nullptr ||
        FiltedData == nullptr ||
        DiffedData == nullptr ||
        Parameters == nullptr)
        ERR(ERR_BadInputVecs);
    if (AdditionalData[0] == 0 ||
        AdditionalData[7] == 0 ||
        AdditionalData[8] == 0 ||
        AdditionalData[9] == 0 ||
        AdditionalData[11] == 0 ||
        AdditionalData[12] <= 1 ||
        AdditionalData[13] <= 1)
        ERR(ERR_ZeroInputVals);
    if (diagnostics < 0 ||
        diagnostics > 2)
        ERR(ERR_BadDiagNum);
    if (AdditionalData[3] + AdditionalData[4] > 0.9
        || AdditionalData[3] < 0.0
        || AdditionalData[4] < 0.0)
        ERR(ERR_BadCutOff);
    if (AdditionalData[5] < 0 || AdditionalData[5] > 0.9)
        ERR(ERR_BadLinFit);

    {/*********************************** засовываю в блок чтобы goto не ругался **********************************/
        int arrPsize = (int)AdditionalData[12];				// размер входного массива пилы
        int arrSsize = (int)AdditionalData[13];				// размер входного массива сигнала
        int d1, d2, d3;

        vector <myflo> vPila(arrPila, arrPila + arrPsize);
        vector <myflo> vSignal(arrSignal, arrSignal + arrSsize);
        vector <myflo> vAdd(AdditionalData, AdditionalData + 14);

        Plasma_proc_result <myflo> * fdata = new Plasma_proc_result <myflo>;

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
        default:
            __assume(0);
        }

        d1 = fdata->Get_NumberOfSegments(), d2 = fdata->Get_SizeOfSegment(), d3 = fdata->Get_NumberOfParameters();
        try
        {
            for (int i = 0, j = 0, k = 0; i < d1; ++i, j += d2, k += d3)
            {
                memcpy(&OriginalData[j], &fdata->Get_mOriginalData().at(i).at(0), sizeof myflo * d2);
                memcpy(&FiltedData[j], &fdata->Get_mFiltratedData().at(i).at(0), sizeof myflo * d2);
                memcpy(&DiffedData[j], &fdata->Get_mDifferentiatedData().at(i).at(0), sizeof myflo * d2);
                memcpy(&ResultData[j], &fdata->Get_mApproximatedData().at(i).at(0), sizeof myflo * d2);
                memcpy(&Parameters[k], &fdata->Get_mParametersData().at(i).at(0), sizeof myflo * d3);
            }
        }
        catch (...)
        {
            err = ERR_Exception;
        }

        delete fdata;
    }/*************************************************************************************************************/

Error:
    return err;
}

extern "C" __declspec(dllexport) int FindSignalWrapper(_In_ myflo * arrPila,                 // 1д массив с данными пилы
                                                       _In_ myflo * arrSignal,               // 1д массив с данными сигнала
                                                       _In_ myflo * AdditionalData,          // 1д дополнительная информация об импульсе (!размер 14!)
                                                       _Out_ int & DIM1,                     // выходное значение (количество строк в матрице (количество отрезков))
                                                       _Out_ int & DIM2,                     // выходное значение (количество столбиков в матрице (количество точек на отрезк))
                                                       _Out_ myflo * vResP)                  // возвращяемый 1д массив с одним сегментом пилы (X для построения графика)
{
    int err = 0;

    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        AdditionalData == nullptr)
        ERR(ERR_BadInputVecs);
    if (AdditionalData[0] == 0 ||
        AdditionalData[7] == 0 ||
        AdditionalData[8] == 0 ||
        AdditionalData[9] == 0 ||
        AdditionalData[11] == 0 ||
        AdditionalData[12] <= 1 ||
        AdditionalData[13] <= 1)
        ERR(ERR_ZeroInputVals);
    if (AdditionalData[3] + AdditionalData[4] > 0.9
        || AdditionalData[3] < 0.0
        || AdditionalData[4] < 0.0)
        ERR(ERR_BadCutOff);
    if (AdditionalData[5] < 0 || AdditionalData[5] > 0.9)
        ERR(ERR_BadLinFit);
    
    {/*********************************** засовываю в блок чтобы goto не ругался **********************************/
        st_time_end_time[0] = AdditionalData[1];		    // время начала обработки (если этот параметр не выбран: -1)
        st_time_end_time[1] = AdditionalData[2];		    // время конца обработки (если этот параметр не выбран: -1)
        leftP = AdditionalData[3];					        // часть точек отсечки слева
        rightP = AdditionalData[4];					        // часть точек отсечки справа
        freqP = (int)AdditionalData[7];					    // частота пилы
        int resistance = (int)AdditionalData[8];			// сопротивление на зонде
        int coefPila = (int)AdditionalData[9];				// коэффициент усиления пилы
        int arrPsize = (int)AdditionalData[12];				// размер входного массива пилы
        int arrSsize = (int)AdditionalData[13];				// размер входного массива сигнала

        int numSegments = 0;

        vector <myflo> vPila(arrPila, arrPila + arrPsize);
        vector <myflo> vSignal(arrSignal, arrSignal + arrSsize);
        vector <myflo> vSegPila;
        vector <int> vStartSegIndxs;

        /* домножаем пилу на коэффициент усиления */
        vectormult(vPila, coefPila);
        /* переворачиваем ток чтобы смотрел вверх(если нужно), и делим на сопротивление */
        vectordiv(vSignal, resistance);

        if (is_invalid(vPila[0])
            || is_invalid(vPila.back())
            || is_invalid(vSignal[0])
            || is_invalid(vSignal.back()))
            ERR(ERR_BadFactorizing);

        ERR(find_signal_and_make_pila(vPila, vSignal, vSegPila, vStartSegIndxs));
        numSegments = vStartSegIndxs.size();
        
        if (vStartSegIndxs.size() == 0
            || vStartSegIndxs.empty()
            || is_invalid(vSegPila[0])
            || is_invalid(vSegPila.back())
            || is_invalid(vStartSegIndxs[0])
            || is_invalid(vStartSegIndxs.back()))
            ERR(ERR_BadNoise);

        DIM1 = numSegments;
        DIM2 = vSegPila.size();
        memcpy(vResP, vSegPila.data(), sizeof myflo * vSegPila.size());
    }/*************************************************************************************************************/
    
Error:
    if (err < 0)
        MessageBoxA(NULL, ERR_GetErrorDescription(err).c_str(), "Error!", MB_ICONINFORMATION | MB_OK);
    return err;
}

extern "C" __declspec(dllexport) int SetUpPila(_In_ int diagnostics,           // diagnostics type (zond::0|setka::1|cilind::2)
                                               _In_ myflo * arrPila,           // 1д массив с данными пилы
                                               _In_ myflo * arrSignal,         // 1д массив с данными сигнала
                                               _In_ int arrPsize,              // размер пилы
                                               _In_ int arrSsize,              // размер сигнала
                                               _Out_ myflo * vResP)            // возвращяемый 1д массив с пилой для построения графика
{
    int err = 0;

    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        vResP == nullptr)
        ERR(ERR_BadInputVecs);
    if (arrPsize <= 1 ||
        arrSsize <= 1)
        ERR(ERR_ZeroInputVals);

    {/*********************************** засовываю в блок чтобы goto не ругался **********************************/
        st_time_end_time[0] = -1;		    // время начала обработки (если этот параметр не выбран: -1)
        st_time_end_time[1] = -1;		    // время конца обработки (если этот параметр не выбран: -1)
        leftP = 0.0;					    // часть точек отсечки слева
        rightP = 0.0;					    // часть точек отсечки справа

        const vector <myflo> vPila(arrPila, arrPila + arrPsize);
        vector <myflo> vSignal(arrSignal, arrSignal + arrSsize);
        vector <myflo> vSegPila, IndHandler;
        vector <int> vStartSegIndxs;
        vector <int> vnIndices;
        myflo vless = 0.0, vmore = 0.0;
        
        ERR(find_signal_and_make_pila(vPila, vSignal, vSegPila, vStartSegIndxs));

        if (vStartSegIndxs.size() == 0
            || vStartSegIndxs.empty()
            || is_invalid(vSegPila[0])
            || is_invalid(vSegPila.back())
            || is_invalid(vStartSegIndxs[0])
            || is_invalid(vStartSegIndxs.back()))
            ERR(ERR_BadNoise);

        /* находим все участки на пиле при помощи пиков */
        err = PeakFinder::findPeaks(vPila, vnIndices, false);
        if (vnIndices.size() <= 3)
            ERR(ERR_TooFewSegs);

        /* приводим индексы начал отрезков к одной фазе с током */
        match_pila_and_signal(vPila, vSignal, vnIndices);

#ifdef KLUDGE
        /* ВРЕМЕННЫЙ КОСТЫЛЬ */
        if (diagnostics == 0 || diagnostics == 1)
        {
            IndHandler.assign(vnIndices.begin(), vnIndices.end());
            vectormult(IndHandler, 1.001f);
            vnIndices.assign(IndHandler.begin(), IndHandler.end());
            /*проверяем превышение крайнего индекса*/
            while (vnIndices.back() > vSignal.size())
                vnIndices.pop_back();
        }
        /* ВРЕМЕННЫЙ КОСТЫЛЬ */
#endif

        try
        {
            /* заполняем основную часть массива с пилой */
            for (int i = 0; i < vnIndices.size() - 1; ++i)
                if(vnIndices.at(i) + vSegPila.size() < arrSsize)
                    memcpy(&vResP[vnIndices.at(i)], vSegPila.data(), sizeof myflo * vSegPila.size());
        }
        catch (...)
        {
            err = ERR_Exception;
        }

        try
        {
            /* заполняем начало массива с пилой */
            if (vSegPila.size() >= vnIndices.at(0))
                memcpy(&vResP[0], &vSegPila.at(vSegPila.size() - vnIndices.at(0)), sizeof myflo * vnIndices.at(0));
            else
            {
                memcpy(&vResP[vnIndices.at(0) - vSegPila.size()], vSegPila.data(), sizeof myflo * vSegPila.size());
                memcpy(&vResP[0], &vSegPila.at(vSegPila.size() - (vnIndices.at(0) - vSegPila.size())), sizeof myflo * (vnIndices.at(0) - vSegPila.size()));
            }
        }
        catch (...)
        {
            err = ERR_Exception;
        }

        /* заполняем конец массива с пилой */
        if (vSegPila.size() >= (arrSsize - vnIndices.back()))
            memcpy(&vResP[vnIndices.back()], vSegPila.data(), sizeof myflo * (arrSsize - vnIndices.back()));
        else
        {
            memcpy(&vResP[vnIndices.back()], vSegPila.data(), sizeof myflo * vSegPila.size());
            memcpy(&vResP[vnIndices.back() + vSegPila.size()], vSegPila.data(), sizeof myflo * (arrSsize - vnIndices.back() - vSegPila.size()));
        }


        /* убираю нули если они есть на всем протяжении массива */
        if (vSegPila.front() < vSegPila.back()) // для случая '/'
        {
            vless = vSegPila.front();
            vmore = vSegPila.back();
        }
        else                                    // для случая '\'
        {
            vless = vSegPila.back();
            vmore = vSegPila.front();
        }
        for (int i = 1; i < arrSsize; ++i)
            if (vResP[i] == 0 || vResP[i] < vless || vResP[i] > vmore)
                vResP[i] = vResP[i - 1];
    }/*************************************************************************************************************/

Error:
    return err;
}

extern "C" __declspec(dllexport) int OriginAll(_In_ int diagnostics,               // diagnostics type (zond::0|setka::1|cilind::2)
                                               _In_ double* arrPila,               // 1д массив с данными пилы
                                               _In_ double* arrSignal,             // 1д массив с данными сигнала
                                               _In_ double* AdditionalData,        // 1д дополнительная информация об импульсе (!размер 14!)
                                               _Out_ double* OriginalData,         // возвращаемый 1д массив с разделенными на отрезки оригинальными данными
                                               _Out_ double* ResultData,           // возвращяемый 1д массив с результатом обработки
                                               _Out_ double* FiltedData,           // возвращяемый 1д массив с отфильтрованными данными (если эта функция выбрана)
                                               _Out_ double* DiffedData,		   // возвращяемый 1д массив с данными производной (если эта функция выбрана)
                                               _Out_ double* Parameters)           // возвращаемый 1д массив с параметрами резултата обработки (температура, плотность и тд))
{
    int err = 0;

    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        AdditionalData == nullptr ||
        OriginalData == nullptr ||
        ResultData == nullptr ||
        FiltedData == nullptr ||
        DiffedData == nullptr ||
        Parameters == nullptr)
        ERR(ERR_BadInputVecs);
    if (AdditionalData[0] == 0 ||
        AdditionalData[7] == 0 ||
        AdditionalData[8] == 0 ||
        AdditionalData[9] == 0 ||
        AdditionalData[11] == 0 ||
        AdditionalData[12] <= 1 ||
        AdditionalData[13] <= 1)
        ERR(ERR_ZeroInputVals);
    if (diagnostics < 0 ||
        diagnostics > 2)
        ERR(ERR_BadDiagNum);
    if (AdditionalData[3] + AdditionalData[4] > 0.9
        || AdditionalData[3] < 0.0
        || AdditionalData[4] < 0.0)
        ERR(ERR_BadCutOff);
    if (AdditionalData[5] < 0 || AdditionalData[5] > 0.9)
        ERR(ERR_BadLinFit);

    {/*********************************** засовываю в блок чтобы goto не ругался **********************************/
        int arrPsize = (int)AdditionalData[12];				// размер входного массива пилы
        int arrSsize = (int)AdditionalData[13];				// размер входного массива сигнала
        int d1, d2, d3;

        vector <double> vPila(arrPila, arrPila + arrPsize);
        vector <double> vSignal(arrSignal, arrSignal + arrSsize);
        vector <double> vAdd(AdditionalData, AdditionalData + 14);

        vector <myflo> vP(vPila.begin(), vPila.end());
        vector <myflo> vS(vSignal.begin(), vSignal.end());
        vector <myflo> vA(vAdd.begin(), vAdd.end());

        Plasma_proc_result <myflo> * fdata = new Plasma_proc_result <myflo>;

        switch (diagnostics)
        {
        case 0: // Zond
            ERR(Zond(vP, vS, vA, *fdata)); // вызов обработчика зонда
            break;
        case 1: // Setka
            ERR(Setka(vP, vS, vA, *fdata)); // вызов обработчика сетки
            break;
        case 2: // Cilinder|Magnit
            ERR(Cilinder(vP, vS, vA, *fdata)); // вызов обработчика цилиндра|магнита
            break;
        default:
            __assume(0);
        }

        d1 = fdata->Get_NumberOfSegments(), d2 = fdata->Get_SizeOfSegment(), d3 = fdata->Get_NumberOfParameters();
        try
        {
            vector <double> v;
            for (int i = 0, j = 0, k = 0; i < d1; ++i, j += d2, k += d3)
            {
                v.assign(fdata->Get_mOriginalData().at(i).begin(), fdata->Get_mOriginalData().at(i).end());
                memcpy(&OriginalData[j], v.data(), sizeof(double) * d2);
                
                v.assign(fdata->Get_mApproximatedData().at(i).begin(), fdata->Get_mApproximatedData().at(i).end());
                memcpy(&ResultData[j], v.data(), sizeof (double) * d2);
                
                v.assign(fdata->Get_mFiltratedData().at(i).begin(), fdata->Get_mFiltratedData().at(i).end());
                memcpy(&FiltedData[j], v.data(), sizeof (double) * d2);

                v.assign(fdata->Get_mDifferentiatedData().at(i).begin(), fdata->Get_mDifferentiatedData().at(i).end());
                memcpy(&DiffedData[j], v.data(), sizeof (double) * d2);
                
                v.assign(fdata->Get_mParametersData().at(i).begin(), fdata->Get_mParametersData().at(i).end());
                memcpy(&Parameters[k], v.data(), sizeof (double) * d3);
            }
        }
        catch (...)
        {
            err = ERR_Exception;
        }

        delete fdata;
    }/*************************************************************************************************************/

Error:
    return err;
}

extern "C" __declspec(dllexport) int OriginFindSignal(_In_ int diagnostics,                 // diagnostics type (zond::0|setka::1|cilind::2)
                                                      _In_ double* arrPila,                 // 1д массив с данными пилы
                                                      _In_ double* arrSignal,               // 1д массив с данными сигнала
                                                      _In_ double* AdditionalData,          // 1д дополнительная информация об импульсе (!размер 14!)
                                                      _Out_ int& DIM1,                      // выходное значение (количество строк в матрице (количество отрезков))
                                                      _Out_ int& DIM2,                      // выходное значение (количество столбиков в матрице (количество точек на отрезк))
                                                      _Out_ double* vResP,                  // возвращяемый 1д массив с одним сегментом пилы (X для построения графика)
                                                      _Out_ int* vSsI)                      // возвращаемы 1д массив с индексами начал сегментов
{
    int err = 0;

    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        AdditionalData == nullptr)
        ERR(ERR_BadInputVecs);
    if (AdditionalData[0] == 0 ||
        AdditionalData[7] == 0 ||
        AdditionalData[8] == 0 ||
        AdditionalData[9] == 0 ||
        AdditionalData[11] == 0 ||
        AdditionalData[12] <= 1 ||
        AdditionalData[13] <= 1)
        ERR(ERR_ZeroInputVals);
    if (AdditionalData[3] + AdditionalData[4] > 0.9
        || AdditionalData[3] < 0.0
        || AdditionalData[4] < 0.0)
        ERR(ERR_BadCutOff);
    if (AdditionalData[5] < 0 || AdditionalData[5] > 0.9)
        ERR(ERR_BadLinFit);

    {/*********************************** засовываю в блок чтобы goto не ругался **********************************/
        st_time_end_time[0] = AdditionalData[1];		    // время начала обработки (если этот параметр не выбран: -1)
        st_time_end_time[1] = AdditionalData[2];		    // время конца обработки (если этот параметр не выбран: -1)
        leftP = AdditionalData[3];					        // часть точек отсечки слева
        rightP = AdditionalData[4];					        // часть точек отсечки справа
        freqP = (int)AdditionalData[7];					    // частота пилы
        int resistance = (int)AdditionalData[8];			// сопротивление на зонде
        int coefPila = (int)AdditionalData[9];				// коэффициент усиления пилы
        int arrPsize = (int)AdditionalData[12];				// размер входного массива пилы
        int arrSsize = (int)AdditionalData[13];				// размер входного массива сигнала

        int numSegments = 0;

        vector <double> vPila(arrPila, arrPila + arrPsize);
        vector <double> vSignal(arrSignal, arrSignal + arrSsize);
        vector <double> vSegPila;
        vector <int> vStartSegIndxs;

        /* домножаем пилу на коэффициент усиления */
        vectormult(vPila, coefPila);
        /* переворачиваем ток чтобы смотрел вверх(если нужно), и делим на сопротивление */
        vectordiv(vSignal, resistance);

        if (is_invalid(vPila[0])
            || is_invalid(vPila.back())
            || is_invalid(vSignal[0])
            || is_invalid(vSignal.back()))
            ERR(ERR_BadFactorizing);
        
        vector <myflo> vP(vPila.begin(), vPila.end());
        vector <myflo> vS(vSignal.begin(), vSignal.end());
        vector <myflo> vSP;

        ERR(find_signal_and_make_pila(vP, vS, vSP, vStartSegIndxs));
        numSegments = vStartSegIndxs.size();

        vPila.assign(vP.begin(), vP.end());
        vSignal.assign(vS.begin(), vS.end());
        vSegPila.assign(vSP.begin(), vSP.end());
        
        if (vStartSegIndxs.size() == 0
            || vStartSegIndxs.empty()
            || is_invalid(vSegPila[0])
            || is_invalid(vSegPila.back())
            || is_invalid(vStartSegIndxs[0])
            || is_invalid(vStartSegIndxs.back()))
            ERR(ERR_BadNoise);

#ifdef KLUDGE
        /* ВРЕМЕННЫЙ КОСТЫЛЬ */
        if (diagnostics == 0 || diagnostics == 1)
        {
            vector <double> IndHandler(vStartSegIndxs.begin(), vStartSegIndxs.end());
            vectormult(IndHandler, 1.001f);
            vStartSegIndxs.assign(IndHandler.begin(), IndHandler.end());
            /*проверяем превышение крайнего индекса*/
            while (vStartSegIndxs.back() > vSignal.size())
                vStartSegIndxs.pop_back();
        }
        /* ВРЕМЕННЫЙ КОСТЫЛЬ */
#endif
        
        DIM1 = numSegments;
        DIM2 = vSegPila.size();
        memcpy(vResP, vSegPila.data(), sizeof (double) * vSegPila.size());
        memcpy(vSsI, vStartSegIndxs.data(), sizeof (int) * vStartSegIndxs.size());
    }/*************************************************************************************************************/

Error:
    return err;
}

extern "C" __declspec(dllexport) int OriginMakeOne(_In_ int diagnostics,                    // diagnostics type (zond::0|setka::1|cilind::2)
                                                   _In_ double* arrPila,                    // 1д массив с данными пилы
                                                   _In_ double* arrSignal,                  // 1д массив с данными сигнала
                                                   _In_ double* AdditionalData,             // 1д дополнительная информация об импульсе (!размер 6!)
                                                   _Out_ double* arrres,			        // array to be filled with the result
                                                   _Out_ double* arrfilt,			        // array to be filled with the filtration
                                                   _Out_ double* arrdiff,			        // array to be filled with the differentiation
                                                   _Out_ double* arrcoeffs)		            // additional coeffs/results vector
{
    int err = 0;

    if (arrPila == nullptr ||
        arrSignal == nullptr ||
        AdditionalData == nullptr)
        ERR(ERR_BadInputVecs);
    if (AdditionalData[0] == 0 ||
        AdditionalData[4] == 0 ||
        AdditionalData[5] <= 1)
        ERR(ERR_ZeroInputVals);
    if (diagnostics < 0 ||
        diagnostics > 2)
        ERR(ERR_BadDiagNum);
    if (AdditionalData[1] < 0 || AdditionalData[1] > 0.9)
        ERR(ERR_BadLinFit);
    
    {/*********************************** засовываю в блок чтобы goto не ругался **********************************/
        S = AdditionalData[0];
        linfitP = AdditionalData[1];
        filtS = AdditionalData[2];
        fuel = (int)AdditionalData[3];
        Num_iter = (int)AdditionalData[4];
        int arrsize = (int)AdditionalData[5];

        vector <double> vPila(arrPila, arrPila + arrsize);
        vector <double> vSignal(arrSignal, arrSignal + arrsize);
        vector <double> vres, vfilt, vdiff, vcoeffs;

        vector <myflo> vP(vPila.begin(), vPila.end());
        vector <myflo> vS(vSignal.begin(), vSignal.end());
        vector <myflo> vR, vF, vD, vC = { S, linfitP, filtS , (myflo)fuel, (myflo)Num_iter };

        ERR(make_one_segment(diagnostics, vP, vS, vR, vF, vD, vC));

        vres.assign(vR.begin(), vR.end());
        vfilt.assign(vF.begin(), vF.end());
        vdiff.assign(vD.begin(), vD.end());
        vcoeffs.assign(vC.begin(), vC.end());

        memcpy(arrres, vres.data(), sizeof (double) * vres.size());
        memcpy(arrfilt, vfilt.data(), sizeof (double) * vfilt.size());
        memcpy(arrdiff, vdiff.data(), sizeof (double) * vdiff.size());
        memcpy(arrcoeffs, vcoeffs.data(), sizeof (double) * vcoeffs.size());
    }/*************************************************************************************************************/

Error:
    return err;
}