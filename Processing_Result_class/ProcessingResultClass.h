#pragma once

#ifndef PROC_RES
#define PROC_RES

#include <iostream>     // for std::
#include <vector>       // for vector

using namespace std;

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
typedef float myflo;
#endif

template <typename T>
class Plasma_proc_result
{
public:

	Plasma_proc_result():
	NumberOfSegments(0),
	SizeOfSegment(0),
	NumberOfParameters(0)
	{
	}

	~Plasma_proc_result()
	{
		vPila.clear();

		for (int i = 0; i < NumberOfSegments; ++i)
		{
			mOriginalData[i].clear();
			mFeltrationData[i].clear();
			mApproximatedData[i].clear();
			mParametersData[i].clear();
		}

		mOriginalData.clear();
		mFeltrationData.clear();
		mApproximatedData.clear();
		mParametersData.clear();
	}

	Plasma_proc_result(int Number_Of_Segments, int Size_Of_Segment, int Number_Of_Parameters)
	{
		NumberOfSegments = Number_Of_Segments;
		SizeOfSegment = Size_Of_Segment;
		NumberOfParameters = Number_Of_Parameters;

		mOriginalData.resize(Number_Of_Segments);
		mFeltrationData.resize(Number_Of_Segments);
		mApproximatedData.resize(Number_Of_Segments);
		mParametersData.resize(Number_Of_Segments);
	}

	void SetSegmentsNumber(int Number_Of_Segments)
	{
		NumberOfSegments = Number_Of_Segments;

		mOriginalData.resize(Number_Of_Segments);
		mFeltrationData.resize(Number_Of_Segments);
		mApproximatedData.resize(Number_Of_Segments);
		mParametersData.resize(Number_Of_Segments);
	}

	void SetSegmentsSize(int Size_Of_Segment)
	{
		SizeOfSegment = Size_Of_Segment;
	}

	void SetParamsNumber(int Number_Of_Parameters)
	{
		NumberOfParameters = Number_Of_Parameters;
	}

	void SetPila(vector <T> & v)
	{
		vPila = v;
	}

	int SetOriginSegment(vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mOriginalData[i] = v;
	}

	int SetFiltedSegment(vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mFeltrationData[i] = v;
	}

	int SetApproxSegment(vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mApproximatedData[i] = v;
	}

	int SetParamsSegment(vector <T>& v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mParametersData[i] = v;
	}

	vector < T > & const Get_vPila()
	{
		return vPila;
	}

	vector < vector <T> > & const Get_mOriginalData()
	{
		return mOriginalData;
	}

	vector < vector <T> > & const Get_mFeltrationData()
	{
		return mFeltrationData;
	}

	vector < vector <T> > & const Get_mApproximatedData()
	{
		return mApproximatedData;
	}

	vector < vector <T> > & const Get_mParametersData()
	{
		return mParametersData;
	}

	int Get_NumberOfSegments() const
	{
		return NumberOfSegments;
	}

	int Get_SizeOfSegment() const
	{
		return SizeOfSegment;
	}

	int Get_NumberOfParameters() const
	{
		return NumberOfParameters;
	}

private:

	vector <T> vPila;

	vector < vector <T> > mOriginalData;
	vector < vector <T> > mFeltrationData;
	vector < vector <T> > mApproximatedData;
	vector < vector <T> > mParametersData;

	int NumberOfSegments;
	int SizeOfSegment;
	int NumberOfParameters;
};
#endif
