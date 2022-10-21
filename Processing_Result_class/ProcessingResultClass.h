#pragma once

#ifndef PROC_RES
#define PROC_RES

#include <iostream>     // for std::
#include <vector>       // for vector

using namespace std;

#ifndef WIN32_LEAN_AND_MEAN
typedef float myflo;
#endif

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

	void SetPila(vector <myflo> & v)
	{
		vPila = v;
	}

	int SetOriginSegment(vector <myflo> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mOriginalData[i] = v;
	}

	int SetFiltedSegment(vector <myflo> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mFeltrationData[i] = v;
	}

	int SetApproxSegment(vector <myflo> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mApproximatedData[i] = v;
	}

	int SetParamsSegment(vector <myflo>& v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mParametersData[i] = v;
	}

	vector < myflo > & const Get_vPila()
	{
		return vPila;
	}

	vector < vector <myflo> > & const Get_mOriginalData()
	{
		return mOriginalData;
	}

	vector < vector <myflo> > & const Get_mFeltrationData()
	{
		return mFeltrationData;
	}

	vector < vector <myflo> > & const Get_mApproximatedData()
	{
		return mApproximatedData;
	}

	vector < vector <myflo> > & const Get_mParametersData()
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

	vector <myflo> vPila;

	vector < vector <myflo> > mOriginalData;
	vector < vector <myflo> > mFeltrationData;
	vector < vector <myflo> > mApproximatedData;
	vector < vector <myflo> > mParametersData;

	int NumberOfSegments;
	int SizeOfSegment;
	int NumberOfParameters;
};
#endif
