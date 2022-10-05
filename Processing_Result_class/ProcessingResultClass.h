#pragma once

#ifndef PROC_RES
#define PROC_RES

#include <iostream>     // for std::
#include <vector>       // for vector

using namespace std;

typedef double myflo;

class Plasma_proc_result
{
public:

	Plasma_proc_result():
	NumberOfSegments(0),
	SizeOfSegment(0),
	NumberOfParameters(0)
	{
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

	void SetSegmentsNumber(int Number_Of_Parameters)
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

	vector < myflo > & Get_vPila()
	{
		return vPila;
	}

	vector < vector <myflo> > & Get_mOriginalData()
	{
		return mOriginalData;
	}

	vector < vector <myflo> > & Get_mFeltrationData()
	{
		return mFeltrationData;
	}

	vector < vector <myflo> > & Get_mApproximatedData()
	{
		return mApproximatedData;
	}

	vector < vector <myflo> > & Get_mParametersData()
	{
		return mParametersData;
	}

	int Get_NumberOfSegments()
	{
		return NumberOfSegments;
	}

	int Get_SizeOfSegment()
	{
		return SizeOfSegment;
	}

	int Get_NumberOfParameters()
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
