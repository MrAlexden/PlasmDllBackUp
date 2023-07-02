#pragma once

#ifndef PROC_RES
#define PROC_RES

#include <vector>       // for vector

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
		vRamp.clear();

		for (int i = 0; i < NumberOfSegments; ++i)
		{
			mOriginalData[i].clear();
			mFeltrationData[i].clear();
			mApproximatedData[i].clear();
			mDifferentiatedData[i].clear();
			mParametersData[i].clear();
		}

		mOriginalData.clear();
		mFeltrationData.clear();
		mApproximatedData.clear();
		mDifferentiatedData.clear();
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
		mDifferentiatedData.resize(Number_Of_Segments);
		mParametersData.resize(Number_Of_Segments);
	}
	void SetSegmentsNumber(int Number_Of_Segments)
	{
		NumberOfSegments = Number_Of_Segments;

		mOriginalData.resize(Number_Of_Segments);
		mFeltrationData.resize(Number_Of_Segments);
		mApproximatedData.resize(Number_Of_Segments);
		mDifferentiatedData.resize(Number_Of_Segments);
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


	/* DATA SETTING */
	void SetRamp(std::vector <T> & v)
	{
		vRamp = v;
	}
	int SetOriginSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mOriginalData[i] = v;
	}
	int SetFiltedSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mFeltrationData[i] = v;
	}
	int SetApproxSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mApproximatedData[i] = v;
	}
	int SetDiffedSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mDifferentiatedData[i] = v;
	}
	int SetParamsSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mParametersData[i] = v;
	}


	/* DATA GETTING */
	std::vector < T > & const Get_vRamp()
	{
		return vRamp;
	}
	std::vector < std::vector <T> > & const Get_mOriginalData()
	{
		return mOriginalData;
	}
	std::vector < std::vector <T> > & const Get_mFiltratedData()
	{
		return mFeltrationData;
	}
	std::vector < std::vector <T> > & const Get_mApproximatedData()
	{
		return mApproximatedData;
	}
	std::vector < std::vector <T> > & const Get_mDifferentiatedData()
	{
		return mDifferentiatedData;
	}
	std::vector < std::vector <T> > & const Get_mParametersData()
	{
		return mParametersData;
	}


	/* DIMENSION GETTING */
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

	std::vector <T> vRamp;

	std::vector < std::vector <T> > mOriginalData;
	std::vector < std::vector <T> > mFeltrationData;
	std::vector < std::vector <T> > mApproximatedData;
	std::vector < std::vector <T> > mDifferentiatedData;
	std::vector < std::vector <T> > mParametersData;

	int NumberOfSegments;
	int SizeOfSegment;
	int NumberOfParameters;
};
#endif
