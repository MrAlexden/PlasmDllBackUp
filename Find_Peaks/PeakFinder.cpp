#include "PeakFinder.h"

void diff(vector <myflo>& in, vector <myflo>& out, int k)
{
	out.resize(in.size());

	for (int i = 1; i < in.size(); ++i)
		out[i - 1] = (in[i] - in[i - 1]) * k;
	out[out.size() - 1] = out[out.size() - 2];
}

void vectorElementsProduct(vector <myflo>& a, vector <myflo>& b, vector <myflo>& out)
{
	out.resize(a.size());

	for (int i = 0; i < a.size(); ++i)
		out[i] = a[i] * b[i];
}

void findIndicesLessThan(vector <myflo>& in, myflo threshold, vector <int>& indices)
{
	for (int i = 0; i < in.size(); ++i)
		if (in[i] < threshold)
			indices.push_back(i + 1);
}

void selectElementsFromIndices(vector <myflo>& in, vector <int>& indices, vector <myflo>& out)
{
	out.resize(indices.size());

	for (int i = 0; i < indices.size(); ++i)
		out[i] = in[indices[i]];
}

void selectElementsFromIndices(vector <int>& in, vector <int>& indices, vector <int>& out)
{
	out.resize(indices.size());

	for (int i = 0; i < indices.size(); ++i)
		out[i] = in[indices[i]];
}

void signVector(vector <myflo>& in, vector<int>& out)
{
	out.resize(in.size());

	for (int i = 0; i < in.size(); ++i)
	{
		if (in[i] > 0)
			out[i] = 1;
		else if (in[i] < 0)
			out[i] = -1;
		else
			out[i] = 0;
	}
}

void vectormult(vector <myflo>& in, int scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] *= scalar;
}
void vectormult(vector <myflo>& in, myflo scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] *= scalar;
}
void vectormult(vector <int>& in, int scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] *= scalar;
}
void vectormult(vector <int>& in, myflo scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] *= scalar;
}

void vectordiv(vector <myflo>& in, int scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] /= scalar;
}
void vectordiv(vector <myflo>& in, myflo scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] /= scalar;
}
void vectordiv(vector <int>& in, int scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] /= scalar;
}
void vectordiv(vector <int>& in, myflo scalar)
{
#pragma omp parallel for schedule(static, 1) 
	for (int i = 0; i < in.size(); ++i)
		in[i] /= scalar;
}

void PeakFinder::findPeaks(vector <myflo> & x0, vector <int> & peakInds, bool includeEndpoints, int extrema)
{
	int minIdx = distance(x0.begin(), min_element(x0.begin(), x0.end()));
	int maxIdx = distance(x0.begin(), max_element(x0.begin(), x0.end()));

	myflo sel = (x0[maxIdx] - x0[minIdx]) / 4.0;
	int len0 = x0.size();

	if (extrema == -1) 
		vectormult(x0, extrema);

	vector <myflo> dx;
	diff(x0, dx, 1);
	replace(dx.begin(), dx.end(), (myflo)0.0, -PeakFinder::EPS);
	vector <myflo> dx0(dx.begin(), dx.end() - 1);
	vector <myflo> dx0_1(dx.begin() + 1, dx.end());
	vector <myflo> dx0_2;

	vectorElementsProduct(dx0, dx0_1, dx0_2);

	vector<int> ind;
	findIndicesLessThan(dx0_2, 0, ind); // Find where the derivative changes sign	
	vector <myflo> x;
	myflo leftMin;
	int minMagIdx;
	myflo minMag;
	
	if (includeEndpoints)
	{
		//x = [x0(1);x0(ind);x0(end)];	
		selectElementsFromIndices(x0, ind, x);
		x.insert(x.begin(), x0[0]);
		x.insert(x.end(), x0[x0.size() - 1]);
		//ind = [1;ind;len0];
		ind.insert(ind.begin(), 1);
		ind.insert(ind.end(), len0);
		minMagIdx = distance(x.begin(), min_element(x.begin(), x.end()));
		minMag = x[minMagIdx];
		//std::cout<<"Hola"<<std::endl;
		leftMin = minMag;
	}
	else
	{
		selectElementsFromIndices(x0, ind, x);
		if (x.size() > 2)
		{
			minMagIdx = distance(x.begin(), min_element(x.begin(), x.end()));
			minMag = x[minMagIdx];
			leftMin = x[0] < x0[0] ? x[0] : x0[0];
		}
	}

	int len = x.size();

	if (len > 2)
	{
		myflo tempMag = minMag;
		bool foundPeak = false;
		int ii;
		
		if (includeEndpoints)
		{
			// Deal with first point a little differently since tacked it on
			// Calculate the sign of the derivative since we tacked the first
			//  point on it does not neccessarily alternate like the rest.
			vector <myflo> xSub0(x.begin(), x.begin() + 3);//tener cuidado subvector
			vector <myflo> xDiff;//tener cuidado subvector
			diff(xSub0, xDiff, 1);

			vector<int> signDx;
			signVector(xDiff, signDx);

			if (signDx[0] <= 0) // The first point is larger or equal to the second
			{
				if (signDx[0] == signDx[1]) // Want alternating signs
				{
					x.erase(x.begin() + 1);
					ind.erase(ind.begin() + 1);
					len = len - 1;
				}
			}
			else // First point is smaller than the second
			{
				if (signDx[0] == signDx[1]) // Want alternating signs
				{
					x.erase(x.begin());
					ind.erase(ind.begin());
					len = len - 1;
				}
			}
		}
		
		//Skip the first point if it is smaller so we always start on maxima
		if (x[0] >= x[1])
			ii = 0;
		else
			ii = 1;

		//Preallocate max number of maxima
		myflo maxPeaks = ceil((myflo)len / 2.0);
		vector<int> peakLoc(maxPeaks, 0);
		vector <myflo> peakMag(maxPeaks, 0.0);
		int cInd = 1;
		int tempLoc;
		
		while (ii < len)
		{
			ii = ii + 1;//This is a peak
			//Reset peak finding if we had a peak and the next peak is bigger
			//than the last or the left min was small enough to reset.
			if (foundPeak)
			{
				tempMag = minMag;
				foundPeak = false;
			}

			//Found new peak that was lager than temp mag and selectivity larger
			//than the minimum to its left.

			if (x[ii - 1] > tempMag && x[ii - 1] > leftMin + sel)
			{
				tempLoc = ii - 1;
				tempMag = x[ii - 1];
			}

			//Make sure we don't iterate past the length of our vector
			if (ii == len)
				break; //We assign the last point differently out of the loop

			ii = ii + 1; // Move onto the valley

			//Come down at least sel from peak
			if (!foundPeak && tempMag > sel + x[ii - 1])
			{
				foundPeak = true; //We have found a peak
				leftMin = x[ii - 1];
				peakLoc[cInd - 1] = tempLoc; // Add peak to index
				peakMag[cInd - 1] = tempMag;
				cInd = cInd + 1;
			}
			else if (x[ii - 1] < leftMin) // New left minima
				leftMin = x[ii - 1];

		}

		// Check end point
		if (includeEndpoints)
		{
			if (x[x.size() - 1] > tempMag && x[x.size() - 1] > leftMin + sel)
			{
				peakLoc[cInd - 1] = len - 1;
				peakMag[cInd - 1] = x[x.size() - 1];
				cInd = cInd + 1;
			}
			else if (!foundPeak && tempMag > minMag)// Check if we still need to add the last point
			{
				peakLoc[cInd - 1] = tempLoc;
				peakMag[cInd - 1] = tempMag;
				cInd = cInd + 1;
			}
		}
		else if (!foundPeak)
		{
			myflo minAux = x0[x0.size() - 1] < x[x.size() - 1] ? x0[x0.size() - 1] : x[x.size() - 1];
			if (x[x.size() - 1] > tempMag && x[x.size() - 1] > leftMin + sel)
			{
				peakLoc[cInd - 1] = len - 1;
				peakMag[cInd - 1] = x[x.size() - 1];
				cInd = cInd + 1;
			}
			else if (!tempMag > minAux + sel)// Check if we still need to add the last point
			{
				peakLoc[cInd - 1] = tempLoc;
				peakMag[cInd - 1] = tempMag;
				cInd = cInd + 1;
			}
		}

		//Create output
		if (cInd > 0)
		{
			try
			{
				if (cInd < peakLoc.size())
				{
					vector<int> peakLocTmp(peakLoc.begin(), peakLoc.begin() + cInd - 1);
					selectElementsFromIndices(ind, peakLocTmp, peakInds);
				}
				else
				{
					vector<int> peakLocTmp(peakLoc.begin(), peakLoc.begin() + peakLoc.size() - 1);
					selectElementsFromIndices(ind, peakLocTmp, peakInds);
				}
			}
			catch (...)
			{
				//MessageBoxA(NULL, "Index is out of range. Programm continued running", "Error!", MB_ICONWARNING | MB_OK);
			}
		}
	}
}