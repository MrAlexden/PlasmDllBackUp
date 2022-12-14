#include "Lev-Marq.h"

/* calculate the error function (chi-squared) */
inline myflo error_func(_In_ const vector <myflo>& vec_X, 
						_In_ const vector <myflo>& vec_Y, 
						_In_ const vector <myflo>& vParams, 
						_In_ myflo(*fx)(myflo, const vector <myflo>&))
{
	myflo res, err = 0;

	for (int x = 0 + (ef_koef - 1); x < vec_Y.size() - (ef_koef - 1); x += ef_koef)
		/* ef_koef - 1 ????? ???? 0 ???? ef_koef == 0*/
	{
		res = fx(vec_X[x], vParams) - vec_Y[x];
		err += res * res * ef_koef;
	}

	return err;
}

/* solve the equation Ax=b for a symmetric positive-definite matrix A,
   using the Cholesky decomposition A=LL^T.  The matrix L is passed in "ch".
   Elements above the diagonal are ignored. */
inline void solve_axb_cholesky(_In_ Matrix& ch, 
							   _Out_ vector <myflo>& delta, 
							   _In_ const vector <myflo>& drvtv)
{
	int i, j, npar = ch.getCols();
	myflo sum;

	/* solve (ch)*y = drvtv for y (where delta[] is used to store y) */ 
	for (i = 0; i < npar; ++i)
	{
		sum = 0;
		for (j = 0; j < i; ++j)
			sum += ch(i, j) * delta[j];
		delta[i] = (drvtv[i] - sum) / ch(i, i);
	}

	/* solve (ch)^T*delta = y for delta (where delta[] is used to store both y and delta) */
	for (i = npar - 1; i >= 0; --i)
	{
		sum = 0;
		for (j = i + 1; j < npar; ++j)
			sum += ch(j, i) * delta[j];
		delta[i] = (delta[i] - sum) / ch(i, i);
	}
}

/* This function takes a symmetric, positive-definite matrix "Hessian" and returns
   its (lower-triangular) Cholesky factor in "ch".  Elements above the
   diagonal are neither used nor modified.  The same array may be passed
   as both ch and Hessian, in which case the decomposition is performed in place. */
inline bool cholesky_decomp(_Inout_ Matrix& ch, 
							_In_ Matrix& Hessian)
{
	int i, j, k, npar = ch.getCols();
	myflo sum;

	for (i = 0; i < npar; ++i)
	{
		for (j = 0; j < i; ++j)
		{
			sum = 0;
			for (k = 0; k < j; ++k)
				sum += ch(i, k) * ch(j, k);
			ch(i, j) = (Hessian(i, j) - sum) / ch(j, j);
		}

		sum = 0;
		for (k = 0; k < i; ++k)
			sum += ch(i, k) * ch(i, k);
		sum = Hessian(i, i) - sum;
		if (sum < TOL) return true; /* not positive-definite */
		ch(i, i) = sqrt(sum);
	}
	return false;
}

/* calculate partial derivative of given func in particular point */
myflo partial_derivative(_In_ myflo x, 
						 _In_ const vector <myflo> & vParams, 
						 _In_ myflo(*fx)(myflo, const vector <myflo> &), 
						 _In_ int numofparam)
{
	myflo h = 0.01f, curfx = fx(x, vParams);
	vector <myflo> newvParams = vParams;

	if (numofparam == -1)
		x += h;
	else
		newvParams[numofparam] += h;

	return ((fx(x, newvParams)) - (curfx)) / h;
}

/* make linear approximation of given data */
vector <myflo> linear_fit(_In_ const vector <myflo> & vec_X, 
						  _In_ const vector <myflo> & vec_Y)
{
	vector <myflo> vParams = { 0.0f, 0.0f };
	vector <bool> vFixed = { false, false };

	levmarq(vec_X, vec_Y, vParams, vFixed, Num_iter, fx_LINE);

	return vParams;
}

void levmarq(_In_ const vector <myflo> & vec_X,					// independent data
			 _In_ const vector <myflo> & vec_Y,					// dependent data
			 _Inout_ vector <myflo> & vParams,					// vector of parameters we are searching for
			 _In_ vector <bool> & vFixed,						// vector of param's fixed status 
			 _In_ unsigned int niter,							// number of max iterations this method does
			 _In_ myflo (*fx)(myflo, const vector <myflo> &))	// the function that the data will be approximated by
{
	int x, i, j, it, npar = 0;
	myflo lambda = 0.01f, up = 10, down = 1 / up, mult, err = 0, newerr = 0, derr = 0, target_derr = TOL;

	/* check for input mistakes */
	if (vFixed.empty())
		for (i = 0; i < vParams.size(); ++i)
			vFixed[i] = false;
	if (vParams.size() != vFixed.size())
			vFixed.resize(vParams.size());
	if (niter == NULL || niter >= 400)
		niter = 400;

	/* params counting, fixed params are not calculated */
	for (i = 0; i < vFixed.size(); ++i)
		if (vFixed[i] == false)
			++npar;

	Matrix Hessian(npar, npar, 0.0), ch(npar, npar, 0.0);
	vector <myflo> grad(npar), drvtv(npar), delta(npar, 0.0), newvParams = vParams;

	/* calculate the initial error ("chi-squared") */
	err = error_func(vec_X, vec_Y, vParams, fx);

	/* main iteration */
	for (it = 0; it < niter; ++it)
	{
		/* zeroing */
		for (i = 0; i < npar; ++i)
		{
			drvtv[i] = 0.0;
			for (j = 0; j < npar; ++j)
			{
				Hessian(i, j) = 0;
			}
		}

		/* calculate the approximation to the Hessian and the "derivative" drvtv */
		for (x = 0 + (ef_koef - 1); x < vec_Y.size() - (ef_koef - 1); x += ef_koef)
		{
			/* calculate gradient */
			for (i = 0, j = 0; i < vFixed.size(); ++i)
			{
				if (vFixed[i] == false)
				{
					grad[j] = partial_derivative(vec_X[x], vParams, fx, i);
					++j;
				}
			}

			for (i = 0; i < npar; ++i)
			{
				drvtv[i] += (vec_Y[x] - fx(vec_X[x], vParams)) * grad[i] * ef_koef;

				for (j = 0; j < npar; ++j)
				{
					Hessian(i, j) += grad[i] * grad[j] * ef_koef;
				}
			}
		}

		/*  make a step "delta."  If the step is rejected, increase lambda and try again */
		mult = 1 + lambda;
		bool ill = true; /* ill-conditioned? */
		while (ill && (it < niter))
		{
			for (i = 0; i < npar; ++i)
				Hessian(i, i) = Hessian(i, i) * mult;

			ill = cholesky_decomp(ch, Hessian);

			if (!ill)
			{
				solve_axb_cholesky(ch, delta, drvtv);

				for (i = 0, j = 0; i < vFixed.size(); ++i)
				{
					if (vFixed[i] == false)
					{
						newvParams[i] = vParams[i] + delta[j];
						++j;
					}
				}

				newerr = error_func(vec_X, vec_Y, newvParams, fx);
				derr = newerr - err;
				ill = (derr > 0);
			}

			if (ill)
			{
				mult = (1 + lambda * up) / (1 + lambda);
				lambda *= up;
				++it;
			}
		}

		for (i = 0; i < vFixed.size(); ++i)
			if (vFixed[i] == false)
				vParams[i] = newvParams[i];

		err = newerr;
		lambda *= down;

		if ((!ill) && (/*???? -derr<...,?????:*/abs(err) < target_derr)) break;
	}
}