#include "Lev-Marq.h"

/* calculate the error function (chi-squared) */
myflo error_func(vector <myflo>& vec_X, vector <myflo>& vec_Y, vector <myflo>& vParams, myflo(*fx)(myflo, vector <myflo>&))
{
	myflo res, e = 0;
	for (int x = 0 + (ef_koef - 1); x < vec_Y.size() - (ef_koef - 1); x += ef_koef)
	{
		res = fx(vec_X[x], vParams) - vec_Y[x];
		e += res * res * ef_koef;
	}
	return e;
}

/* solve the equation Ax=b for a symmetric positive-definite matrix A,
   using the Cholesky decomposition A=LL^T.  The matrix L is passed in "ch".
   Elements above the diagonal are ignored. */
void solve_axb_cholesky(Matrix& ch, vector <myflo>& delta, vector <myflo>& drvtv)
{
	int i, j, npar = ch.getCols();
	myflo sum;

	/* solve (ch)*y = drvtv for y (where delta[] is used to store y) */
#pragma omp parallel for schedule(static, 1) 
	for (i = 0; i < npar; ++i)
	{
		sum = 0;
		for (j = 0; j < i; ++j)
			sum += ch(i, j) * delta[j];
		delta[i] = (drvtv[i] - sum) / ch(i, i);
	}

	/* solve (ch)^T*delta = y for delta (where delta[] is used to store both y and delta) */
#pragma omp parallel for schedule(static, 1) 
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
bool cholesky_decomp(Matrix& ch, Matrix& Hessian)
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
myflo partial_derivative(myflo x, vector <myflo>& vParams, myflo(*fx)(myflo, vector <myflo>&), int numofparam)
{
	myflo h = 0.01, curfx = fx(x, vParams);
	vector <myflo> newvParams = vParams;

	if (numofparam == -1)
		x += h;
	else
		newvParams[numofparam] += h;

	return ((fx(x, newvParams)) - (curfx)) / h;
}

/* make linear approximation of given data */
vector <myflo> linear_fit(vector <myflo>& vec_X, vector <myflo>& vec_Y)
{
	vector <myflo> vParams = { 0, 0 };
	vector <bool> vFixed = { false, false };

	levmarq(vec_X, vec_Y, vParams, vFixed, Num_iter, fx_LINE);

	return vParams;
}

void levmarq(vector <myflo> & vec_X,		// independent data
	vector <myflo> & vec_Y,					// dependent data
	vector <myflo> & vParams,				// vector of parameters we are searching for
	vector <bool> & vFixed,					// vector of param's fixed status 
	unsigned int niter,						// number of max iterations this method does
	myflo (*fx)(myflo, vector <myflo>&))    // the function that the data will be approximated by
{
	int x, i, j, it, npar = 0;
	myflo lambda = 0.01, up = 10, down = 1 / up, mult, err = 0, newerr = 0, derr = 0, target_derr = TOL;

	/* check for input mistakes */
	if (vFixed.empty())
		for (i = 0; i < vParams.size(); ++i)
			vFixed[i] = false;
	if (niter == NULL)
		niter = 100;

	/* params counting, fixed params are not calculated */
	if (vParams.size() == vFixed.size())
	{
		for (i = 0; i < vFixed.size(); ++i)
			if (vFixed[i] == false)
				npar++;
	}
	else return;

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
					j++;
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
						j++;
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
				it++;
			}
		}

		for (i = 0; i < vFixed.size(); ++i)
			if (vFixed[i] == false)
				vParams[i] = newvParams[i];

		err = newerr;
		lambda *= down;

		if ((!ill) && (/*было -derr<...,стало:*/abs(err) < target_derr)) break;
	}
}