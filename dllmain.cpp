// dllmain.cpp : Определяет точку входа для приложения DLL.
#include "pch.h"

#define TOL 1e-30 /* smallest value allowed in cholesky_decomp() */
#define Pi 3.14159265

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

void levmarq(vector <double> &vec_Y, vector <double> &vec_X, vector <double> &vParams, UINT Num_iter)
{
	//double A = vParams[0], B = vParams[1], C = vParams[2], D = vParams[3];
	UINT x, i, j, it, ill, npar = ceil(vParams.size() % 2);
	double lambda = 0.001, up = 10, down = 1 / up, mult, err = 0, newerr = 0, derr = 0, target_derr = 10 ^ -12, new_C, new_D;

	Matrix Hessian(npar, npar, 0.0), ch(npar, npar, 0.0);
	vector <double> grad(npar), drvtv(npar), delta(npar);

	/* zeroing */
	for (i = 0; i < npar; i++)
	{
		grad[i] = 0.0;
		delta[i] = 0.0;
	}

	/* calculate the initial error ("chi-squared") */
	err = error_func(vec_Y, vec_X, vParams);

	/* main iteration */
	for (it = 0; it < Num_iter; it++)
	{
		/* zeroing */
		for (i = 0; i < npar; i++)
		{
			drvtv[i] = 0.0;
			for (j = 0; j < npar; j++)
			{
				Hessian(i, j) = 0;
				ch(i, j) = 0;
			}
		}

		/* calculate the approximation to the Hessian and the "derivative" drvtv */
		for (x = vec_Y.size() * 0.5; x < vec_Y.size(); x++)
		{
			grad[0] = dfxdC_STEP(vec_X[x], vParams[3]);
			grad[1] = dfxdD_STEP(vec_X[x], vParams[2], vParams[3]);

			for (i = 0; i < npar; i++)
			{
				drvtv[i] += (vec_Y[x] - fx_STEP(vec_X[x], vParams[0], vParams[1], vParams[2], vParams[3])) * grad[i];

				for (j = 0; j < npar; j++)
				{
					Hessian(i, j) += grad[i] * grad[j];
				}
			}
		}

		/*  make a step "delta."  If the step is rejected, increase lambda and try again */
		mult = 1 + lambda;
		ill = 1; /* ill-conditioned? */
		while (ill && (it < Num_iter))
		{
			for (i = 0; i < npar; i++) Hessian[i][i] = Hessian[i][i] * mult;

			ill = cholesky_decomp(ch, Hessian);

			if (!ill)
			{
				solve_axb_cholesky(ch, delta, drvtv);
				new_C = C + delta[0];
				new_D = D + delta[1];
				newerr = error_func(vec_Y, vec_X, A, B, new_C, new_D);
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
		C = new_C;
		D = new_D;
		err = newerr;
		lambda *= down;

		if ((!ill) && (-derr < target_derr)) break;
	}
	vParams[2] = C;
	vParams[3] = D;
}

/* calculate the error function (chi-squared) */
double error_func(vector <double> &vec_Y, vector <double> &vec_X, vector <double> &vParams)
{
	double res, e = 0;
	for (int i = 0; i < vec_Y.size(); i+=2)
	{
		res = fx_STEP(vec_X[i], vParams[0], vParams[1], vParams[2], vParams[3]) - vec_Y[i];
		e += res * res;
	}
	return e*2;
}

/* solve the equation Ax=b for a symmetric positive-definite matrix A,
   using the Cholesky decomposition A=LL^T.  The matrix L is passed in "ch".
   Elements above the diagonal are ignored. */
void solve_axb_cholesky(Matrix &ch, vector <double> &delta, vector <double> &drvtv)
{
	int i, j, npar = ch.getCols();
	double sum;

	/* solve L*y = b for y (where delta[] is used to store y) */

	for (i = 0; i < npar; i++)
	{
		sum = 0;
		for (j = 0; j < i; j++) sum += ch(i, j) * delta[j];
		delta[i] = (drvtv[i] - sum) / ch(i, i);
	}

	/* solve L^T*delta = y for delta (where delta[] is used to store both y and delta) */

	for (i = npar - 1; i >= 0; i--)
	{
		sum = 0;
		for (j = i + 1; j < npar; j++) sum += ch(j, i) * delta[j];
		delta[i] = (delta[i] - sum) / ch(i, i);
	}
}

/* This function takes a symmetric, positive-definite matrix "Hessian" and returns
   its (lower-triangular) Cholesky factor in "ch".  Elements above the
   diagonal are neither used nor modified.  The same array may be passed
   as both ch and Hessian, in which case the decomposition is performed in place. */
bool cholesky_decomp(Matrix &ch, Matrix &Hessian)
{
	int i, j, k, npar = ch.getCols();
	double sum;

	for (i = 0; i < npar; i++)
	{
		for (j = 0; j < i; j++)
		{
			sum = 0;
			for (k = 0; k < j; k++) sum += ch(i, k) * ch(j, k);
			ch(i, j) = (Hessian(i, j) - sum) / ch(j, j);
		}
		sum = 0;
		for (k = 0; k < i; k++) sum += ch(i, k) * ch(i, k);
		sum = Hessian(i, i) - sum;
		if (sum < TOL) return -1; /* not positive-definite */
		ch(i, i) = sqrt(sum);
	}
	return 0;
}

double fx_GAUSS(double x, double y0, double xc, double w, double A)
{
	return y0 + (A / (w * sqrt(Pi / 2))) * exp(-2 * pow(((x - xc) / w), 2));
}

double fx_STEP(double x, double A, double B, double C, double D)
{
	return A + B * x + C * exp(x / D);
}

double dfx_STEP(double x, double B, double C, double D)
{
	return B + (C * exp(x / D)) / D;
}

double dfxdC_STEP(double x, double D)
{
	return exp(x / D);
}

double dfxdD_STEP(double x, double C, double D)
{
	return -(C * x * exp(x / D)) / (D * D);
}
