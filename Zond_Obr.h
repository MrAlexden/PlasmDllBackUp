#pragma once

void levmarq(vector vec_Y_stepennoi, vector vec_X_stepennoi, vector& vParams, int Num_iter);

double error_func(vector vec_Y_stepennoi, vector vec_X_stepennoi, double A, double B, double C, double D);

void solve_axb_cholesky(Matrix ch, vector& delta, vector drvtv);

bool cholesky_decomp(Matrix& ch, Matrix Hessian);

