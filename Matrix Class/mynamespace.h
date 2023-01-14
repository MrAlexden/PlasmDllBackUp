#pragma once

#include <vector>

namespace myspace
{
#define derivative_step 0.01

#define sqr(a) (a)*(a)

#define MATRIX std::vector<std::vector <T>>
#define VECTOR std::vector<T>

	/*********************************************************/
	/*         Operators Extension for std::vector           */
	/*********************************************************/

	/* vector + vector */
	template <typename T>
	VECTOR operator + (const VECTOR & inl, const VECTOR & inr)
	{
		if (inl.size() != inr.size())
		{
			VECTOR v;
			return v;
		}

		VECTOR v(inl);

		for (int i = 0; i < inl.size(); ++i)
			v[i] += inr[i];

		return v;
	};

	/* vector - vector */
	template <typename T>
	VECTOR operator - (const VECTOR & inl, const VECTOR & inr)
	{
		if (inl.size() != inr.size())
		{
			VECTOR v;
			return v;
		}

		VECTOR v(inl);

		for (int i = 0; i < inl.size(); ++i)
			v[i] -= inr[i];

		return v;
	};

	/* vector += vector */
	template <typename T>
	bool operator += (VECTOR & inl, const VECTOR & inr)
	{
		if (inl.size() != inr.size())
			return false;

		for (int i = 0; i < inl.size(); ++i)
			inl[i] += inr[i];

		return true;
	};

	/* vector -= vector */
	template <typename T>
	bool operator -= (VECTOR & inl, const VECTOR & inr)
	{
		if (inl.size() != inr.size())
			return false;

		for (int i = 0; i < inl.size(); ++i)
			inl[i] -= inr[i];

		return true;
	};

	/* vector * scalar */
	template <typename T, typename TT>
	VECTOR operator * (const VECTOR & inl, const TT scalar)
	{
		VECTOR v(inl);

		for (int i = 0; i < inl.size(); ++i)
			v[i] *= (T)scalar;

		return v;
	};

	/* vector / scalar */
	template <typename T, typename TT>
	VECTOR operator / (const VECTOR & inl, const TT scalar)
	{
		VECTOR v(inl);

		for (int i = 0; i < inl.size(); ++i)
			v[i] /= (T)scalar;

		return v;
	};


	/*********************************************************/
	/*        My Matrix Class (child of std::vector)         */
	/*********************************************************/


	template <typename T>
	class matrix : public MATRIX
	{
	public:


		/*********************************************************/
		/*                    Constructors                       */
		/*********************************************************/

		/* дефолтный конструктор */
		matrix() {};

		/* конструктор копирования из матрицы */
		matrix(const matrix & in) :
			MATRIX(in.rows),
			rows(in.rows),
			cols(in.cols)
		{
			*this = in;
		};

		/* конструктор копирования из вектора */
		matrix(const VECTOR & in) :
			MATRIX(1),
			rows(1),
			cols(in.size())
		{
			this->SetSize(1, in.size());

			(*this)[0].assign(in.begin(), in.end());
		};

		/* конструктор с заданным размером */
		matrix(const size_t nrows, const size_t ncols, const T default_value = NULL) :
			MATRIX(nrows),
			rows(nrows),
			cols(ncols)
		{
			for (auto & outer : *this)
				outer.resize(ncols, default_value);
		};


		/*********************************************************/
		/*                     Destructors                       */
		/*********************************************************/

		/* дефолтный деструктор */
		~matrix()
		{
		};


		/*********************************************************/
		/*                   Public Methods                      */
		/*********************************************************/

		/* устанавливает размер матрицы */
		inline void SetSize(const size_t nrows, const size_t ncols)
		{
			SetRowsNumber(nrows);
			SetColsNumber(ncols);
			NullProperties();
		};

		/* возвращает текущую транспонированную матрицу (A^T) */
		matrix transpose() const
		{
			matrix m(cols, rows);

			for (int i = 0; i < cols; ++i)
				for (int j = 0; j < rows; ++j)
					m[i][j] = (*this)[j][i];

			return m;
		};

		/* возвращает определитель текущей матрицы */
		T det() const
		{
			/* текущая должна быть квадратной */
			if (rows == cols && rows > 0)
			{
				T determinant = NULL;

				switch (rows)
				{
				case 1:
					return (*this)[0][0];
					break;
				case 2:
					return (*this)[0][0] * (*this)[1][1] - (*this)[1][0] * (*this)[0][1];
					break;
				default:
					for (int i = 0, k = 1; i < rows; ++i)
					{
						determinant += k * (*this)[0][i] * (this->TakeAllExcept(0, i)).det();
						k = -k;
					}
					break;
				}

				return determinant;
			}

			return NULL; // определитель == 0 -> матрица вырожденная (невозможно посчитать обратную)
		};

		/* возвращает текущую перевернутую матрицу (A^-1) */
		matrix inverse() const
		{
			/* текущая должна быть квадратной */
			if (rows != cols || rows < 1)
			{
				matrix m;
				return m;
			}

			T determinant = this->det();

			/* текущая должна быть невырожденной */
			if (determinant == NULL)
			{
				matrix m;
				return m;
			}

			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = sqr(-1.0, i + j + 2) * this->TakeAllExcept(j, i).det() / determinant;

			return m;
		};

		/* возвращает значения запрошенной колонки в векторе */
		VECTOR GetColumn(const size_t num) const
		{
			if (num >= cols)
			{
				VECTOR v;
				return v;
			}

			return (this->transpose()).at(num);
		};

		/* возвращает диагональную матрицу данной матрицы */
		matrix diag()
		{
			/* текущая должна быть квадратной */
			if (rows != cols || rows < 1)
			{
				matrix m;
				return m;
			}

			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				m[i][i] = (*this)[i][i];

			return m;
		};


		/*********************************************************/
		/*                  Public Operators                     */
		/*********************************************************/

		/* matrix = matrix */
		matrix & operator = (const matrix & in)
		{ 
			this->SetSize(in.rows, in.cols);

			for (int i = 0; i < in.rows; ++i)
				(*this)[i].assign(in[i].begin(), in[i].end());

			CopyProperties(in);

			return *this;
		};

		/* matrix + matrix */
		matrix operator + (const matrix & in) const
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
			{
				matrix m;
				return m;
			}

			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] + in[i][j];

			return m;
		};

		/* matrix - matrix */
		matrix operator - (const matrix & in) const
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
			{
				matrix m;
				return m;
			}

			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] - in[i][j];

			return m;
		};

		/* matrix * matrix */
		matrix operator * (const matrix & in) const
		{
			/* кол-во колонок левой матрицы должно быть равно кол-ву строк правой */
			if (cols != in.rows)
			{
				matrix m;
				return m;
			}

			matrix m(rows, in.cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < in.cols; ++j)
					for (int k = 0; k < cols; ++k)
						m[i][j] += (*this)[i][k] * in[k][j];

			return m;
		};

		/* matrix / matrix */
		matrix operator / (const matrix & in) const
		{
			/* кол-во строк правой матрицы должно быть равно кол-ву колонок левой */
			if (in.rows != cols)
			{
				matrix m;
				return m;
			}

			matrix m(rows, in.cols);

			m = (*this) * in.inverse();

			return m;
		};

		/* matrix += matrix */
		matrix operator += (const matrix & in)
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
			{
				matrix m;
				return m;
			}

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					(*this)[i][j] += in[i][j];

			return *this;
		};

		/* matrix -= matrix */
		matrix operator -= (const matrix & in)
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
			{
				matrix m;
				return m;
			}

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					(*this)[i][j] -= in[i][j];

			return *this;
		};

		/* matrix *= matrix */
		matrix operator *= (const matrix & in)
		{
			/* кол-во колонок левой матрицы должно быть равно кол-ву строк правой */
			if (cols != in.rows)
			{
				matrix m;
				return m;
			}

			matrix m(rows, in.cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < in.cols; ++j)
					for (int k = 0; k < cols; ++k)
						m[i][j] += (*this)[i][k] * in[k][j];

			*this = m;

			return *this;
		};

		/* matrix /= matrix */
		matrix operator /= (const matrix & in)
		{
			/* кол-во строк правой матрицы должно быть равно кол-ву колонок левой */
			if (in.rows != cols)
			{
				matrix m;
				return m;
			}

			*this *= in.inverse();

			return *this;
		};

		/* matrix + scalar */
		template <typename TT>
		matrix operator + (const TT in) const
		{
			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] + (T)in;

			return m;
		};

		/* matrix - scalar */
		template <typename TT>
		matrix operator - (const TT in) const
		{
			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] - (T)in;

			return m;
		};

		/* matrix * scalar */
		template <typename TT>
		matrix operator * (const TT in) const
		{
			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] * (T)in;

			return m;
		};

		/* matrix / scalar */
		template <typename TT>
		matrix operator / (const TT in) const
		{
			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] / (T)in;

			return m;
		};

	private:


		/*********************************************************/
		/*                  Private Methods                      */
		/*********************************************************/

		/* устанавливает количество строк */
		inline void SetRowsNumber(const size_t nrows)
		{
			rows = nrows;
			this->resize(nrows);
		};

		/* устанавливает количество столбцов */
		inline void SetColsNumber(const size_t ncols)
		{
			cols = ncols;
			for (auto & outer : *this)
				outer.resize(ncols);
		};

		/* обнуляет сойства текущей матрицы */
		inline void NullProperties()
		{
		};

		/* копирует свойства из входной матрицы в текущую */
		inline void CopyProperties(const matrix & in)
		{
		};

		/* возвращает всю матрицу кроме данной строки и столбца */
		matrix TakeAllExcept(const size_t row, const size_t col) const
		{
			matrix m(rows - 1, cols - 1);

			for (int i = 0, ki = 0; i < rows; ++i)
			{
				if (i != row)
				{
					for (int j = 0, kj = 0; j < cols; ++j)
					{
						if (j != col)
						{
							m[ki][kj] = (*this)[i][j];
							kj++;
						}
					}
					ki++;
				}
			}

			return m;
		};

	public:

		size_t rows = 0;
		size_t cols = 0;

	private:
		
	};


	/*********************************************************/
	/*            Operators Extension for matrix             */
	/*********************************************************/

	/* vector * matrix */
	template <typename T>
	VECTOR operator * (const VECTOR & inl, const matrix <T> & inr)
	{
		/* размер вектора слева должен быть равен кол-ву строк матрицы справа */
		if (inl.size() != inr.rows)
		{
			VECTOR v;
			return v;
		}

		VECTOR v(inr.cols);

		for (int j = 0; j < inr.cols; ++j)
			for (int k = 0; k < inl.size(); ++k)
				v[j] += inl[k] * inr[k][j];

		return v;
	};


	/*********************************************************/
	/*                  Usefull Functions                    */
	/*********************************************************/

	/* возвращает единичную квадратную матрицу заданного размера */
	template <typename T>
	matrix <T> IdentityMatrix(const size_t dimension)
	{
		if (dimension <= NULL)
		{
			matrix <T> m;
			return m;
		}

		matrix <T> m(dimension, dimension);

		for (int i = 0; i < dimension; ++i)
			for (int j = 0; j < dimension; ++j)
				if (i == j) m[i][j] = (T)1;

		return m;
	};

	/* возвращает вектор размера size, в точку с индексом point кладет значение value, остальное нули */
	template <typename T>
	VECTOR unit_vector(const size_t size, const size_t point, const T value = 1)
	{
		VECTOR v(size, 0.0);
		v[point] = value;

		return v;
	};

	/* возвращает значение частной производной данной функции по данному параметру */
	template <typename T>
	T partial_derivative(_In_ const VECTOR & vParams,
						 _In_ const T x,
						 _In_ T (*fx)(const T, const VECTOR &),
						 _In_ int numofparam)
	{
		T f = fx(x, vParams), newx = x;
		VECTOR newvParams = vParams;

		switch (numofparam)
		{
		case -1:
			newx += derivative_step;
			break;
		default:
			if (numofparam >= 0)
				newvParams[numofparam] += derivative_step;
			break;
		}
		
		return (fx(newx, newvParams) - f) / derivative_step;
	};

	/* возвращает градиет функции по тем переменным, которые обьявлены как незафиксированные в vFixed */
	template <typename T>
	VECTOR gradient(_In_ const VECTOR & vParams,
					_In_ const std::vector <bool> & vFixed,
					_In_ const T x,
				    _In_ T (*fx)(const T, const VECTOR &))
	{
		if (vParams.size() != vFixed.size())
		{
			VECTOR v;
			return v;
		}

		int i, j, npar = 0;

		/* считаем по скольки параметрам создаем градиент */
		for (i = 0; i < vFixed.size(); ++i)
			if (!vFixed[i])
				++npar;

		VECTOR v(npar);

		/* заполняем вектор градиента */
		for (i = 0, j = 0; i < vFixed.size(); ++i)
			if (!vFixed[i])
				v[j++] = partial_derivative(vParams, x, fx, i);

		return v;
	};

	/* check whether a value valid or not */
	template <typename T>
	bool is_invalid(T val)
	{
		if (val == NAN ||
			val == INFINITY ||
			val == -INFINITY ||
			val != val)
			return true;

		return false;
	};

	/* возвращает норму вектора */
	template <typename T>
	T norm(const VECTOR & in, const unsigned int boost = 1)
	{
		T sum = 0;

		for (int i = 0 + (boost - 1); i < in.size() - (boost - 1); i += boost)
			sum += sqr(in[i]) * boost;

		return sqrt(sum);
	};

	/* возвращает квадратичное отклонение функции fx от сета данных vx, vy */
	template <typename T>
	T Chi_sqr(_In_ const VECTOR & vx,
			  _In_ const VECTOR & vy,
			  _In_ const VECTOR & vParams,
			  _In_ T (*fx)(const T, const VECTOR &),
			  _In_opt_ const unsigned int boost = 1) // коэффициент ускорения вычисления
	{
		if (vx.size() != vy.size() || vParams.empty())
			return (T)-1;

		T err = NULL, f;

		for (int i = 0 + (boost - 1); i < vy.size() - (boost - 1); i += boost)
		{	/* boost - 1 чтобы было 0 если boost == 1 */
			f = fx(vx[i], vParams) - vy[i];
			err += sqr(f) * boost;
		}

		return err;
	};

	/* возвращает ТОЧНУЮ нижне-треугольную матрицу Гессе (вторых частных производных) для заданной функции fx.
	учитываются только те параметры, которые обьявлены как незафиксированные в vFixed */
	template <typename T>
	matrix <T> HessianMatrix_accurate(_In_ const VECTOR & vParams,
									  _In_ const std::vector <bool> & vFixed,
									  _In_ const T x,	
									  _In_ T (*fx)(const T, const VECTOR &))
	{
		if (vParams.size() != vFixed.size())
		{
			matrix <T> m;
			return m;
		}

		int i, j, npar = 0;
		T f = fx(x, vParams);

		/* считаем по скольки параметрам создаем матрицу Гессе */
		for (i = 1; i < vFixed.size(); ++i)
			if (!vFixed[i])
				++npar;

		matrix <T> hessian(npar, npar);

		for (i = 0; i < npar; ++i)
		{
			for (j = i; j < npar; ++j)
			{
				if (i == j)
				{
					VECTOR vneg = vParams,
						vpos = vParams;

					vneg[j] -= derivative_step;
					vpos[j] += derivative_step;

					hessian[i][j] = (fx(x, vneg) - 2 * f + fx(x, vpos)) / (derivative_step * derivative_step);
				}
				else
				{
					VECTOR v1 = vParams,
						   v2 = vParams,
						   v3 = vParams,
						   v4 = vParams;

					v1[i] += derivative_step;
					v2[i] += derivative_step;
					v3[i] -= derivative_step;
					v4[i] -= derivative_step;

					v1[j] += derivative_step;
					v2[j] -= derivative_step;
					v3[j] += derivative_step;
					v4[j] -= derivative_step;

					hessian[i][j] = (fx(x, v1) - fx(x, v2) - fx(x, v3) + fx(x, v4)) / (4 * derivative_step * derivative_step);
					hessian[j][i] = hessian[i][j];
				}
			}
		}

		return hessian;
	};

	/* возвращает ПРИБЛИЖЕННУЮ нижне-треугольную матрицу Гессе (вторых частных производных) для заданной функции fx.
	учитываются только те параметры, которые обьявлены как незафиксированные в vFixed */
	template <typename T>
	matrix <T> HessianMatrix_approximate(_In_ const VECTOR & vParams,
										 _In_ const std::vector <bool> & vFixed,
										 _In_ const T x,
										 _In_ T (*fx)(const T, const VECTOR &))
	{
		if (vParams.size() != vFixed.size())
		{
			matrix <T> m;
			return m;
		}

		matrix <T> jacobian = gradient(vParams, vFixed, x, fx);

		return jacobian.transpose() * jacobian;
	};

	/* заполняет матрицу Гессе и вектор частных производных ошибки для метода LMfit */
	inline int make_matrixs_for_LM(_In_ const std::vector <double> & vx,
								   _In_ const std::vector <double> & vy,
								   _In_ const std::vector <double> & vParams,
								   _In_ const std::vector <bool> & vFixed,
								   _In_ const unsigned int npar,
								   _In_ double (*fx)(const double, const std::vector <double> &),
								   _Out_ matrix <double> & hessian,
								   _Out_ std::vector <double> & errpartdrvtv,
								   _In_opt_ const unsigned int boost = 1)
	{
		hessian.clear();
		errpartdrvtv.clear();

		hessian.SetSize(npar, npar);
		errpartdrvtv.resize(npar);

		/* Этот цикл пропускает создание матрицы Якоби размера "vx.size() X npar" и сразу
			заполняет матрицу Гессе произведением J^T * J */
		for (int x = 0 + (boost - 1); x < vx.size() - (boost - 1); x += boost)
		{
			/* заполнение матрицы Якоби производной в каждой точке */
			std::vector <double> jacobian = gradient(vParams, vFixed, vx[x], fx);

			/* заполнение матрицы Гессе J^T * J */
			for (int i = 0; i < npar; ++i)
			{
				for (int j = 0; j < npar; ++j)
					hessian[i][j] += jacobian[i] * jacobian[j] * boost;
				
				/* сразу же заполняем вектор (y - y(p)) * J */
				errpartdrvtv[i] += (vy[x] - fx(vx[x], vParams)) * jacobian[i] * boost;
			}
		}
		
		/* заполнение вектора частных производных Chi_sqr */
		/*for (int i = 0, j = 0; i < vFixed.size(); ++i)//работает
			if (!vFixed[i])
				errpartdrvtv[j++] = (curr_err -
					Chi_sqr(vx, vy, vParams + unit_vector<double>(vParams.size(), i, derivative_step), fx, boost))
					/ derivative_step;*/

		return 0;
	};

	/* аппроксимирует данные vx и vy функцией fx методом Левенберга-Марквардта и возвращает вектор с аппроксимированными переменными vParams */
	inline int LMfit(_In_ const std::vector <double> & vx,							// independent data
					 _In_ const std::vector <double> & vy,							// dependent data
					 _Inout_ std::vector <double> & vParams,						// vector of parameters we are searching for
					 _In_ const std::vector <bool> & vFixed,						// vector of param's fixed status 
					 _In_ const unsigned int niter,									// number of max iterations this method does
					 _In_ double (*fx)(const double, const std::vector <double> &),	// the function that the data will be approximated by
					 _In_opt_ const unsigned int boost = 10)							// koefficient which speeds up the computing in its value times
	{
		/* check for input mistakes */
		if (vx.size() <= boost
			|| vy.size() <= boost)
			return -1;
		if (vFixed.empty())
			return -1;
		if (vParams.size() != vFixed.size())
			return -1;
		if (niter <= NULL)
			return -1;

		unsigned int i, j, npar, outer;
		double lambda = 0.01,
			   nu = 10,
			   params_err = 1E-5,
			   minimum_err_step = 1E-12,
			   curr_err = Chi_sqr(vx, vy, vParams, fx, boost),
			   new_err = 1,
			   err_step = 1;
		bool stop = false;

		/* params counting, fixed params are not calculated */
		for (i = 0, npar = 0; i < vFixed.size(); ++i)
		{
			if (!vFixed[i])
				++npar;
			if (is_invalid(vParams[i]))
				vParams[i] = 1;
		}

		matrix <double> hessian, I = IdentityMatrix<double>(npar);
		std::vector <double> errpartdrvtv;

		for (outer = 0; outer < niter; ++outer, stop = false)
		{
			/* создание матрицы Гессе на данном шаге. Матрица Гессе заменена на J^T * J ->
			произведение транспонированной матрицы Якоби (градиента в данной точке сета vx)
			на нормальную матрицу Якоби */
			make_matrixs_for_LM(vx, vy, vParams, vFixed, npar, fx, hessian, errpartdrvtv, boost);

			for (; outer < niter && !stop;)
			{
				for (i = 0; i < npar; ++i)
					hessian[i][i] *= 1 + lambda;

				/* находим приращение к искомым параметрам */
				//std::vector <double> deltaParams = errpartdrvtv * (hessian + hessian.diag() * (1 + lambda)).inverse(); // Levenberg-Marquardt Method
				std::vector <double> deltaParams = errpartdrvtv * hessian.inverse();
				//std::vector <double> deltaParams = errpartdrvtv * hessian.inverse(); // Gauss-Newton Method
				//std::vector <double> deltaParams = errpartdrvtv * lambda; // Gradient Descent Method

				/* проверяем, были ли ошибки при расчете приращения к параметрам */
				if (deltaParams.empty())
					break;

				/* кладем прибавку к параметрам в вектор оригинальной длины, на случай если какие-то параметры были фиксированными */
				std::vector <double> deltaParamshandler(vFixed.size());
				for (i = 0, j = 0; i < vFixed.size(); ++i)
					if (!vFixed[i])
						deltaParamshandler[i] = deltaParams[j++];

				/* считаем отклонение от данных с новыми параметрами */
				new_err = Chi_sqr(vx, vy, vParams + deltaParamshandler, fx, boost);
				err_step = new_err - curr_err;

				if (norm(deltaParamshandler) <= params_err)
					stop = true;
				else
				{
					if (err_step < 0)
					{
						/* присваиваем параметрам их уточненное значение */
						vParams += deltaParamshandler;

						/* пересчитываем матрицу Гессе, Якоби и градиент с новыми параметрами */
						//make_matrixs_for_LM(vx, vy, vParams, vFixed, npar, fx, hessian, errpartdrvtv, boost);

						/* выходим из цикла так как приблизились к решению */
						stop = true;

						/* присваиваем текущую ошибку */
						curr_err = new_err;
						
						++outer;
					}
					else
					{
						/* увеличиваем лямбду т.к. мы отдалились от решения */
						//lambda *= (lambda < 1 / params_err) ? nu : 1;
						lambda *= nu;

						/* остаемся в цикле */
						stop = false;
					}
				}
			}
			
			/* уменьшаем лямбду т.к. мы приблизились к решению */
			//lambda /= (lambda > params_err) ? nu : 1;
			lambda /= nu;

			if (abs(err_step) <= minimum_err_step && stop)
				break;
		}

		MessageBoxA(NULL, ("Iterations: " + std::to_string(outer) + "\nError: " + std::to_string(-err_step)).c_str(), "Error!", MB_ICONINFORMATION | MB_OK);
		return 0;// 0.947 ... 0.641 // 1 ... 2 ... -0.066
	};

#undef MATRIX
#undef VECTOR
};