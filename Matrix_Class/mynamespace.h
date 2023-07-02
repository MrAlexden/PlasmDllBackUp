#pragma once

#ifndef MYNAMESPACE
#define MYNAMESPACE

#include <vector>
#include <complex>
#include <thread>
#include <numbers> // std::numbers
#include "FFTW/fftw3.h"

namespace myspace
{

#ifndef M_PI
	#define M_PI std::numbers::pi
#endif 

#define sqr(a) (a)*(a)

#define MATRIX std::vector<std::vector <T>>
#define VECTOR std::vector<T>

#pragma region Vector&Matrix


	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	/*    При возврате std::vector из функции по значению
				      используется std::move, 
		   значит никаких указателей и ссылок не нужно       */
	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


	/*********************************************************/
	/*         Operators Extension for std::vector           */
	/*********************************************************/

	/* vector + vector */
	template <typename T>
		requires std::is_convertible_v<T, double>
	VECTOR operator + (const VECTOR & inl, const VECTOR & inr)
	{
		if (inl.size() != inr.size())
			return VECTOR();

		VECTOR v(inl);

		for (int i = 0; i < inl.size(); ++i)
			v[i] += inr[i];

		return v;
	};

	/* vector - vector */
	template <typename T>
		requires std::is_convertible_v<T, double>
	VECTOR operator - (const VECTOR & inl, const VECTOR & inr)
	{
		if (inl.size() != inr.size())
			return VECTOR();

		VECTOR v(inl);

		for (int i = 0; i < inl.size(); ++i)
			v[i] -= inr[i];

		return v;
	};

	/* vector += vector */
	template <typename T>
		requires std::is_convertible_v<T, double>
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
		requires std::is_convertible_v<T, double>
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
		requires std::is_convertible_v<T, double>
	VECTOR operator * (const VECTOR & inl, const TT scalar)
	{
		VECTOR v(inl);

		for (int i = 0; i < inl.size(); ++i)
			v[i] *= (T)scalar;

		return v;
	};

	/* vector / scalar */
	template <typename T, typename TT>
		requires std::is_convertible_v<T, double>
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
		requires std::is_convertible_v<T, double>
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
			for (int i = 0; i < this->size(); ++i)
				(*this)[i].resize(ncols, default_value);
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
		inline matrix transpose() const
		{
			matrix m(cols, rows);

			for (int i = 0; i < cols; ++i)
				for (int j = 0; j < rows; ++j)
					m[i][j] = (*this)[j][i];

			return m;
		};

		/* возвращает определитель текущей матрицы */
		inline T det() const
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
		inline matrix inverse() const
		{
			/* текущая должна быть квадратной */
			if (rows != cols || rows < 1)
				return matrix();

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
					m[i][j] = pow(-1.0, i + j + 2) * this->TakeAllExcept(j, i).det() / determinant;

			return m;
		};

		/* возвращает значения запрошенной колонки в векторе */
		inline VECTOR GetColumn(const size_t num) const
		{
			if (num >= cols)
				return VECTOR();

			return (this->transpose())[num];
		};

		/* возвращает диагональную матрицу данной матрицы */
		inline matrix diag()
		{
			/* текущая должна быть квадратной */
			if (rows != cols || rows < 1)
				return matrix();

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

		/* matrix + matrix
		матрицы должны быть одинакового размера */
		matrix operator + (const matrix & in) const
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
				return matrix();

			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] + in[i][j];

			return m;
		};

		/* matrix - matrix
		матрицы должны быть одинакового размера */
		matrix operator - (const matrix & in) const
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
				return matrix();

			matrix m(rows, cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					m[i][j] = (*this)[i][j] - in[i][j];

			return m;
		};

		/* matrix * matrix
		кол-во колонок левой матрицы должно быть равно кол-ву строк правой */
		matrix operator * (const matrix & in) const
		{
			/* кол-во колонок левой матрицы должно быть равно кол-ву строк правой */
			if (cols != in.rows)
				return matrix();

			matrix m(rows, in.cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < in.cols; ++j)
					for (int k = 0; k < cols; ++k)
						m[i][j] += (*this)[i][k] * in[k][j];

			return m;
		};

		/* matrix / matrix
		кол-во строк правой матрицы должно быть равно кол-ву колонок левой */
		matrix operator / (const matrix & in) const
		{
			/* кол-во строк правой матрицы должно быть равно кол-ву колонок левой */
			if (in.rows != cols)
				return matrix();

			matrix m(rows, in.cols);

			m = (*this) * in.inverse();

			return m;
		};

		/* matrix += matrix */
		matrix & operator += (const matrix & in)
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
				return matrix();

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					(*this)[i][j] += in[i][j];

			return *this;
		};

		/* matrix -= matrix
		матрицы должны быть одинакового размера */
		matrix & operator -= (const matrix & in)
		{
			/* матрицы должны быть одинакового размера */
			if (in.rows != rows || in.cols != cols)
				return matrix();

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					(*this)[i][j] -= in[i][j];

			return *this;
		};

		/* matrix *= matrix
		кол-во колонок левой матрицы должно быть равно кол-ву строк правой */
		matrix & operator *= (const matrix & in)
		{
			/* кол-во колонок левой матрицы должно быть равно кол-ву строк правой */
			if (cols != in.rows)
				return matrix();

			matrix m(rows, in.cols);

			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < in.cols; ++j)
					for (int k = 0; k < cols; ++k)
						m[i][j] += (*this)[i][k] * in[k][j];

			*this = m;

			return *this;
		};

		/* matrix /= matrix
		кол-во строк правой матрицы должно быть равно кол-ву колонок левой */
		matrix & operator /= (const matrix & in)
		{
			/* кол-во строк правой матрицы должно быть равно кол-ву колонок левой */
			if (in.rows != cols)
				return matrix();

			*this *= in.inverse();

			return *this;
		};

		/* matrix + scalar */
		template <typename TT>
			requires std::is_convertible_v<TT, T>
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
			requires std::is_convertible_v<TT, T>
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
			requires std::is_convertible_v<TT, T>
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
			requires std::is_convertible_v<TT, T>
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
			//TODO
		};

		/* копирует свойства из входной матрицы в текущую */
		inline void CopyProperties(const matrix & in)
		{
			//TODO
		};

		/* возвращает всю матрицу кроме данной строки и столбца */
		inline matrix TakeAllExcept(const size_t row, const size_t col) const
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

	/* vector * matrix
	размер вектора должен быть равен кол-ву строк матрицы */
	template <typename T>
		requires std::is_convertible_v<T, double>
	VECTOR operator * (const VECTOR & inl, const matrix <T> & inr)
	{
		/* размер вектора слева должен быть равен кол-ву строк матрицы справа */
		if (inl.size() != inr.rows)
			return VECTOR();

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
		requires std::is_convertible_v<T, double>
	matrix <T> IdentityMatrix(const size_t dimension)
	{
		if (dimension <= NULL)
			return matrix<T>();

		matrix <T> m(dimension, dimension);

		for (int i = 0; i < dimension; ++i)
			for (int j = 0; j < dimension; ++j)
				if (i == j) m[i][j] = (T)1;

		return m;
	};

	/* возвращает вектор размера size. В точку с индексом point кладет значение value, остальное нули */
	template <typename T>
		requires std::is_convertible_v<T, double>
	VECTOR unit_vector(const size_t size, const size_t point, const T value = 1)
	{
		VECTOR v(size, 0.0);
		v[point] = value;

		return v;
	};

#pragma endregion

#pragma region UsefullFuncs

#define derivative_step 0.001

	/* takes sub vector from iterator "begin" to iterator "end",
	accepts std::vector ONLY*/
	template<typename Iterator,
			 typename T = Iterator::value_type>
		requires (std::input_iterator<Iterator>
				  && std::is_same_v<Iterator, typename std::vector<T>::iterator>
				  && std::is_convertible_v<T, double>)
	inline std::vector <T> between(const Iterator begin, const Iterator end) // из итератора нельзя достать тип контейнера - тильт:(
	{
		if (begin >= end)
			return std::vector <T>();

		std::vector <T> out(end - begin);

		std::memmove(out.data(), &*begin, sizeof(T) * out.size());

		return out;
	};

	/* возвращает значение частной производной данной функции по данному параметру */
	template <typename T>
		requires std::is_convertible_v<T, double>
	inline T partial_derivative(_In_ const T x,
								_In_ const VECTOR& vP,
								_In_ T(*fx)(const T x, const VECTOR & params),
								_In_ int numofparam)
	{
		T f = std::move(fx(x, vP)), newx = x;
		VECTOR newvP = vP;

		switch (numofparam)
		{
		case -1:
			newx += derivative_step;
			break;
		default:
			if (numofparam >= 0)
				newvP[numofparam] += derivative_step;
			break;
		}

		return (fx(newx, newvP) - f) / derivative_step;
	};

	/* check whether a value valid or not */
	template <typename T>
		requires std::is_convertible_v<T, double>
	inline bool is_invalid(T val)
	{
		if (val == NAN ||
			val == INFINITY ||
			val == -INFINITY ||
			val != val)
			return true;

		return false;
	};

	/* возвращает градиет функции по тем переменным, которые обьявлены как незафиксированные в vFixed */
	/* TOO SLOW */
	template <typename T>
		requires std::is_convertible_v<T, double>
	inline VECTOR gradient(_In_ const T x,
						   _In_ const VECTOR& vP,
						   _In_ T(*fx)(const T, const VECTOR&),
						   _In_ const std::vector <bool>& vFixed = {})
	{
		int i, j, npar = 0;

		if (vFixed.empty())
			for (i = 0; i < vP.size(); ++i)
				vFixed[i] = false;
		if (vP.size() != vFixed.size())
			return VECTOR();

		/* считаем по скольки параметрам создаем градиент */
		for (i = 0; i < vFixed.size(); ++i)
			if (!vFixed[i])
				++npar;

		VECTOR v(npar);

		/* заполняем вектор градиента */
		for (i = 0, j = 0; i < vFixed.size(); ++i)
			if (!vFixed[i])
				v[j++] = partial_derivative(x, vP, fx, i);

		return v;
	};

	/* возвращает норму вектора */
	template <typename T>
		requires std::is_convertible_v<T, double>
	inline T norm(const VECTOR& in, const unsigned int myboost = 1)
	{
		T sum = 0;

		for (int i = 0 + (myboost - 1); i < in.size() - (myboost - 1); i += myboost)
			sum += sqr(in[i]) * myboost;

		return sqrt(sum);
	};

	/* возвращает квадратичное отклонение функции fx от сета данных vx, vy */
	/* TOO SLOW */
	template <typename T>
		requires std::is_convertible_v<T, double>
	inline T Chi_sqr(_In_ const VECTOR& vx,
					 _In_ const VECTOR& vy,
					 _In_ const VECTOR& vP,
					 _In_ T(*fx)(const T, const VECTOR&),
					 _In_opt_ const unsigned int myboost = 1) // коэффициент ускорения вычисления
	{
		if (vx.size() != vy.size() || vP.empty())
			return (T)-1;

		T err = NULL, f;

		for (size_t i = 0 + (myboost - 1); i < vy.size() - (myboost - 1); i += myboost)
			/* myboost -1 чтобы было 0 если myboost == 1 */
			f = fx(vx[i], vP) - vy[i], err += sqr(f) * myboost;

		return err;
	};

	template <typename T>
		requires std::is_convertible_v<T, double>
	inline VECTOR atan2(const VECTOR& r, const VECTOR& l)
	{
		if (r.size() != l.size())
			return VECTOR();

		VECTOR out(r.size());

		for (size_t i = 0; i < out.size(); ++i)
			out[i] = std::atan2(r[i], l[i]);

		return out;
	};

	template <typename T>
		requires std::is_convertible_v<T, double>
	inline VECTOR atan2(const std::vector<std::complex<T>> & in)
	{
		VECTOR out(in.size());

		for (size_t i = 0; i < out.size(); ++i)
			out[i] = std::atan2(in[i].imag(), in[i].real());

		return out;
	};

#pragma endregion

#pragma region MySlowLM

#define NOMYLM
#ifdef MYLM

	/* возвращает ТОЧНУЮ нижне-треугольную матрицу Гессе (вторых частных производных) для заданной функции fx.
	учитываются только те параметры, которые обьявлены как незафиксированные в vFixed */
	/* TOO SLOW */
	template <typename T>
	matrix <T> HessianMatrix_accurate(_In_ const VECTOR & vP,
									  _In_ const std::vector <bool> & vFixed,
									  _In_ const T x,	
									  _In_ T(*fx)(const T, const VECTOR &))
	{
		if (vP.size() != vFixed.size())
		{
			matrix <T> m;
			return m;
		}

		int i, j, npar = 0;
		T f = fx(x, vP);

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
					VECTOR vneg = vP,
						vpos = vP;

					vneg[j] -= derivative_step;
					vpos[j] += derivative_step;

					hessian[i][j] = (fx(x, vneg) - 2 * f + fx(x, vpos)) / (derivative_step * derivative_step);
				}
				else
				{
					VECTOR v1 = vP,
						   v2 = vP,
						   v3 = vP,
						   v4 = vP;

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
	/* TOO SLOW */
	template <typename T>
	matrix <T> HessianMatrix_approximate(_In_ const VECTOR & vP,
										 _In_ const std::vector <bool> & vFixed,
										 _In_ const T x,
										 _In_ T (*fx)(const T, const VECTOR &))
	{
		if (vP.size() != vFixed.size())
		{
			matrix <T> m;
			return m;
		}

		matrix <T> jacobian = gradient(x, vP, fx, vFixed);

		return jacobian.transpose() * jacobian;
	};

	/* заполняет матрицу Гессе и вектор частных производных ошибки для метода LMfit */
	/* TOO SLOW */
	inline int make_matrixs_for_LM(_In_ const std::vector <double> & vx,
								   _In_ const std::vector <double> & vy,
								   _In_ const std::vector <double> & vP,
								   _In_ const std::vector <bool> & vFixed,
								   _In_ const unsigned int npar,
								   _In_ double (*fx)(const double, const std::vector <double> &),
								   _Out_ matrix <double> & hessian,
								   _Out_ std::vector <double> & errpartdrvtv,
								   _In_opt_ const unsigned int myboost = 1)
	{
		hessian.clear();
		errpartdrvtv.clear();

		hessian.SetSize(npar, npar);
		errpartdrvtv.resize(npar, 0.0);

		/* Этот цикл пропускает создание матрицы Якоби размера "vx.size() X npar" и сразу
			заполняет матрицу Гессе произведением J^T * J */
		for (int x = 0 + (myboost - 1); x < vx.size() - (myboost - 1); x += myboost)
		{
			/* заполнение матрицы Якоби производной в каждой точке */
			std::vector <double> jacobian = gradient(vx[x], vP, fx, vFixed);

			/* заполнение матрицы Гессе J^T * J */
			for (int i = 0; i < npar; ++i)
			{
				for (int j = 0; j < npar; ++j)
					hessian[i][j] += jacobian[i] * jacobian[j] * myboost;
				
				/* сразу же заполняем вектор (y - y(p)) * J */
				errpartdrvtv[i] += (vy[x] - fx(vx[x], vP)) * jacobian[i] * myboost;
			}
		}
		
		/* заполнение вектора частных производных Chi_sqr */
		/*for (int i = 0, j = 0; i < vFixed.size(); ++i)//работает
			if (!vFixed[i])
				errpartdrvtv[j++] = (curr_err -
					Chi_sqr(vx, vy, vP + unit_vector<double>(vP.size(), i, derivative_step), fx, myboost))
					/ derivative_step;*/

		return 0;
	};

	/* аппроксимирует данные vx и vy функцией fx методом Левенберга-Марквардта и возвращает вектор с аппроксимированными переменными vP */
	/* TOO SLOW */
	inline int LMfit(_In_ const std::vector <double> & vx,							// independent data
					 _In_ const std::vector <double> & vy,							// dependent data
					 _Inout_ std::vector <double> & vP,						// vector of parameters we are searching for
					 _In_ const std::vector <bool> & vFixed,						// vector of param's fixed status 
					 _In_ const unsigned int niter,									// number of max iterations this method does
					 _In_ double (*fx)(const double, const std::vector <double> &),	// the function that the data will be approximated by
					 _In_opt_ const unsigned int myboost = 10)						// koefficient which speeds up the computing in its value times
	{
		/* check for input mistakes */
		if (vx.size() <= myboost
			|| vy.size() <= myboost)
			return -1;
		if (vFixed.empty())
			return -1;
		if (vP.size() != vFixed.size())
			return -1;
		if (niter <= NULL)
			return -1;

		unsigned int i, j, npar, outer;
		double lambda = 0.01,
			   nu = 10,
			   params_err = 1E-5,
			   minimum_err_step = 1E-12,
			   curr_err = Chi_sqr(vx, vy, vP, fx, myboost),
			   new_err = 1,
			   err_step = 1;
		bool stop = false;

		/* params counting, fixed params are not calculated */
		for (i = 0, npar = 0; i < vFixed.size(); ++i)
		{
			if (!vFixed[i])
				++npar;
			if (is_invalid(vP[i]))
				vP[i] = 1;
		}

		matrix <double> hessian, I = IdentityMatrix<double>(npar);
		std::vector <double> errpartdrvtv;

		for (outer = 0; outer < niter; ++outer, stop = false)
		{
			/* создание матрицы Гессе на данном шаге. Матрица Гессе заменена на J^T * J ->
			произведение транспонированной матрицы Якоби (градиента в данной точке сета vx)
			на нормальную матрицу Якоби */
			make_matrixs_for_LM(vx, vy, vP, vFixed, npar, fx, hessian, errpartdrvtv, myboost);

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
				new_err = Chi_sqr(vx, vy, vP + deltaParamshandler, fx, myboost);
				err_step = new_err - curr_err;

				if (norm(deltaParamshandler) <= params_err)
					stop = true;
				else
				{
					if (err_step < 0)
					{
						/* присваиваем параметрам их уточненное значение */
						vP += deltaParamshandler;

						/* пересчитываем матрицу Гессе, Якоби и градиент с новыми параметрами */
						//make_matrixs_for_LM(vx, vy, vP, vFixed, npar, fx, hessian, errpartdrvtv, myboost);

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

		//MessageBoxA(NULL, ("Iterations: " + std::to_string(outer) + "\nError: " + std::to_string(-err_step)).c_str(), "Error!", MB_ICONINFORMATION | MB_OK);
		return 0;// 0.947 ... 0.641 // 1 ... 2 ... -0.066
	};

#endif

#pragma endregion

#pragma region MyComplex

#define NOMYCOMPLEX
#ifdef MYCOMPLEX

	/****************************************************************/
	/* My complex number (unfortunately same speed as std::complex) */
	/****************************************************************/

	struct Complex
	{
		Complex() : r(0), i(0) {}
		Complex(double _r) : r(_r), i(0) {}
		Complex(double _r, double _i) : r(_r), i(_i) {}
		Complex(std::complex<double> _c) : r(_c.real()), i(_c.imag()) {}

		/* Complex = Complex */
		Complex & operator = (const Complex& right)
		{
			this->r = right.r;
			this->i = right.i;

			return *this;
		};
		/* Complex = Complex */
		Complex& operator = (Complex&& right)
		{
			*this = std::move(right);

			return *this;
		};

		/* Complex * Complex */
		Complex operator * (const Complex& right) const
		{
			return Complex(this->r * right.r - this->i * right.i, this->r * right.i + this->i * right.r);
		};
		/* Complex *= Complex */
		Complex & operator *= (const Complex& right)
		{
			double _r = this->r, _i = this->i;

			this->r = _r * right.r - _i * right.i;
			this->i = _r * right.i + _i * right.r;

			return *this;
		};

		/* Complex / Complex */
		Complex operator / (const Complex& right) const
		{
			double denominator = right.r * right.r + right.i * right.i;

			return Complex((this->r * right.r + this->i * right.i) / denominator, (this->i * right.r - this->r * right.i) / denominator);
		};
		/* Complex / double */
		Complex operator / (const double right) const
		{
			return Complex(this->r / right, this->i / right);
		};

		/* Complex += Complex */
		Complex & operator += (const Complex& right)
		{
			this->r = this->r + right.r;
			this->i = this->i + right.i;

			return *this;
		};

		/* Complex - Complex */
		Complex operator - (const Complex& right) const
		{
			return Complex(this->r - right.r, this->i - right.i);
		};

		Complex conj()
		{
			return Complex(this->r, this->i * -1);
		};

		double r = 0;
		double i = 0;
	};

#endif

#pragma endregion

#pragma region MyFFT

#define NOMYFFT
#ifdef MYFFT

	/****************************************************************/
	/*             My FFT (10 times slower than FFTW)               */
	/****************************************************************/

	class FFT
	{
#define TIMER0
#ifndef MYCOMPLEX
		typedef std::complex<double> Complex;
#endif
		
	public:

#ifdef MYCOMPLEX
		inline void fft(std::vector<std::complex<double>>& vec)
		{
			size_t n = vec.size();
			if (n <= 1)
				return;

			std::vector<Complex> vec_complex; vec_complex.assign(vec.begin(), vec.end());

			if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec_complex.data(), vec_complex.size());
			else
				fft_transform_bluestein(vec_complex.data(), vec_complex.size());

			for (size_t i = 0; i < n; ++i)
				vec[i] = std::move(std::complex<double>(vec_complex[i].r, vec_complex[i].i));
		};
#endif

		inline void fft(std::vector<Complex>& vec, bool make_half = false)
		{
			half = make_half;
			size_t n = vec.size();
			if (n <= 1)
				return;
			else if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec.data(), vec.size());
			else
				fft_transform_bluestein(vec.data(), vec.size());

			if (half) vec.resize(n / 2 + 1);
		};

		inline void fft(std::vector<double>& vec, bool make_half = false)
		{
			half = make_half;
			size_t n = vec.size();
			if (n <= 1)
				return;

			std::vector<Complex> vec_complex(vec.begin(), vec.end());

			if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec_complex.data(), vec_complex.size());
			else
				fft_transform_bluestein(vec_complex.data(), vec_complex.size());

			if (half) vec.resize(n / 2 + 1);

#ifdef MYCOMPLEX
			for (size_t i = 0; i < n; ++i)
				vec[i] = vec_complex[i].r;
#else
			for (size_t i = 0; i < vec.size(); ++i)
				vec[i] = vec_complex[i].real();
#endif
		};

		inline void fft(double * vec, size_t n)
		{
			if (n <= 1)
				return;

			std::vector<Complex> vec_complex(vec, vec + n);

			if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec_complex.data(), vec_complex.size());
			else
				fft_transform_bluestein(vec_complex.data(), vec_complex.size());

#ifdef MYCOMPLEX
			for (size_t i = 0; i < n; ++i)
				vec[i] = vec_complex[i].r;
#else
			for (size_t i = 0; i < n; ++i)
				vec[i] = vec_complex[i].real();
#endif
	};

#ifdef MYCOMPLEX
		inline void ifft(std::vector<std::complex<double>>& vec)
		{
			size_t n = vec.size();
			if (n <= 1)
				return;

			std::vector<Complex> vec_complex; vec_complex.assign(vec.begin(), vec.end());

			if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec_complex.data(), vec_complex.size(), true);
			else
				fft_transform_bluestein(vec_complex.data(), vec_complex.size(), true);

			for (size_t i = 0; i < n; ++i)
				vec[i] = std::move(std::complex<double>(vec_complex[i].r, vec_complex[i].i));
		};
#endif

		inline void ifft(std::vector<Complex>& vec)
		{
			size_t n = vec.size();
			if (n <= 1)
				return;
			else if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec.data(), vec.size(), true);
			else
				fft_transform_bluestein(vec.data(), vec.size(), true);
		};

		inline void ifft(std::vector<double>& vec)
		{
			size_t n = vec.size();
			if (n <= 1)
				return;

			std::vector<Complex> vec_complex(vec.begin(), vec.end());

			if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec_complex.data(), vec_complex.size(), true);
			else
				fft_transform_bluestein(vec_complex.data(), vec_complex.size(), true);

#ifdef MYCOMPLEX
			for (size_t i = 0; i < n; ++i)
				vec[i] = vec_complex[i].r;
#else
			for (size_t i = 0; i < n; ++i)
				vec[i] = vec_complex[i].real();
#endif
		};

		inline void ifft(double* vec, size_t n)
		{
			if (n <= 1)
				return;

			std::vector<Complex> vec_complex(vec, vec + n);

			if ((n & (n - 1)) == 0)  // Power of 2
				fft_transform_radix2(vec_complex.data(), vec_complex.size(), true);
			else
				fft_transform_bluestein(vec_complex.data(), vec_complex.size(), true);

#ifdef MYCOMPLEX
			for (size_t i = 0; i < n; ++i)
				vec[i] = vec_complex[i].r;
#else
			for (size_t i = 0; i < n; ++i)
				vec[i] = vec_complex[i].real();
#endif
		};

	private:

		/* Что-то не пойму почему, но без inline не компилится */
		inline void fft_transform_radix2(Complex* in, size_t n, bool invert = false, size_t* rev = nullptr)
		{
#ifdef TIMER1
			// Замер времени выполнения
			auto start = std::chrono::high_resolution_clock::now();
#endif
			Complex* wlen_pw = new Complex[n];

			if (rev == nullptr)
				// Permute vec by reversing the bits of addresses
				for (size_t i = 1, levels = log2(n); i < n; ++i)
				{
					// Reverse the bits of i
					size_t j = 0, ii = i;
					for (int k = 0; k < levels; ++k)
					{
						j = (j << 1) | (ii & 1);
						ii >>= 1;
					}

					if (j > i) std::swap(in[i], in[j]);
				}
			else
			{
				for (int i = 0; i < n; ++i)
					if (i < rev[i])
						std::swap(in[i], in[rev[i]]);
			}

#ifdef TIMER1
			auto end = std::chrono::high_resolution_clock::now();
			std::cout << "Время выполнения bit reverse: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " микросекунд." << std::endl;

			// Замер времени выполнения
			start = std::chrono::high_resolution_clock::now();
#endif 

			for (int len = 2, amp = invert ? 2 : -2; len <= n; len <<= 1)
			{
				double ang = amp * M_PI / len;
				int len2 = len >> 1;

				wlen_pw[0] = 1.0;
				for (int i = 1; i < len2; ++i)
					wlen_pw[i] = std::move(wlen_pw[i - 1] * Complex(cos(ang), sin(ang)));

				for (int i = 0; i < n; i += len)
				{
					Complex t,
						* pu = in + i,
						* pv = in + i + len2,
						* pu_end = in + i + len2,
						* pw = wlen_pw;
					for (; pu != pu_end; ++pu, ++pv, ++pw)
					{
						t = *pv * *pw;
						*pv = *pu - t;
						*pu += t;
					}
				}
			}

			delete[] wlen_pw;

#ifdef TIMER1
			end = std::chrono::high_resolution_clock::now();
			std::cout << "Время выполнения cooley-tukey: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " микросекунд." << std::endl;
#endif 
		};

		inline void fft_transform_bluestein(Complex* in, size_t n, bool invert = false)
		{
			double PI = invert ? M_PI : -M_PI;
			size_t m = 1 << ((int)log2(n << 1) + 1); // Find m = 2^k such that m >= 2 * n

			Complex* avec = new Complex[m], * bvec = new Complex[m], * omega = new Complex[n];
			size_t* rev = new size_t[m];

			std::thread T1([&]()
				{
					*avec = *in;
					*bvec = 1.0;
					*omega = 1.0;
					for (size_t i = 1; i < n; ++i)
					{
						omega[i] = std::move(std::polar(1.0, PI * ((i * i) % (n << 1)) / n));
						avec[i] = in[i] * omega[i];
#ifndef MYCOMPLEX
						bvec[i] = bvec[m - i] = std::move(std::conj(omega[i]));
#else
						bvec[i] = bvec[m - i] = std::move(omega[i].conj());
#endif
					}
				});
			std::thread T2([&]()
				{
					for (int i = 0, log_n = log2(m); i < m; ++i)
					{
						rev[i] = 0;
						for (int j = 0; j < log_n; ++j)
							if (i & (1 << j))
								rev[i] |= 1 << (log_n - 1 - j);
					}
				});
			T1.join(); T2.join();
			T1 = std::thread(&FFT::fft_transform_radix2, this, avec, m, false, rev);
			T2 = std::thread(&FFT::fft_transform_radix2, this, bvec, m, false, rev);
			T1.join(); T2.join();
			for (size_t i = 0; i < m; ++i)
				avec[i] *= bvec[i];
			fft_transform_radix2(avec, m, true, rev);

#ifdef TIMER2
			// Замер времени выполнения
			auto start = std::chrono::high_resolution_clock::now();
#endif

			n = half ? n : n >> 1 + 1;
			for (size_t i = 0; i < n; ++i)
				in[i] = std::move(avec[i] * omega[i] / double(m));

#ifdef TIMER2
			auto end = std::chrono::high_resolution_clock::now();
			std::cout << "Время выполнения последнего цикла: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " микросекунд." << std::endl;
#endif 

			delete[] avec, bvec, rev;
		};

		bool half = false;
	};

#endif

#pragma endregion

#pragma region FFTW

#define REAL 0
#define IMAG 1

	/* Computes the 1-D fast Fourier transform (BY FFTW) */
	inline void FFTW(double* in_r, size_t n)
	{
		fftw_complex* in_c = new fftw_complex[n];
		for (int i = 0; i < n; ++i) in_c[i][REAL] = in_r[i], in_c[i][IMAG] = 0.0;

		// create a DFT plan
		fftw_plan plan = fftw_plan_dft_1d(n, in_c, in_c, FFTW_FORWARD, FFTW_ESTIMATE);
		// execute the plan
		fftw_execute(plan);
		// do some cleaning
		fftw_destroy_plan(plan);
		fftw_cleanup();

		for (int i = 0; i < n; ++i) in_r[i] = in_c[i][REAL];

		delete[] in_c;
	};
	inline void FFTW(std::complex<double>* in_c, size_t n)
	{
		fftw_complex* in_c_ = new fftw_complex[n];
		for (int i = 0; i < n; ++i) in_c_[i][REAL] = in_c[i].real(), in_c_[i][IMAG] = in_c[i].imag();

		// create a DFT plan
		fftw_plan plan = fftw_plan_dft_1d(n, in_c_, in_c_, FFTW_FORWARD, FFTW_ESTIMATE);
		// execute the plan
		fftw_execute(plan);
		// do some cleaning
		fftw_destroy_plan(plan);
		fftw_cleanup();

		for (int i = 0; i < n; ++i) in_c[i] = std::complex<double>(in_c_[i][REAL], in_c_[i][IMAG]);

		delete[] in_c_;
	};

	/* Computes the 1-D inverse fast Fourier transform (BY FFTW) */
	inline void IFFTW(std::complex<double>* in_c, size_t n)
	{
		fftw_complex* in_c_ = new fftw_complex[n];
		for (int i = 0; i < n; ++i) in_c_[i][REAL] = in_c[i].real(), in_c_[i][IMAG] = in_c[i].imag();

		// create a DFT plan
		fftw_plan plan = fftw_plan_dft_1d(n, in_c_, in_c_, FFTW_BACKWARD, FFTW_ESTIMATE);
		// execute the plan
		fftw_execute(plan);
		// do some cleaning
		fftw_destroy_plan(plan);
		fftw_cleanup();

		for (int i = 0; i < n; ++i) in_c[i] = std::complex<double>(in_c_[i][REAL] / n, in_c_[i][IMAG] / n);

		delete[] in_c_;
	};

#undef REAL
#undef IMAG

#pragma endregion

#pragma region HilbertTransform

	/* Function to perform Hilber transform */
	inline void Hilbert(std::complex<double>* signal, const size_t length)
	{
		myspace::FFTW(signal, length);

		size_t h_length = length >> 1;

		for (std::complex<double>* it = signal + 1; it != signal + h_length; ++it)
			*it = std::complex<double>((*it).real() * 2, (*it).imag() * 2);

		if (length % 2 != 0 && length > 1)
			signal[h_length] = std::complex<double>(signal[h_length].real() * 2, signal[h_length].imag() * 2);

		for (std::complex<double>* it = signal + h_length + 1; it != signal + length; ++it)
			*it = std::complex<double>(0, 0);

		myspace::IFFTW(signal, length);
	};

	/* Unwraps the Phase array by eliminating discontinuities
		whose absolute values exceed either pi or 180 */
	template <typename T>
		requires std::is_convertible_v<T, double>
	inline void unwrap_phase(VECTOR& in)
	{
		for (size_t i = 1; i < in.size(); ++i)
			in[i] = in[i] - floor((in[i] - in[i - 1]) / (2 * M_PI) + 0.5) * 2 * M_PI;
	};

#pragma endregion

#pragma region MyFastLM

	/* solve the equation Ax=b for a symmetric positive-definite matrix A,
   using the Cholesky decomposition A=LL^T.  The matrix L is passed in "ch".
   Elements above the diagonal are ignored. */
	template <typename T>
		requires std::is_convertible_v<T, double>
	static inline void solve_axb_cholesky(_In_ const matrix<T>& ch,
										  _Out_ VECTOR& delta,
										  _In_ const VECTOR& drvtv)
	{
		int i, j, npar = ch.cols; // important not to be unsigned
		T sum;

		/* solve (ch)*y = drvtv for y (where delta[] is used to store y) */
		for (i = 0; i < npar; ++i)
		{
			sum = 0;
			for (j = 0; j < i; ++j)
				sum += ch[i][j] * delta[j];
			delta[i] = (drvtv[i] - sum) / ch[i][i];
		}

		/* solve (ch)^T*delta = y for delta (where delta[] is used to store both y and delta) */
		for (i = npar - 1; i >= 0; --i)
		{
			sum = 0;
			for (j = i + 1; j < npar; ++j)
				sum += ch[j][i] * delta[j];
			delta[i] = (delta[i] - sum) / ch[i][i];
		}
	};

	/* This function takes a symmetric, positive-definite matrix "Hessian" and returns
    its (lower-triangular) Cholesky factor in "ch".  Elements above the
    diagonal are neither used nor modified.  The same array may be passed
    as both ch and Hessian, in which case the decomposition is performed in place. */
	template <typename T>
		requires std::is_convertible_v<T, double>
	static inline bool cholesky_decomp(_Inout_ matrix<T>& ch,
									   _In_ const matrix<T>& Hessian)
	{
		int i, j, k, npar = ch.cols; // important not to be unsigned
		T sum;

		for (i = 0; i < npar; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				sum = 0;
				for (k = 0; k < j; ++k)
					sum += ch[i][k] * ch[j][k];
				ch[i][j] = (Hessian[i][j] - sum) / ch[j][j];
			}

			sum = 0;
			for (k = 0; k < i; ++k)
				sum += ch[i][k] * ch[i][k];
			sum = Hessian[i][i] - sum;
			if (sum < 10E-30) return true; /* not positive-definite */
			ch[i][i] = sqrt(sum);
		}
		return false;
	};

	template <typename T>
		requires std::is_convertible_v<T, double>
	inline int LevenbergMarquardt(_In_ const VECTOR& vX,							// independent data
								  _In_ const VECTOR& vY,							// dependent data
								  _Inout_ VECTOR& vP,								// vector of parameters we are searching for
								  _In_ T(*fx)(const T, const VECTOR&),				// the function that the data will be approximated by
								  _In_ const std::vector <bool>& vF = 
											std::vector <bool>(),					// vector of param's fixed status 
								  _In_ const unsigned int niter = 400,				// number of max iterations this method does
								  _In_opt_ const unsigned int myboost = 1)
	{
		int x, i, j, it, npar = 0; std::vector <bool> vFixed(vF);
		T lambda = 0.01f, up = 10, down = 1 / up, mult, err = 0, newerr = 0, derr = 0, target_derr = 1E-12;

		/* check for input mistakes */
		if (vX.size() != vY.size()
			|| vX.size() <= myboost
			|| vY.size() <= myboost)
			return -1;
		if (vFixed.empty())
			for (i = 0; i < vP.size(); ++i)
				vFixed.push_back(false);
		if (vP.size() != vFixed.size())
			//vFixed.resize(vP.size());
			return -2;
		if (niter == NULL)
			return -3;

		/* params counting, fixed params are not calculated */
		for (i = 0; i < vFixed.size(); ++i)
		{
			if (!vFixed[i])
				++npar;
			if (is_invalid(vP[i]))
				vP[i] = 1;
		}

		matrix <T> Hessian(npar, npar, 0.0), ch(npar, npar, 0.0);
		VECTOR grad(npar), drvtv(npar), delta(npar, 0.0), newvP = vP;

		/* calculate the initial error ("chi-squared") */
		err = Chi_sqr(vX, vY, vP, fx, myboost);

		/* main iteration */
		for (it = 0; it < niter; ++it)
		{
			/* zeroing */
			for (i = 0; i < npar; ++i)
			{
				drvtv[i] = 0.0;
				for (j = 0; j < npar; ++j)
					Hessian[i][j] = 0;
			}

			/* calculate the approximation to the Hessian and the "derivative" drvtv */
			for (x = 0 + (myboost - 1); x < vY.size() - (myboost - 1); x += myboost)
			{
				/* calculate gradient */
				for (i = 0, j = 0; i < vFixed.size(); ++i)
					if (!vFixed[i])
						grad[j++] = partial_derivative(vX[x], vP, fx, i);

				for (i = 0; i < npar; ++i)
				{
					drvtv[i] += (vY[x] - fx(vX[x], vP)) * grad[i] * myboost;

					for (j = 0; j < npar; ++j)
						Hessian[i][j] += grad[i] * grad[j] * myboost;
				}
			}

			/*  make a step "delta."  If the step is rejected, increase lambda and try again */
			mult = 1 + lambda;
			bool ill = true; /* ill-conditioned? */
			while (ill && (it < niter))
			{
				for (i = 0; i < npar; ++i)
					Hessian[i][i] = Hessian[i][i] * mult;

				ill = cholesky_decomp(ch, Hessian);

				if (!ill)
				{
					solve_axb_cholesky(ch, delta, drvtv);

					for (i = 0, j = 0; i < vFixed.size(); ++i)
						if (!vFixed[i])
							newvP[i] = vP[i] + delta[j++];

					newerr = Chi_sqr(vX, vY, newvP, fx, myboost);
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
				if (!vFixed[i])
					vP[i] = newvP[i];

			err = newerr;
			lambda *= down;

			if ((!ill) && (-derr < target_derr))
				//if ((!ill) && (abs(err) < target_derr))
				break;
		}

		return 0;
	};

#pragma endregion

#pragma region SavitskyGolay

	//! calculate savitzky golay coefficients.
	static inline void sg_coeff(_In_ const matrix<double>& c,
								_Out_ std::vector <double>& res,
								_In_ const size_t window,
								_In_ const size_t deg,
								_In_ const matrix<double>& A)
	{
		res.resize(window);
		for (int i = 0; i < window; ++i)
		{
			res[i] = c[0][0];
			for (int j = 1; j <= deg; ++j)
				res[i] += c[j][0] * A[i][j];
		}
		return;
	};

	/*! \brief savitzky golay smoothing.
	 *
	 * This method means fitting a polynome of degree 'deg' to a sliding window
	 * of width 2w+1 throughout the data.  The needed coefficients are
	 * generated dynamically by doing a least squares fit on a "symmetric" unit
	 * vector of size 2w+1, e.g. for w=2 b=(0,0,1,0,0). evaluating the polynome
	 * yields the sg-coefficients.  at the border non symmectric vectors b are
	 * used. */
	inline void sg_smooth(_Inout_ std::vector <double>& v_in,
						  _In_ const int width,
						  _In_ const int deg)
	{
		int i, j, k;

		if ((width < 1) || (deg <= 0) || (v_in.size() < (2 * width + 2)))
			return;

		const int window = 2 * width + 1, endidx = v_in.size() - 1, endidxv = v_in.size() + window - 1;
		std::vector <double> v(v_in.size() + window);

		// realising reflect method
		for (i = 0; i < width; ++i)
			v[i] = v_in[width - i], v[endidxv - i] = v_in[endidx - (width - i)];
		memcpy(&v[width], v_in.data(), sizeof(double) * v_in.size());
		memset(v_in.data(), 0.0, sizeof(double) * v_in.size());

		// generate input matrix for least squares fit
		matrix<double> A(window, deg + 1);
		for (i = 0; i < A.rows; ++i)
			for (j = 0; j < A.cols; ++j)
				A[i][j] = pow(i, j);

		matrix <double> a_transposed = A.transpose(), b(window, 1, 0.0), m; b[width][0] = 1.0; m = (a_transposed * A).inverse() * a_transposed * b;
		std::vector <double> c;
		sg_coeff(m, c, window, deg, A);
		
		//now loop over rest of data. reusing the "symmetric" coefficients.
		for (i = 0, k = 0; k < v_in.size(); ++i, ++k)
			for (j = 0; j < c.size(); ++j)
				v_in[k] += c[j] * v[i + j];

		return;
	};

#pragma endregion

#undef MATRIX
#undef VECTOR
};

#endif