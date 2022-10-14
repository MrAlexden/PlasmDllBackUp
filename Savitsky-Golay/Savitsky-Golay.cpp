#include "Savitsky-Golay.h"

//! comfortable array of doubles
using float_vect = vector <myflo>;
//! comfortable array of ints;
using int_vect = vector <int>;

/*! matrix class.
 *
 * This is a matrix class derived from a vector of float_vects.  Note that
 * the matrix elements indexed [row][column] with indices starting at 0 (c
 * style). Also note that because of its design looping through rows should
 * be faster than looping through columns.
 *
 * \brief two dimensional floating point array
 */
class float_mat : public std::vector<float_vect> {
private:
    //! disable the default constructor
    explicit float_mat() {};
    //! disable assignment operator until it is implemented.(не работает пока типа)
    float_mat& operator =(const float_mat&) { return *this; };
public:
    //! constructor with sizes
    float_mat(const size_t rows, const size_t cols, const myflo def = 0.0);
    //! copy constructor for matrix
    float_mat(const float_mat& m);
    //! copy constructor for vector
    float_mat(const float_vect& v);

    //! use default destructor
    // ~float_mat() {};

    //! get size
    size_t nr_rows(void) const { return size(); };
    //! get size
    size_t nr_cols(void) const { return front().size(); };
};



// constructor with sizes
float_mat::float_mat(const size_t rows, const size_t cols, const myflo defval)
    : std::vector<float_vect>(rows) {
    int i;
    for (i = 0; i < rows; ++i) {
        (*this)[i].resize(cols, defval);
    }
    if ((rows < 1) || (cols < 1)) {
        char buffer[1024];

        //sprintf_s(buffer, "cannot build matrix with %d rows and %d columns\n", rows, cols);
        //sgs_error(buffer);
    }
}

// copy constructor for matrix
float_mat::float_mat(const float_mat& m) : std::vector<float_vect>(m.size()) {

    float_mat::iterator inew = begin();
    float_mat::const_iterator iold = m.begin();
    for (/* empty */; iold < m.end(); ++inew, ++iold) {
        const size_t oldsz = iold->size();
        inew->resize(oldsz);
        const float_vect oldvec(*iold);
        *inew = oldvec;
    }
}

// copy constructor for vector
float_mat::float_mat(const float_vect& v)
    : std::vector<float_vect>(1) {

    const size_t oldsz = v.size();
    front().resize(oldsz);
    front() = v;
}

//////////////////////
// Helper functions //
//////////////////////

//! permute() orders the rows of A to match the integers in the index array.
void permute(float_mat& A, int_vect& idx)
{
    int_vect i(idx.size());
    int j, k;

    for (j = 0; j < A.nr_rows(); ++j) {
        i[j] = j;
    }

    // loop over permuted indices
    for (j = 0; j < A.nr_rows(); ++j) {
        if (i[j] != idx[j]) {

            // search only the remaining indices
            for (k = j + 1; k < A.nr_rows(); ++k) {
                if (i[k] == idx[j]) {
                    std::swap(A[j], A[k]); // swap the rows and
                    i[k] = i[j];     // the elements of
                    i[j] = idx[j];   // the ordered index.
                    break; // next j
                }
            }
        }
    }
}

/*! \brief Implicit partial pivoting.
 *
 * The function looks for pivot element only in rows below the current
 * element, A[idx[row]][column], then swaps that row with the current one in
 * the index map. The algorithm is for implicit pivoting (i.e., the pivot is
 * chosen as if the max coefficient in each row is set to 1) based on the
 * scaling information in the vector scale. The map of swapped indices is
 * recorded in swp. The return value is +1 or -1 depending on whether the
 * number of row swaps was even or odd respectively. */
static int partial_pivot(float_mat& A, const size_t row, const size_t col,
    float_vect& scale, int_vect& idx, myflo tol)
{
    if (tol <= 0.0)
        tol = TOL;

    int swapNum = 1;

    // default pivot is the current position, [row,col]
    int pivot = row;
    myflo piv_elem = fabs(A[idx[row]][col]) * scale[idx[row]];

    // loop over possible pivots below current
    int j;
    for (j = row + 1; j < A.nr_rows(); ++j) {

        const myflo tmp = fabs(A[idx[j]][col]) * scale[idx[j]];

        // if this elem is larger, then it becomes the pivot
        if (tmp > piv_elem) {
            pivot = j;
            piv_elem = tmp;
        }
    }

    if (pivot > row) {           // bring the pivot to the diagonal
        j = idx[row];           // reorder swap array
        idx[row] = idx[pivot];
        idx[pivot] = j;
        swapNum = -swapNum;     // keeping track of odd or even swap
    }
    return swapNum;
}

/*! \brief Perform backward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is upper
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the lower triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
static void lu_backsubst(float_mat & A, float_mat & a, bool diag = false)
{
    int r, c, k;

    for (r = (A.nr_rows() - 1); r >= 0; --r) {
        for (c = (A.nr_cols() - 1); c > r; --c) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] -= A[r][c] * a[c][k];
            }
        }
        if (!diag) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] /= A[r][r];
            }
        }
    }
}

/*! \brief Perform forward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is lower
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the upper triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
static void lu_forwsubst(float_mat & A, float_mat & a, bool diag = true)
{
    int r, k, c;
    for (r = 0; r < A.nr_rows(); ++r) {
        for (c = 0; c < r; ++c) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] -= A[r][c] * a[c][k];
            }
        }
        if (!diag) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] /= A[r][r];
            }
        }
    }
}

/*! \brief Performs LU factorization in place.
 *
 * This is Crout's algorithm (cf., Num. Rec. in C, Section 2.3).  The map of
 * swapped indeces is recorded in idx. The return value is +1 or -1
 * depending on whether the number of row swaps was even or odd
 * respectively.  idx must be preinitialized to a valid set of indices
 * (e.g., {1,2, ... ,A.nr_rows()}). */
static int lu_factorize(float_mat & A, int_vect & idx, myflo tol = TOL)
{
    if (tol <= 0.0)
        tol = TOL;

    if ((A.nr_rows() == 0) || (A.nr_rows() != A.nr_cols())) {
        //sgs_error("lu_factorize(): cannot handle empty "
        //           "or nonsquare matrices.\n");

        return -1;
    }

    float_vect scale(A.nr_rows());  // implicit pivot scaling
    int i, j;
    for (i = 0; i < A.nr_rows(); ++i) {
        myflo maxval = 0.0;
        for (j = 0; j < A.nr_cols(); ++j) {
            if (fabs(A[i][j]) > maxval)
                maxval = fabs(A[i][j]);
        }
        if (maxval == 0.0) {
            //sgs_error("lu_factorize(): zero pivot found.\n");
            return -1;
        }
        scale[i] = 1.0 / maxval;
    }

    int swapNum = 1;
    int c, r;
    for (c = 0; c < A.nr_cols(); ++c) {            // loop over columns
        swapNum *= partial_pivot(A, c, c, scale, idx, tol); // bring pivot to diagonal
        for (r = 0; r < A.nr_rows(); ++r) {      //  loop over rows
            int lim = (r < c) ? r : c;
            for (j = 0; j < lim; ++j) {
                A[idx[r]][c] -= A[idx[r]][j] * A[idx[j]][c];
            }
            if (r > c)
                A[idx[r]][c] /= A[idx[c]][c];
        }
    }
    permute(A, idx);
    return swapNum;
}

/*! \brief Solve a system of linear equations.
 * Solves the inhomogeneous matrix problem with lu-decomposition. Note that
 * inversion may be accomplished by setting a to the identity_matrix. */
void lin_solve(const float_mat & A, const float_mat & a, float_mat & b)
{
    int i, j;
    myflo tol = TOL;
    float_mat B(A);

    b.resize(a.nr_rows());
    for (i = 0; i < b.size(); ++i)
    {
        b[i].resize(a.nr_cols());
        for (j = 0; j < b[i].size(); ++j)
        {
            b[i][j] = a[i][j];
        }
    }

    int_vect idx(B.nr_rows());

    for (j = 0; j < B.nr_rows(); ++j) {
        idx[j] = j;  // init row swap label array
    }
    lu_factorize(B, idx, tol); // get the lu-decomp.
    permute(b, idx);          // sort the inhomogeneity to match the lu-decomp
    lu_forwsubst(B, b);       // solve the forward problem
    lu_backsubst(B, b);       // solve the backward problem
    return;
}
static float_mat lin_solve(const float_mat& A, const float_mat& a, myflo tol = TOL)
{
    float_mat B(A);
    float_mat b(a);
    int_vect idx(B.nr_rows());
    int j;

    for (j = 0; j < B.nr_rows(); ++j) {
        idx[j] = j;  // init row swap label array
    }
    lu_factorize(B, idx, tol); // get the lu-decomp.
    permute(b, idx);          // sort the inhomogeneity to match the lu-decomp
    lu_forwsubst(B, b);       // solve the forward problem
    lu_backsubst(B, b);       // solve the backward problem
    return b;
}

///////////////////////
// related functions //
///////////////////////

//! Returns the inverse of a matrix using LU-decomposition.
void invert(float_mat & A, float_mat & res)
{
    const int n = A.size();
    float_mat E(n, n, 0.0);
    float_mat B(A);
    int i;

    for (i = 0; i < n; ++i) {
        E[i][i] = 1.0;
    }

    lin_solve(B, E, res);

    return;
}
static float_mat invert(const float_mat& A)
{
    const int n = A.size();
    float_mat E(n, n, 0.0);
    float_mat B(A);
    int i;

    for (i = 0; i < n; ++i) {
        E[i][i] = 1.0;
    }

    return lin_solve(B, E);
}

//! returns the transposed matrix.
void transpose(float_mat & a, float_mat & res)
{
    int i, j;

    res.resize(a.nr_cols());
    for (i = 0; i < res.size(); ++i)
    {
        res[i].resize(a.nr_rows());
        for (j = 0; j < res[i].size(); ++j)
        {
            res[i][j] = a[j][i];
        }
    }
    return;
}
static float_mat transpose(const float_mat& a)
{
    float_mat res(a.nr_cols(), a.nr_rows());
    int i, j;

    for (i = 0; i < a.nr_rows(); ++i) {
        for (j = 0; j < a.nr_cols(); ++j) {
            res[j][i] = a[i][j];
        }
    }
    return res;
}

//! matrix multiplication.
float_mat operator *(const float_mat & a, const float_mat & b)
{
    float_mat res(a.nr_rows(), b.nr_cols());
    if (a.nr_cols() != b.nr_rows()) {
        //sgs_error("incompatible matrices in multiplication\n");
        return res;
    }

    int i, j, k;

    for (i = 0; i < a.nr_rows(); ++i) {
        for (j = 0; j < b.nr_cols(); ++j) {
            myflo sum(0.0);
            for (k = 0; k < a.nr_cols(); ++k) {
                sum += a[i][k] * b[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}

void matrixmult(float_mat & a, float_mat & b, float_mat & res)
{
    int i, j, k;

    if (a.nr_cols() != b.nr_rows()) {
        //sgs_error("incompatible matrices in multiplication\n");
        return;
    }

    res.resize(a.nr_rows());
#pragma omp parallel for schedule(static, 1)
    for (i = 0; i < res.size(); ++i)
    {
        res[i].resize(b.nr_cols());
#pragma omp parallel for schedule(static, 1)
        for (j = 0; j < res[i].size(); ++j)
        {
            myflo sum(0.0);
            for (k = 0; k < a.nr_cols(); ++k) 
            {
                sum += a[i][k] * b[k][j];
            }
            res[i][j] = sum;
        }
    }
    return;
}

//! calculate savitzky golay coefficients.
void sg_coeff(float_mat & c, vector <myflo> & res, const size_t window, const size_t deg, float_mat & A)
{
    res.resize(window);

#pragma omp parallel for schedule(static, 1) 
    for (int i = 0; i < window; ++i) {
        res[i] = c[0][0];
        for (int j = 1; j <= deg; ++j) {
            res[i] += c[j][0] * A[i][j];
        }
    }
    return;
}

/*! \brief savitzky golay smoothing.
 *
 * This method means fitting a polynome of degree 'deg' to a sliding window
 * of width 2w+1 throughout the data.  The needed coefficients are
 * generated dynamically by doing a least squares fit on a "symmetric" unit
 * vector of size 2w+1, e.g. for w=2 b=(0,0,1,0,0). evaluating the polynome
 * yields the sg-coefficients.  at the border non symmectric vectors b are
 * used. */
void sg_smooth(vector <myflo>& vorig, vector <myflo>& res, const int width, int deg)
{
    int i, j, k;

    if ((width < 1) || (deg < 0) || (vorig.size() < (2 * width + 2)))
        return;

    res.resize(vorig.size(), 0.0);

    const int window = 2 * width + 1;
    const int endidx = vorig.size() - 1;

    vector <myflo> v(vorig.size() + window);
    const int endidxv = v.size() - 1;
#pragma omp parallel for schedule(static, 1) private (i)
    for (i = 0; i < width; ++i)
    {
        v[i] = vorig[width - i];
        v[endidxv - i] = vorig[endidx - (width - i)];
    }
    memcpy(&v[width], vorig.data(), sizeof myflo * vorig.size());

    // do a regular sliding window average
    if (deg <= 0) {
        // handle border cases first because we need different coefficients
        for (i = 0; i < width; ++i) {
            const myflo scale = 1.0 / myflo(i + 1);
            const float_vect c1(width, scale);
            for (j = 0; j <= i; ++j) {
                res[i] += c1[j] * v[j];
                res[endidx - i] += c1[j] * v[endidx - j];
            }
        }

        // now loop over rest of data. reusing the "symmetric" coefficients.
        const myflo scale = 1.0 / myflo(window);
        const  float_vect c2(window, scale);

        for (i = 0; i <= (v.size() - window); ++i) {
            for (j = 0; j < window; ++j) {
                res[i + width] += c2[j] * v[i + j];
            }
        }
        return;
    }

    const size_t rows(window);
    const size_t cols(deg + 1);
    float_mat A(rows, cols);

    // generate input matrix for least squares fit
#pragma omp parallel for schedule(static, 1) private (i, j)
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            A[i][j] = pow(i, j);
        }
    }

    //float_mat c(invert(transpose(A) * A) * (transpose(A) * transpose(b)));
    float_mat AT(0, 0, 0.0), ATA_mult(0, 0, 0.0), ATA_inv(0, 0, 0.0), Ar(0, 0, 0.0);

    transpose(A, AT);
    matrixmult(AT, A, ATA_mult);
    invert(ATA_mult, ATA_inv);
    matrixmult(ATA_inv, AT, Ar);

    //handle border cases first because we need different coefficients
    //for (i = 0; i < width; ++i) 
    //{
    //    float_mat mc(0, 0, 0.0), b1(window, 1, 0.0);
    //    b1[i][0] = 1.0;
    //    vector <myflo> c1;

    //    matrixmult(Ar, b1, mc);

    //    sg_coeff(mc, c1, window, deg, A);

    //    for (j = 0; j < c1.size(); ++j) 
    //    {
    //        res[i] += c1[j] * v[j];
    //        res[endidx - i] += c1[j] * v[endidx - j];
    //    }
    //}

    // now loop over rest of data. reusing the "symmetric" coefficients.
    float_mat b2(window, 1, 0.0), mc(0, 0, 0.0), ATbT_mult(0, 0, 0.0);
    b2[width][0] = 1.0;
    vector <myflo> c2;

    matrixmult(Ar, b2, mc);

    sg_coeff(mc, c2, window, deg, A);

    for (i = 0, k = 0; k < res.size(); ++i, ++k)
    {
        for (j = 0; j < c2.size(); ++j)
        {
            res[k] += c2[j] * v[i + j];
        }
    }

    return;
}