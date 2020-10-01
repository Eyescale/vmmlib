/*
 * Copyright (c) 2006-2016, Visualization and Multimedia Lab,
 *                          University of Zurich <http://vmml.ifi.uzh.ch>,
 *                          Eyescale Software GmbH,
 *                          Blue Brain Project, EPFL
 *
 * This file is part of VMMLib <https://github.com/VMML/vmmlib/>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution.  Neither the name of the Visualization and Multimedia
 * Lab, University of Zurich nor the names of its contributors may be used to
 * endorse or promote products derived from this software without specific prior
 * written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __VMML__MATRIX__HPP__
#define __VMML__MATRIX__HPP__

#include <vmmlib/types.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace vmml
{
/** Matrix with R rows and C columns of element type T */
template <size_t R, size_t C, typename T>
class Matrix
{
public:
    /**
     * Construct a zero-initialized matrix.
     * Square matrices are initialized as identity.
     */
    Matrix();

    /**
     * Construct a matrix with default values.
     * Missing data is zero-initialized. Additional data is ignored.
     */
    Matrix(const T* begin, const T* end);

    /**
     * Construct a matrix with default values.
     * Missing data is zero-initialized. Additional data is ignored.
     */
    explicit Matrix(const std::vector<T>& data);

    /**
     * Copy-construct a matrix.
     */
    Matrix(const Matrix<R, C, T>& source);

    /**
     * Copy-construct a matrix.
     * Missing data is zero-initialized. Additional data is ignored.
     */
    template <size_t P, size_t Q>
    Matrix(const Matrix<P, Q, T>& source);

    /**
     * Construct a new 4x4 transformation matrix from a rotation quaternion and
     * a translation vector.
     */
    template <size_t S>
    Matrix(
        const Quaternion<T>& rotation, const vector<S, T>& translation,
        typename std::enable_if<R == S + 1 && C == S + 1 && S == 3>::type* = 0);

    /**
     * Construct a new transformation matrix from an eye position, lookat
     * position and up vector, following convention from gluLookAt().
     */
    template <size_t S>
    Matrix(
        const vector<S, T>& eye, const vector<S, T>& lookat,
        const vector<S, T>& up,
        typename std::enable_if<R == S + 1 && C == S + 1 && S == 3>::type* = 0);

    /** @name matrix-matrix operations */
    //@{
    /** @return true if both matrices are equal. */
    bool operator==(const Matrix& other) const;

    /** @return true if both matrices are not equal. */
    bool operator!=(const Matrix& other) const;

    /** @return true if both matrices are equal within the given tolerance. */
    bool equals(const Matrix& other,
                T tolerance = std::numeric_limits<T>::epsilon()) const;

    /** Set this to the product of the two matrices (left_RxP * right_PxC) */
    template <size_t P>
    Matrix<R, C, T>& multiply(const Matrix<R, P, T>& left,
                              const Matrix<P, C, T>& right);

    /** @return matrix_RxP = (this) matrix * other matrix_CxP; */
    template <size_t P>
    Matrix<R, P, T> operator*(const Matrix<C, P, T>& other) const;

    /** Multiply two square matrices */
    template <size_t O, size_t P, typename = typename std::enable_if<
                                      R == C && O == P && R == O>::type>
    Matrix<R, C, T>& operator*=(const Matrix<O, P, T>& right);

    /** Element-wise addition of two matrices */
    Matrix operator+(const Matrix& other) const;

    /** Element-wise substraction of two matrices */
    Matrix operator-(const Matrix& other) const;

    /** Element-wise addition of two matrices */
    void operator+=(const Matrix& other);

    /** Element-wise substraction of two matrices */
    void operator-=(const Matrix& other);
    //@}

    /** @name matrix-vector operations */
    //@{
    /** Transform column vector by matrix ( res = matrix * vec ) */
    vector<R, T> operator*(const vector<C, T>& other) const;
    //@}

    /** @name Data access */
    //@{
    /** @return the element at the given row and column. */
    T& operator()(size_t rowIndex, size_t colIndex);

    /** @return the element at the given row and column. */
    T operator()(size_t rowIndex, size_t colIndex) const;

    /** @return the pointer to the data storage in column-major order. */
    const T* data() const { return array; }
    /** @return the sub matrix of size OxP at the given start indices. */
    template <size_t O, size_t P>
    Matrix<O, P, T> getSubMatrix(
        size_t rowOffset, size_t colOffset,
        typename std::enable_if<O <= R && P <= C>::type* = 0) const;

    /** Set the sub matrix of size OxP at the given start indices. */
    template <size_t O, size_t P>
    typename std::enable_if<O <= R && P <= C>::type* setSubMatrix(
        const Matrix<O, P, T>& sub_matrix, size_t rowOffset, size_t colOffset);

    /** Assign the given matrix. */
    const Matrix& operator=(const Matrix<R, C, T>& source_);

    /**
     * Assign the given matrix.
     *
     * Remaining data is zero-filled. Additional data is ignored.
     */
    template <size_t P, size_t Q>
    const Matrix& operator=(const Matrix<P, Q, T>& source_);

    /**
     * Assign the given data.
     *
     * Remaining data is zero-filled. Additional data is ignored.
     */
    void operator=(const std::vector<T>& data);

    /** @return the negated matrix of this matrix. */
    Matrix<R, C, T> operator-() const;

    /** @return a vector of the given column. */
    vector<R, T> getColumn(size_t columnIndex) const;

    /** Set the given column. */
    void setColumn(size_t index, const vector<R, T>& column);

    /** @return a vector of the given row. */
    vector<C, T> getRow(size_t index) const;

    /** Set the given row. */
    void setRow(size_t index, const vector<C, T>& row);

    /** @return the translation vector (of a 3x3 or 4x4 matrix) */
    vector<C - 1, T> getTranslation() const;

    /** Set the translation vector (of a 3x3 or 4x4 matrix) */
    Matrix<R, C, T>& setTranslation(const vector<C - 1, T>& t);

    /**
     * Decompose a 4x4 transformation matrix to eye position, lookAt position
     * and up vector.
     */
    template <size_t S>
    void getLookAt(vector<S, T>& eye, vector<S, T>& lookAt, vector<S, T>& up,
                   typename std::enable_if<R == S + 1 && C == S + 1 &&
                                           S == 3>::type* = 0) const;
    //@}

    /**
     * Compute and return the inverted matrix of this matrix.
     *
     * Sets values to quiet_NaN if not invertible.
     */
    Matrix<R, C, T> inverse() const;

    template <size_t O, size_t P>
    typename std::enable_if<O == P && R == C && O == R && R >= 2>::type*
        getAdjugate(Matrix<O, P, T>& adjugate) const;

    template <size_t O, size_t P>
    typename std::enable_if<O == P && R == C && O == R && R >= 2>::type*
        getCofactors(Matrix<O, P, T>& cofactors) const;

    // returns the determinant of a square matrix R-1, C-1
    template <size_t O, size_t P>
    T getMinor(Matrix<O, P, T>& minor_, size_t row_to_cut, size_t col_to_cut,
               typename std::enable_if<O == R - 1 && P == C - 1 && R == C &&
                                       R >= 2>::type* = 0) const;

    /** @name Transformations on 4*4 matrices */
    //@{
    template <typename TT>
    Matrix<R, C, T>& rotate_x(
        TT angle, typename std::enable_if<R == C && R == 4, TT>::type* = 0);

    template <typename TT>
    Matrix<R, C, T>& rotate_y(
        TT angle, typename std::enable_if<R == C && R == 4, TT>::type* = 0);

    template <typename TT>
    Matrix<R, C, T>& rotate_z(
        TT angle, typename std::enable_if<R == C && R == 4, TT>::type* = 0);

    template <typename TT>
    Matrix<R, C, T>& pre_rotate_x(
        TT angle, typename std::enable_if<R == C && R == 4, TT>::type* = 0);

    template <typename TT>
    Matrix<R, C, T>& pre_rotate_y(
        TT angle, typename std::enable_if<R == C && R == 4, TT>::type* = 0);

    template <typename TT>
    Matrix<R, C, T>& pre_rotate_z(
        TT angle, typename std::enable_if<R == C && R == 4, TT>::type* = 0);

    template <typename TT>
    Matrix<R, C, T>& scale(
        const vector<3, TT>& scale_,
        typename std::enable_if<R == C && R == 4, TT>::type* = 0);

    template <typename TT>
    Matrix<R, C, T>& scaleTranslation(
        const vector<3, TT>& scale_,
        typename std::enable_if<R == C && R == 4, TT>::type* = 0);
    //@}

    friend std::ostream& operator<<(std::ostream& os,
                                    const Matrix<R, C, T>& matrix)
    {
        const std::ios::fmtflags flags = os.flags();
        const int prec = os.precision();

        os.setf(std::ios::right, std::ios::adjustfield);
        os.precision(7);

        for (size_t rowIndex = 0; rowIndex < R; ++rowIndex)
        {
            os << "|";
            for (size_t colIndex = 0; colIndex < C; ++colIndex)
                os << std::setw(10) << matrix(rowIndex, colIndex) << " ";
            os << "|" << std::endl;
        }
        os.precision(prec);
        os.setf(flags);
        return os;
    };

    T array[R * C]; //!< column by column storage
};
}

#include <vmmlib/quaternion.hpp>
#include <vmmlib/vector.hpp>

namespace vmml
{
/** @name Free matrix functions */
//@{
template <typename T>
inline T computeDeterminant(const Matrix<1, 1, T>& matrix_)
{
    return matrix_.array[0];
}

template <typename T>
inline T computeDeterminant(const Matrix<2, 2, T>& matrix_)
{
    return matrix_(0, 0) * matrix_(1, 1) - matrix_(0, 1) * matrix_(1, 0);
}

template <typename T>
inline T computeDeterminant(const Matrix<3, 3, T>& m_)
{
    return m_(0, 0) * (m_(1, 1) * m_(2, 2) - m_(1, 2) * m_(2, 1)) +
           m_(0, 1) * (m_(1, 2) * m_(2, 0) - m_(1, 0) * m_(2, 2)) +
           m_(0, 2) * (m_(1, 0) * m_(2, 1) - m_(1, 1) * m_(2, 0));
}

template <typename T>
inline T computeDeterminant(const Matrix<4, 4, T>& m)
{
    return m(0, 3) * m(1, 2) * m(2, 1) * m(3, 0) -
           m(0, 2) * m(1, 3) * m(2, 1) * m(3, 0) -
           m(0, 3) * m(1, 1) * m(2, 2) * m(3, 0) +
           m(0, 1) * m(1, 3) * m(2, 2) * m(3, 0) +
           m(0, 2) * m(1, 1) * m(2, 3) * m(3, 0) -
           m(0, 1) * m(1, 2) * m(2, 3) * m(3, 0) -
           m(0, 3) * m(1, 2) * m(2, 0) * m(3, 1) +
           m(0, 2) * m(1, 3) * m(2, 0) * m(3, 1) +
           m(0, 3) * m(1, 0) * m(2, 2) * m(3, 1) -
           m(0, 0) * m(1, 3) * m(2, 2) * m(3, 1) -
           m(0, 2) * m(1, 0) * m(2, 3) * m(3, 1) +
           m(0, 0) * m(1, 2) * m(2, 3) * m(3, 1) +
           m(0, 3) * m(1, 1) * m(2, 0) * m(3, 2) -
           m(0, 1) * m(1, 3) * m(2, 0) * m(3, 2) -
           m(0, 3) * m(1, 0) * m(2, 1) * m(3, 2) +
           m(0, 0) * m(1, 3) * m(2, 1) * m(3, 2) +
           m(0, 1) * m(1, 0) * m(2, 3) * m(3, 2) -
           m(0, 0) * m(1, 1) * m(2, 3) * m(3, 2) -
           m(0, 2) * m(1, 1) * m(2, 0) * m(3, 3) +
           m(0, 1) * m(1, 2) * m(2, 0) * m(3, 3) +
           m(0, 2) * m(1, 0) * m(2, 1) * m(3, 3) -
           m(0, 0) * m(1, 2) * m(2, 1) * m(3, 3) -
           m(0, 1) * m(1, 0) * m(2, 2) * m(3, 3) +
           m(0, 0) * m(1, 1) * m(2, 2) * m(3, 3);
}

template <typename T>
Matrix<1, 1, T> computeInverse(const Matrix<1, 1, T>& m_)
{
    return Matrix<1, 1, T>(std::vector<T>(T(1) / m_(0, 0), 1));
}

template <typename T>
Matrix<2, 2, T> computeInverse(const Matrix<2, 2, T>& m_)
{
    const T det = computeDeterminant(m_);
    if (std::abs(det) < std::numeric_limits<T>::epsilon())
        return Matrix<2, 2, T>(
            std::vector<T>(4, std::numeric_limits<T>::quiet_NaN()));

    Matrix<2, 2, T> inverse;
    m_.getAdjugate(inverse);
    const T detinv = 1 / det;
    inverse(0, 0) *= detinv;
    inverse(0, 1) *= detinv;
    inverse(1, 0) *= detinv;
    inverse(1, 1) *= detinv;
    return inverse;
}

template <typename T>
Matrix<3, 3, T> computeInverse(const Matrix<3, 3, T>& m_)
{
    // Invert a 3x3 using cofactors.  This is about 8 times faster than
    // the Numerical Recipes code which uses Gaussian elimination.
    Matrix<3, 3, T> inverse;
    inverse(0, 0) = m_(1, 1) * m_(2, 2) - m_(1, 2) * m_(2, 1);
    inverse(0, 1) = m_(0, 2) * m_(2, 1) - m_(0, 1) * m_(2, 2);
    inverse(0, 2) = m_(0, 1) * m_(1, 2) - m_(0, 2) * m_(1, 1);
    inverse(1, 0) = m_(1, 2) * m_(2, 0) - m_(1, 0) * m_(2, 2);
    inverse(1, 1) = m_(0, 0) * m_(2, 2) - m_(0, 2) * m_(2, 0);
    inverse(1, 2) = m_(0, 2) * m_(1, 0) - m_(0, 0) * m_(1, 2);
    inverse(2, 0) = m_(1, 0) * m_(2, 1) - m_(1, 1) * m_(2, 0);
    inverse(2, 1) = m_(0, 1) * m_(2, 0) - m_(0, 0) * m_(2, 1);
    inverse(2, 2) = m_(0, 0) * m_(1, 1) - m_(0, 1) * m_(1, 0);

    const T determinant = m_(0, 0) * inverse(0, 0) + m_(0, 1) * inverse(1, 0) +
                          m_(0, 2) * inverse(2, 0);

    if (std::abs(determinant) <= std::numeric_limits<T>::epsilon())
        return Matrix<3, 3, T>(
            std::vector<T>(9, std::numeric_limits<T>::quiet_NaN()));

    const T detinv = static_cast<T>(1.0) / determinant;

    inverse(0, 0) *= detinv;
    inverse(0, 1) *= detinv;
    inverse(0, 2) *= detinv;
    inverse(1, 0) *= detinv;
    inverse(1, 1) *= detinv;
    inverse(1, 2) *= detinv;
    inverse(2, 0) *= detinv;
    inverse(2, 1) *= detinv;
    inverse(2, 2) *= detinv;
    return inverse;
}

template <typename T>
Matrix<4, 4, T> computeInverse(const Matrix<4, 4, T>& m_)
{
    const T* array = m_.array;

    // tuned version from Claude Knaus
    /* first set of 2x2 determinants: 12 multiplications, 6 additions */
    const T t1[6] = {array[2] * array[7] - array[6] * array[3],
                     array[2] * array[11] - array[10] * array[3],
                     array[2] * array[15] - array[14] * array[3],
                     array[6] * array[11] - array[10] * array[7],
                     array[6] * array[15] - array[14] * array[7],
                     array[10] * array[15] - array[14] * array[11]};

    /* first half of comatrix: 24 multiplications, 16 additions */
    Matrix<4, 4, T> inv;
    inv.array[0] = array[5] * t1[5] - array[9] * t1[4] + array[13] * t1[3];
    inv.array[1] = array[9] * t1[2] - array[13] * t1[1] - array[1] * t1[5];
    inv.array[2] = array[13] * t1[0] - array[5] * t1[2] + array[1] * t1[4];
    inv.array[3] = array[5] * t1[1] - array[1] * t1[3] - array[9] * t1[0];
    inv.array[4] = array[8] * t1[4] - array[4] * t1[5] - array[12] * t1[3];
    inv.array[5] = array[0] * t1[5] - array[8] * t1[2] + array[12] * t1[1];
    inv.array[6] = array[4] * t1[2] - array[12] * t1[0] - array[0] * t1[4];
    inv.array[7] = array[0] * t1[3] - array[4] * t1[1] + array[8] * t1[0];

    /* second set of 2x2 determinants: 12 multiplications, 6 additions */
    const T t2[6] = {array[0] * array[5] - array[4] * array[1],
                     array[0] * array[9] - array[8] * array[1],
                     array[0] * array[13] - array[12] * array[1],
                     array[4] * array[9] - array[8] * array[5],
                     array[4] * array[13] - array[12] * array[5],
                     array[8] * array[13] - array[12] * array[9]};

    /* second half of comatrix: 24 multiplications, 16 additions */
    inv.array[8] = array[7] * t2[5] - array[11] * t2[4] + array[15] * t2[3];
    inv.array[9] = array[11] * t2[2] - array[15] * t2[1] - array[3] * t2[5];
    inv.array[10] = array[15] * t2[0] - array[7] * t2[2] + array[3] * t2[4];
    inv.array[11] = array[7] * t2[1] - array[3] * t2[3] - array[11] * t2[0];
    inv.array[12] = array[10] * t2[4] - array[6] * t2[5] - array[14] * t2[3];
    inv.array[13] = array[2] * t2[5] - array[10] * t2[2] + array[14] * t2[1];
    inv.array[14] = array[6] * t2[2] - array[14] * t2[0] - array[2] * t2[4];
    inv.array[15] = array[2] * t2[3] - array[6] * t2[1] + array[10] * t2[0];

    /* determinant: 4 multiplications, 3 additions */
    const T determinant = array[0] * inv.array[0] + array[4] * inv.array[1] +
                          array[8] * inv.array[2] + array[12] * inv.array[3];

    /* division: 16 multiplications, 1 division */
    const T detinv = T(1) / determinant;
    for (size_t i = 0; i != 16; ++i)
        inv.array[i] *= detinv;
    return inv;
}

template <size_t R, size_t C, typename T>
Matrix<R, C, T> computeInverse(const Matrix<R, C, T>&)
{
    throw std::runtime_error("Can't compute inverse of this matrix");
}

/** @return the transposed of a matrix */
template <size_t R, size_t C, typename T>
inline Matrix<C, R, T> transpose(const Matrix<R, C, T>& matrix)
{
    Matrix<C, R, T> transposed;
    for (size_t row = 0; row < R; ++row)
        for (size_t col = 0; col < C; ++col)
            transposed(col, row) = matrix(row, col);
    return transposed;
}
//@}

template <size_t R, size_t C, typename T>
Matrix<R, C, T>::Matrix()
    : array() // http://stackoverflow.com/questions/5602030
{
    if (R == C)
        for (size_t i = 0; i < R; ++i)
            (*this)(i, i) = 1;
}

template <size_t R, size_t C, typename T>
Matrix<R, C, T>::Matrix(const T* begin_, const T* end_)
    : array() // http://stackoverflow.com/questions/5602030
{
    const T* to = std::min(end_, begin_ + R * C);
    ::memcpy(array, begin_, (to - begin_) * sizeof(T));
}

template <size_t R, size_t C, typename T>
Matrix<R, C, T>::Matrix(const std::vector<T>& values)
{
    *this = values;
}

template <size_t R, size_t C, typename T>
Matrix<R, C, T>::Matrix(const Matrix<R, C, T>& source)
{
    *this = source;
}

template <size_t R, size_t C, typename T>
template <size_t P, size_t Q>
Matrix<R, C, T>::Matrix(const Matrix<P, Q, T>& source)
{
    *this = source;
}

template <size_t R, size_t C, typename T>
template <size_t O>
Matrix<R, C, T>::Matrix(
    const Quaternion<T>& rotation, const vector<O, T>& translation,
    typename std::enable_if<R == O + 1 && C == O + 1 && O == 3>::type*)
{
    *this = rotation.getRotationMatrix();
    setTranslation(translation);
    (*this)(3, 0) = 0;
    (*this)(3, 1) = 0;
    (*this)(3, 2) = 0;
    (*this)(3, 3) = 1;
}

template <size_t R, size_t C, typename T>
template <size_t S>
Matrix<R, C, T>::Matrix(
    const vector<S, T>& eye, const vector<S, T>& lookat, const vector<S, T>& up,
    typename std::enable_if<R == S + 1 && C == S + 1 && S == 3>::type*)
    : array() // http://stackoverflow.com/questions/5602030
{
    const vector<3, T> f(vmml::normalize(lookat - eye));
    const vector<3, T> s(vmml::normalize(vmml::cross(f, up)));
    const vector<3, T> u(vmml::cross(s, f));

    (*this)(0, 0) = s.x();
    (*this)(0, 1) = s.y();
    (*this)(0, 2) = s.z();
    (*this)(1, 0) = u.x();
    (*this)(1, 1) = u.y();
    (*this)(1, 2) = u.z();
    (*this)(2, 0) = -f.x();
    (*this)(2, 1) = -f.y();
    (*this)(2, 2) = -f.z();
    (*this)(0, 3) = -vmml::dot(s, eye);
    (*this)(1, 3) = -vmml::dot(u, eye);
    (*this)(2, 3) = vmml::dot(f, eye);
    (*this)(3, 3) = 1;
}

template <size_t R, size_t C, typename T>
inline T& Matrix<R, C, T>::operator()(size_t rowIndex, size_t colIndex)
{
    if (rowIndex >= R || colIndex >= C)
        throw std::runtime_error("matrix( row, column ) index out of bounds");
    return array[colIndex * R + rowIndex];
}

template <size_t R, size_t C, typename T>
inline T Matrix<R, C, T>::operator()(size_t rowIndex, size_t colIndex) const
{
    if (rowIndex >= R || colIndex >= C)
        throw std::runtime_error("matrix( row, column ) index out of bounds");
    return array[colIndex * R + rowIndex];
}

template <size_t R, size_t C, typename T>
bool Matrix<R, C, T>::operator==(const Matrix<R, C, T>& other) const
{
    for (size_t i = 0; i < R * C; ++i)
        if (array[i] != other.array[i])
            return false;
    return true;
}

template <size_t R, size_t C, typename T>
bool Matrix<R, C, T>::operator!=(const Matrix<R, C, T>& other) const
{
    return !operator==(other);
}

template <size_t R, size_t C, typename T>
bool Matrix<R, C, T>::equals(const Matrix<R, C, T>& other,
                             const T tolerance) const
{
    for (size_t i = 0; i < R * C; ++i)
        if (std::abs(array[i] - other.array[i]) > tolerance)
            return false;
    return true;
}

template <size_t R, size_t C, typename T>
const Matrix<R, C, T>& Matrix<R, C, T>::operator=(const Matrix<R, C, T>& source)
{
    ::memcpy(array, source.array, R * C * sizeof(T));
    return *this;
}

template <size_t R, size_t C, typename T>
template <size_t P, size_t Q>
const Matrix<R, C, T>& Matrix<R, C, T>::operator=(const Matrix<P, Q, T>& source)
{
    const size_t minL = std::min(P, R);
    const size_t minC = std::min(Q, C);

    for (size_t i = 0; i < minL; ++i)
        for (size_t j = 0; j < minC; ++j)
            (*this)(i, j) = source(i, j);
    for (size_t i = minL; i < R; ++i)
        for (size_t j = minC; j < C; ++j)
            (*this)(i, j) = 0;
    return *this;
}

template <size_t R, size_t C, typename T>
void Matrix<R, C, T>::operator=(const std::vector<T>& values)
{
    const size_t to = std::min(R * C, values.size());
    ::memcpy(array, values.data(), to * sizeof(T));
    if (to < R * C)
    {
#ifdef _WIN32
        ::memset(array + to, 0, (R * C - to) * sizeof(T));
#else
        ::bzero(array + to, (R * C - to) * sizeof(T));
#endif
    }
}

template <size_t R, size_t C, typename T>
template <size_t P>
Matrix<R, C, T>& Matrix<R, C, T>::multiply(const Matrix<R, P, T>& left,
                                           const Matrix<P, C, T>& right)
{
    // Create copy for multiplication with self
    if (&left == this)
        return multiply(Matrix<R, P, T>(left), right);
    if (&right == this)
        return multiply(left, Matrix<R, P, T>(right));

    for (size_t rowIndex = 0; rowIndex < R; ++rowIndex)
    {
        for (size_t colIndex = 0; colIndex < C; ++colIndex)
        {
            T& component = (*this)(rowIndex, colIndex);
            component = static_cast<T>(0.0);
            for (size_t p = 0; p < P; p++)
                component += left(rowIndex, p) * right(p, colIndex);
        }
    }
    return *this;
}

template <size_t R, size_t C, typename T>
template <size_t P>
Matrix<R, P, T> Matrix<R, C, T>::operator*(const Matrix<C, P, T>& other) const
{
    Matrix<R, P, T> result;
    return result.multiply(*this, other);
}

template <size_t R, size_t C, typename T>
template <size_t O, size_t P, typename>
Matrix<R, C, T>& Matrix<R, C, T>::operator*=(const Matrix<O, P, T>& right)
{
    Matrix<R, C, T> copy(*this);
    multiply(copy, right);
    return *this;
}

template <size_t R, size_t C, typename T>
vector<R, T> Matrix<R, C, T>::operator*(const vector<C, T>& vec) const
{
    vector<R, T> result;

    // this < R, 1 > = < R, P > * < P, 1 >
    for (size_t i = 0; i < R; ++i)
    {
        T tmp(0);
        for (size_t j = 0; j < C; ++j)
            tmp += (*this)(i, j) * vec(j);
        result(i) = tmp;
    }
    return result;
}

template <size_t R, size_t C, typename T>
inline Matrix<R, C, T> Matrix<R, C, T>::operator-() const
{
    Matrix<R, C, T> result;
    for (size_t i = 0; i < R * C; ++i)
        result.array[i] = -array[i];
    return result;
}

template <size_t R, size_t C, typename T>
vector<R, T> Matrix<R, C, T>::getColumn(const size_t index) const
{
    if (index >= C)
        throw std::runtime_error("getColumn() - index out of bounds.");
    vector<R, T> column;
    ::memcpy(&column.array[0], &array[R * index], R * sizeof(T));
    return column;
}

template <size_t R, size_t C, typename T>
void Matrix<R, C, T>::setColumn(size_t index, const vector<R, T>& column)
{
    if (index >= C)
        throw std::runtime_error("setColumn() - index out of bounds.");
    memcpy(array + R * index, column.array, R * sizeof(T));
}

template <size_t R, size_t C, typename T>
vector<C, T> Matrix<R, C, T>::getRow(size_t index) const
{
    if (index >= R)
        throw std::runtime_error("getRow() - index out of bounds.");

    vector<C, T> row;
    for (size_t colIndex = 0; colIndex < C; ++colIndex)
        row(colIndex) = (*this)(index, colIndex);
    return row;
}

template <size_t R, size_t C, typename T>
void Matrix<R, C, T>::setRow(size_t rowIndex, const vector<C, T>& row)
{
    if (rowIndex >= R)
        throw std::runtime_error("setRow() - index out of bounds.");

    for (size_t colIndex = 0; colIndex < C; ++colIndex)
        (*this)(rowIndex, colIndex) = row(colIndex);
}

template <size_t R, size_t C, typename T>
inline Matrix<R, C, T> Matrix<R, C, T>::operator+(
    const Matrix<R, C, T>& other) const
{
    Matrix<R, C, T> result(*this);
    result += other;
    return result;
}

template <size_t R, size_t C, typename T>
void Matrix<R, C, T>::operator+=(const Matrix<R, C, T>& other)
{
    for (size_t i = 0; i < R * C; ++i)
        array[i] += other.array[i];
}

template <size_t R, size_t C, typename T>
inline Matrix<R, C, T> Matrix<R, C, T>::operator-(
    const Matrix<R, C, T>& other) const
{
    Matrix<R, C, T> result(*this);
    result -= other;
    return result;
}

template <size_t R, size_t C, typename T>
void Matrix<R, C, T>::operator-=(const Matrix<R, C, T>& other)
{
    for (size_t i = 0; i < R * C; ++i)
        array[i] -= other.array[i];
}

template <size_t R, size_t C, typename T>
template <size_t O, size_t P>
Matrix<O, P, T> Matrix<R, C, T>::getSubMatrix(
    size_t rowOffset, size_t colOffset,
    typename std::enable_if<O <= R && P <= C>::type*) const
{
    Matrix<O, P, T> result;
    if (O + rowOffset > R || P + colOffset > C)
        throw std::runtime_error("index out of bounds.");

    for (size_t row = 0; row < O; ++row)
        for (size_t col = 0; col < P; ++col)
            result(row, col) = (*this)(rowOffset + row, colOffset + col);

    return result;
}

template <size_t R, size_t C, typename T>
template <size_t O, size_t P>
typename std::enable_if<O <= R && P <= C>::type* Matrix<R, C, T>::setSubMatrix(
    const Matrix<O, P, T>& sub_matrix, size_t rowOffset, size_t colOffset)
{
    for (size_t row = 0; row < O; ++row)
        for (size_t col = 0; col < P; ++col)
            (*this)(rowOffset + row, colOffset + col) = sub_matrix(row, col);
    return 0; // for sfinae
}

template <size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::inverse() const
{
    return computeInverse(*this);
}

template <size_t R, size_t C, typename T>
template <size_t O, size_t P>
typename std::enable_if<O == P && R == C && O == R && R >= 2>::type*
    Matrix<R, C, T>::getAdjugate(Matrix<O, P, T>& adjugate) const
{
    getCofactors(adjugate);
    adjugate = transpose(adjugate);
    return 0;
}

template <size_t R, size_t C, typename T>
template <size_t O, size_t P>
typename std::enable_if<O == P && R == C && O == R && R >= 2>::type*
    Matrix<R, C, T>::getCofactors(Matrix<O, P, T>& cofactors) const
{
    Matrix<R - 1, C - 1, T> minor_;

    const size_t _negate = 1u;
    for (size_t rowIndex = 0; rowIndex < R; ++rowIndex)
        for (size_t colIndex = 0; colIndex < C; ++colIndex)
            if ((rowIndex + colIndex) & _negate)
                cofactors(rowIndex, colIndex) =
                    -getMinor(minor_, rowIndex, colIndex);
            else
                cofactors(rowIndex, colIndex) =
                    getMinor(minor_, rowIndex, colIndex);

    return 0;
}

template <size_t R, size_t C, typename T>
template <size_t O, size_t P>
T Matrix<R, C, T>::getMinor(
    Matrix<O, P, T>& minor_, size_t row_to_cut, size_t col_to_cut,
    typename std::enable_if<O == R - 1 && P == C - 1 && R == C &&
                            R >= 2>::type*) const
{
    size_t rowOffset = 0;
    size_t colOffset = 0;
    for (size_t rowIndex = 0; rowIndex < R; ++rowIndex)
    {
        if (rowIndex == row_to_cut)
            rowOffset = -1;
        else
        {
            for (size_t colIndex = 0; colIndex < R; ++colIndex)
            {
                if (colIndex == col_to_cut)
                    colOffset = -1;
                else
                    minor_(rowIndex + rowOffset, colIndex + colOffset) =
                        (*this)(rowIndex, colIndex);
            }
            colOffset = 0;
        }
    }
    return computeDeterminant(minor_);
}

template <size_t R, size_t C, typename T>
template <typename TT>
Matrix<R, C, T>& Matrix<R, C, T>::rotate_x(
    const TT angle_, typename std::enable_if<R == C && R == 4, TT>::type*)
{
    const T angle = static_cast<T>(angle_);
    const T sine = std::sin(angle);
    const T cosine = std::cos(angle);

    T tmp;

    tmp = array[4] * cosine + array[8] * sine;
    array[8] = -array[4] * sine + array[8] * cosine;
    array[4] = tmp;

    tmp = array[5] * cosine + array[9] * sine;
    array[9] = -array[5] * sine + array[9] * cosine;
    array[5] = tmp;

    tmp = array[6] * cosine + array[10] * sine;
    array[10] = -array[6] * sine + array[10] * cosine;
    array[6] = tmp;

    tmp = array[7] * cosine + array[11] * sine;
    array[11] = -array[7] * sine + array[11] * cosine;
    array[7] = tmp;

    return *this;
}

template <size_t R, size_t C, typename T>
template <typename TT>
Matrix<R, C, T>& Matrix<R, C, T>::rotate_y(
    const TT angle_, typename std::enable_if<R == C && R == 4, TT>::type*)
{
    const T angle = static_cast<T>(angle_);
    const T sine = std::sin(angle);
    const T cosine = std::cos(angle);

    T tmp;

    tmp = array[0] * cosine - array[8] * sine;
    array[8] = array[0] * sine + array[8] * cosine;
    array[0] = tmp;

    tmp = array[1] * cosine - array[9] * sine;
    array[9] = array[1] * sine + array[9] * cosine;
    array[1] = tmp;

    tmp = array[2] * cosine - array[10] * sine;
    array[10] = array[2] * sine + array[10] * cosine;
    array[2] = tmp;

    tmp = array[3] * cosine - array[11] * sine;
    array[11] = array[3] * sine + array[11] * cosine;
    array[3] = tmp;

    return *this;
}

template <size_t R, size_t C, typename T>
template <typename TT>
Matrix<R, C, T>& Matrix<R, C, T>::rotate_z(
    const TT angle_, typename std::enable_if<R == C && R == 4, TT>::type*)
{
    const T angle = static_cast<T>(angle_);
    const T sine = std::sin(angle);
    const T cosine = std::cos(angle);

    T tmp;

    tmp = array[0] * cosine + array[4] * sine;
    array[4] = -array[0] * sine + array[4] * cosine;
    array[0] = tmp;

    tmp = array[1] * cosine + array[5] * sine;
    array[5] = -array[1] * sine + array[5] * cosine;
    array[1] = tmp;

    tmp = array[2] * cosine + array[6] * sine;
    array[6] = -array[2] * sine + array[6] * cosine;
    array[2] = tmp;

    tmp = array[3] * cosine + array[7] * sine;
    array[7] = -array[3] * sine + array[7] * cosine;
    array[3] = tmp;

    return *this;
}

template <size_t R, size_t C, typename T>
template <typename TT>
Matrix<R, C, T>& Matrix<R, C, T>::pre_rotate_x(
    const TT angle_, typename std::enable_if<R == C && R == 4, TT>::type*)
{
    const T angle = static_cast<T>(angle_);
    const T sine = std::sin(angle);
    const T cosine = std::cos(angle);

    T tmp;

    tmp = array[1];
    array[1] = array[1] * cosine + array[2] * sine;
    array[2] = tmp * -sine + array[2] * cosine;

    tmp = array[5];
    array[5] = array[5] * cosine + array[6] * sine;
    array[6] = tmp * -sine + array[6] * cosine;

    tmp = array[9];
    array[9] = array[9] * cosine + array[10] * sine;
    array[10] = tmp * -sine + array[10] * cosine;

    tmp = array[13];
    array[13] = array[13] * cosine + array[14] * sine;
    array[14] = tmp * -sine + array[14] * cosine;

    return *this;
}

template <size_t R, size_t C, typename T>
template <typename TT>
Matrix<R, C, T>& Matrix<R, C, T>::pre_rotate_y(
    const TT angle_, typename std::enable_if<R == C && R == 4, TT>::type*)
{
    const T angle = static_cast<T>(angle_);
    const T sine = std::sin(angle);
    const T cosine = std::cos(angle);

    T tmp;

    tmp = array[0];
    array[0] = array[0] * cosine - array[2] * sine;
    array[2] = tmp * sine + array[2] * cosine;

    tmp = array[4];
    array[4] = array[4] * cosine - array[6] * sine;
    array[6] = tmp * sine + array[6] * cosine;

    tmp = array[8];
    array[8] = array[8] * cosine - array[10] * sine;
    array[10] = tmp * sine + array[10] * cosine;

    tmp = array[12];
    array[12] = array[12] * cosine - array[14] * sine;
    array[14] = tmp * sine + array[14] * cosine;

    return *this;
}

template <size_t R, size_t C, typename T>
template <typename TT>
Matrix<R, C, T>& Matrix<R, C, T>::pre_rotate_z(
    const TT angle_, typename std::enable_if<R == C && R == 4, TT>::type*)
{
    const T angle = static_cast<T>(angle_);
    const T sine = std::sin(angle);
    const T cosine = std::cos(angle);

    T tmp;

    tmp = array[0];
    array[0] = array[0] * cosine + array[1] * sine;
    array[1] = tmp * -sine + array[1] * cosine;

    tmp = array[4];
    array[4] = array[4] * cosine + array[5] * sine;
    array[5] = tmp * -sine + array[5] * cosine;

    tmp = array[8];
    array[8] = array[8] * cosine + array[9] * sine;
    array[9] = tmp * -sine + array[9] * cosine;

    tmp = array[12];
    array[12] = array[12] * cosine + array[13] * sine;
    array[13] = tmp * -sine + array[13] * cosine;

    return *this;
}

template <size_t R, size_t C, typename T>
template <typename TT>
inline Matrix<R, C, T>& Matrix<R, C, T>::scale(
    const vector<3, TT>& scale_,
    typename std::enable_if<R == C && R == 4, TT>::type*)
{
    array[0] *= scale_[0];
    array[1] *= scale_[0];
    array[2] *= scale_[0];
    array[3] *= scale_[0];
    array[4] *= scale_[1];
    array[5] *= scale_[1];
    array[6] *= scale_[1];
    array[7] *= scale_[1];
    array[8] *= scale_[2];
    array[9] *= scale_[2];
    array[10] *= scale_[2];
    array[11] *= scale_[2];
    return *this;
}

template <size_t R, size_t C, typename T>
template <typename TT>
inline Matrix<R, C, T>& Matrix<R, C, T>::scaleTranslation(
    const vector<3, TT>& scale_,
    typename std::enable_if<R == C && R == 4, TT>::type*)
{
    array[12] *= static_cast<T>(scale_[0]);
    array[13] *= static_cast<T>(scale_[1]);
    array[14] *= static_cast<T>(scale_[2]);
    return *this;
}

template <size_t R, size_t C, typename T>
inline Matrix<R, C, T>& Matrix<R, C, T>::setTranslation(
    const vector<C - 1, T>& trans)
{
    for (size_t i = 0; i < C - 1; ++i)
        array[i + R * (C - 1)] = trans[i];
    return *this;
}

template <size_t R, size_t C, typename T>
inline vector<C - 1, T> Matrix<R, C, T>::getTranslation() const
{
    vector<C - 1, T> result;
    for (size_t i = 0; i < C - 1; ++i)
        result[i] = array[i + R * (C - 1)];
    return result;
}

template <size_t R, size_t C, typename T>
template <size_t S>
void Matrix<R, C, T>::getLookAt(
    vector<S, T>& eye, vector<S, T>& lookAt, vector<S, T>& up,
    typename std::enable_if<R == S + 1 && C == S + 1 && S == 3>::type*) const
{
    const Matrix<4, 4, T>& inv = inverse();
    const Matrix<3, 3, T> rotation(transpose(Matrix<3, 3, T>(*this)));

    eye = vector<S, T>(inv * vector<4, T>::zero());
    up = rotation * vector<3, T>::up();
    lookAt = rotation * vector<3, T>::forward();
    lookAt.normalize();
    lookAt = eye + lookAt;
}

} // namespace vmml

#endif
