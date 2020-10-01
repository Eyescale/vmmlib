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

#ifndef __VMML__QUATERNION__HPP__
#define __VMML__QUATERNION__HPP__

#include <vmmlib/types.hpp>
#include <vmmlib/vector.hpp> // used inline

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

namespace vmml
{
template <typename T>
class Quaternion
{
public:
    /** Construct an identity quaternion */
    Quaternion()
        : array()
    {
        array[3] = 1.;
    }
    Quaternion(T x, T y, T z, T w);

    /** Copy constructor for quaternion */
    Quaternion(const Quaternion& other);

    /** Construct a rotation quaternion */
    Quaternion(T angle, vector<3, T> axis);

    // uses the top-left 3x3 part of the supplied matrix as rotation matrix
    template <size_t M>
    Quaternion(const Matrix<M, M, T>& rotation_matrix_,
               typename std::enable_if<M >= 3>::type* = 0);

    /** @name Data Access */
    //@{
    /** @return true if the two quaternion are similar. */
    bool equals(const Quaternion& other,
                T tolerance = std::numeric_limits<T>::epsilon()) const;

    T x() const { return array[0]; }
    T y() const { return array[1]; }
    T z() const { return array[2]; }
    T w() const { return array[3]; }
    /** @return true if both quaternions are equal. */
    bool operator==(const Quaternion& a) const;

    /** @return true if both quaternions are not equal. */
    bool operator!=(const Quaternion& a) const;

    /** @return the negated quaternion of this quaternion. */
    Quaternion operator-() const;

    /** @return multiplicative inverse quaternion */
    Quaternion inverse() const;

    Quaternion getConjugate() const;

    inline vector<3,T> rotate(const vector<3,T>& vec) const;

    T abs() const;
    T absSquare() const;

    /** @return the corresponding 3x3 rotation matrix. */
    Matrix<3, 3, T> getRotationMatrix() const;
    //@}

    void normalize();

    Quaternion& operator=(const Quaternion& other);

    /** @name quaternion/quaternion operations */
    //@{
    Quaternion operator+(const Quaternion<T>& a) const;
    Quaternion operator-(const Quaternion<T>& a) const;

    // caution: a * q != q * a in general
    Quaternion operator*(const Quaternion<T>& a) const;
    void operator+=(const Quaternion<T>& a);
    void operator-=(const Quaternion<T>& a);

    // caution: a *= q != q *= a in general
    void operator*=(const Quaternion<T>& a);
    //@}

    /** @name quaternion/scalar operations */
    //@{
    Quaternion operator*(T a) const;
    Quaternion operator/(T a) const;

    void operator*=(T a);
    void operator/=(T a);
    //@}

    friend std::ostream& operator<<(std::ostream& os, const Quaternion& q)
    {
        const std::ios::fmtflags flags = os.flags();
        const int prec = os.precision();

        os.setf(std::ios::right, std::ios::adjustfield);
        os.precision(5);

        os << "( " << std::setw(10) << q.x() << " " << std::setw(10) << q.y()
           << " " << std::setw(10) << q.z() << " " << std::setw(10) << q.w()
           << " "
           << ")";

        os.precision(prec);
        os.setf(flags);
        return os;
    }

    T array[4];
};
}

#include <vmmlib/matrix.hpp>

namespace vmml
{
/** @name Free quaternion functions */
//@{
template <typename T>
T dot(const Quaternion<T>& p, const Quaternion<T>& q)
{
    return p.array[3] * q.array[3] + p.array[0] * q.array[0] +
           p.array[1] * q.array[1] + p.array[2] * q.array[2];
}

template <typename T>
vector<3, T> cross(const Quaternion<T>& p, const Quaternion<T>& q)
{
    return vector<3, T>(p.array[1] * q.array[2] - p.array[2] * q.array[1],
                        p.array[2] * q.array[0] - p.array[0] * q.array[2],
                        p.array[0] * q.array[1] - p.array[1] * q.array[0]);
}
//@}

template <typename T>
Quaternion<T>::Quaternion(T x_, T y_, T z_, T w_)
{
    array[0] = x_;
    array[1] = y_;
    array[2] = z_;
    array[3] = w_;
}

template <typename T>
Quaternion<T>::Quaternion(const Quaternion& other)
{
    ::memcpy(array, other.array, 4 * sizeof(T));
}

template <typename T>
template <size_t M>
Quaternion<T>::Quaternion(const Matrix<M, M, T>& rot,
                          typename std::enable_if<M >= 3>::type*)
{
    const T trace = rot(0, 0) + rot(1, 1) + rot(2, 2) + T(1);
    static const T TRACE_EPSILON = 1e-5;

    // very small traces may introduce a big numerical error
    if (trace > TRACE_EPSILON)
    {
        T s = 0.5 / std::sqrt(trace);
        array[0] = rot(2, 1) - rot(1, 2);
        array[0] *= s;

        array[1] = rot(0, 2) - rot(2, 0);
        array[1] *= s;

        array[2] = rot(1, 0) - rot(0, 1);
        array[2] *= s;

        array[3] = 0.25 / s;
    }
    else
    {
        const vector<3, T> diag(rot(0, 0), rot(1, 1), rot(2, 2));
        const size_t largest = diag.find_max_index();

        // 0, 0 is largest
        if (largest == 0)
        {
            const T s =
                0.5 / std::sqrt(T(1) + rot(0, 0) - rot(1, 1) - rot(2, 2));
            array[0] = 0.25 / s;

            array[1] = rot(0, 1) + rot(1, 0);
            array[1] *= s;

            array[2] = rot(0, 2) + rot(2, 0);
            array[2] *= s;

            array[3] = rot(1, 2) - rot(2, 1);
            array[3] *= s;
        }
        else if (largest == 1)
        {
            const T s =
                0.5 / std::sqrt(T(1) + rot(1, 1) - rot(0, 0) - rot(2, 2));
            array[0] = rot(0, 1) + rot(1, 0);
            array[0] *= s;

            array[1] = 0.25 / s;

            array[2] = rot(1, 2) + rot(2, 1);
            array[2] *= s;

            array[3] = rot(0, 2) - rot(2, 0);
            array[3] *= s;
        }
        // 2, 2 is largest
        else if (largest == 2)
        {
            const T s =
                0.5 / std::sqrt(T(1) + rot(2, 2) - rot(0, 0) - rot(1, 1));
            array[0] = rot(0, 2) + rot(2, 0);
            array[0] *= s;

            array[1] = rot(1, 2) + rot(2, 1);
            array[1] *= s;

            array[2] = 0.25 / s;

            array[3] = rot(0, 1) - rot(1, 0);
            array[3] *= s;
        }
        else
        {
            throw std::runtime_error("no clue why, but should not get here");
        }
    }
}

template <typename T>
Quaternion<T>::Quaternion(const T angle, vector<3, T> axis)
{
    axis.normalize();
    const T sinAngle2 = std::sin(angle * 0.5);
    array[0] = axis.x() * sinAngle2;
    array[1] = axis.y() * sinAngle2;
    array[2] = axis.z() * sinAngle2;
    array[3] = std::cos(angle * 0.5);
}

template <typename T>
bool Quaternion<T>::equals(const Quaternion& other, const T tolerance) const
{
    return std::abs(array[0] - other.array[0]) <= tolerance &&
           std::abs(array[1] - other.array[1]) <= tolerance &&
           std::abs(array[2] - other.array[2]) <= tolerance &&
           std::abs(array[3] - other.array[3]) <= tolerance;
}

template <typename T>
bool Quaternion<T>::operator==(const Quaternion& rhs) const
{
    return array[0] == rhs.array[0] && array[1] == rhs.array[1] &&
           array[2] == rhs.array[2] && array[3] == rhs.array[3];
}

template <typename T>
bool Quaternion<T>::operator!=(const Quaternion& a) const
{
    return !this->operator==(a);
}

template <typename T>
Quaternion<T> Quaternion<T>::getConjugate() const
{
    return Quaternion<T>(-array[0], -array[1], -array[2], array[3]);
}

template <typename T>
inline vector<3,T> Quaternion<T>::rotate(const vector<3,T>& v) const
{
    const vector<3,T> u(array[0], array[1], array[2]);
    const T s = array[3];
    return 2 * vmml::dot(u,v) * u + (s * s - vmml::dot(u,u)) * v + 2 * s * vmml::cross(u,v);
}

template <typename T>
T Quaternion<T>::abs() const
{
    return std::sqrt(absSquare());
}

template <typename T>
T Quaternion<T>::absSquare() const
{
    return array[0] * array[0] + array[1] * array[1] + array[2] * array[2] +
           array[3] * array[3];
}

template <typename T>
Quaternion<T> Quaternion<T>::inverse() const
{
    const Quaternion<T>& q = getConjugate();
    const T tmp = T(1) / absSquare();
    return q * tmp;
}

template <typename T>
void Quaternion<T>::normalize()
{
    T len = abs();
    if (len == T(0))
        return;
    len = T(1) / len;
    this->operator*=(len);
}

//
// Quaternion/Quaternion operations
//

template <typename T>
Quaternion<T> Quaternion<T>::operator+(const Quaternion<T>& rhs) const
{
    return Quaternion(array[0] + rhs.array[0], array[1] + rhs.array[1],
                      array[2] + rhs.array[2], array[3] + rhs.array[3]);
}

template <typename T>
Quaternion<T> Quaternion<T>::operator-(const Quaternion<T>& rhs) const
{
    return Quaternion(array[0] - rhs.array[0], array[1] - rhs.array[1],
                      array[2] - rhs.array[2], array[3] - rhs.array[3]);
}

// returns Grasssmann product
template <typename T>
Quaternion<T> Quaternion<T>::operator*(const Quaternion<T>& rhs) const
{
    Quaternion<T> ret(*this);
    ret *= rhs;
    return ret;
}

// Grassmann product
template <typename T>
void Quaternion<T>::operator*=(const Quaternion<T>& q)
{
    // optimized version, 7 less mul, but 15 more add/subs
    // after Henrik Engstrom, from a gamedev.net article.
    const T& a_ = array[3];
    const T& b_ = array[0];
    const T& c = array[1];
    const T& d = array[2];
    const T& _x = q.array[3];
    const T& _y = q.array[0];
    const T& _z = q.array[1];
    const T& _w = q.array[2];

    const T tmp_00 = (d - c) * (_z - _w);
    const T tmp_01 = (a_ + b_) * (_x + _y);
    const T tmp_02 = (a_ - b_) * (_z + _w);
    const T tmp_03 = (c + d) * (_x - _y);
    const T tmp_04 = (d - b_) * (_y - _z);
    const T tmp_05 = (d + b_) * (_y + _z);
    const T tmp_06 = (a_ + c) * (_x - _w);
    const T tmp_07 = (a_ - c) * (_x + _w);
    const T tmp_08 = tmp_05 + tmp_06 + tmp_07;
    const T tmp_09 = 0.5 * (tmp_04 + tmp_08);

    array[0] = tmp_01 + tmp_09 - tmp_08;
    array[1] = tmp_02 + tmp_09 - tmp_07;
    array[2] = tmp_03 + tmp_09 - tmp_06;
    array[3] = tmp_00 + tmp_09 - tmp_05;
}

template <typename T>
Quaternion<T> Quaternion<T>::operator-() const
{
    return Quaternion(-array[0], -array[1], -array[2], -array[3]);
}

template <typename T>
void Quaternion<T>::operator+=(const Quaternion<T>& q)
{
    array[0] += q.array[0];
    array[1] += q.array[1];
    array[2] += q.array[2];
    array[3] += q.array[3];
}

template <typename T>
void Quaternion<T>::operator-=(const Quaternion<T>& q)
{
    array[0] -= q.array[0];
    array[1] -= q.array[1];
    array[2] -= q.array[2];
    array[3] -= q.array[3];
}

//
// Quaternion/scalar operations
//

template <typename T>
Quaternion<T> Quaternion<T>::operator*(const T a_) const
{
    return Quaternion(array[0] * a_, array[1] * a_, array[2] * a_,
                      array[3] * a_);
}

template <typename T>
Quaternion<T> Quaternion<T>::operator/(T a_) const
{
    if (a_ == T(0))
        throw std::runtime_error("Division by zero.");

    a_ = T(1) / a_;
    return Quaternion(array[0] * a_, array[1] * a_, array[2] * a_,
                      array[3] * a_);
}

template <typename T>
void Quaternion<T>::operator*=(T q)
{
    array[0] *= q;
    array[1] *= q;
    array[2] *= q;
    array[3] *= q;
}

template <typename T>
void Quaternion<T>::operator/=(T q)
{
    if (q == T(0))
        throw std::runtime_error("Division by zero");

    q = T(1) / q;
    this->operator*=(q);
}

template <typename T>
Matrix<3, 3, T> Quaternion<T>::getRotationMatrix() const
{
    const T w2 = array[3] * array[3];
    const T x2 = array[0] * array[0];
    const T y2 = array[1] * array[1];
    const T z2 = array[2] * array[2];
    const T wx = array[3] * array[0];
    const T wy = array[3] * array[1];
    const T wz = array[3] * array[2];
    const T xy = array[0] * array[1];
    const T xz = array[0] * array[2];
    const T yz = array[1] * array[2];

    Matrix<3, 3, T> matrix;
    matrix(0, 0) = w2 + x2 - y2 - z2;
    matrix(0, 1) = 2. * (xy - wz);
    matrix(0, 2) = 2. * (xz + wy);
    matrix(1, 0) = 2. * (xy + wz);
    matrix(1, 1) = w2 - x2 + y2 - z2;
    matrix(1, 2) = 2. * (yz - wx);
    matrix(2, 0) = 2. * (xz - wy);
    matrix(2, 1) = 2. * (yz + wx);
    matrix(2, 2) = w2 - x2 - y2 + z2;
    return matrix;
}

template <typename T>
Quaternion<T>& Quaternion<T>::operator=(const Quaternion& other)
{
    ::memcpy(array, other.array, 4 * sizeof(T));
    return *this;
}
}
#endif
