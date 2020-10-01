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

#ifndef __VMML__VECTOR__HPP__
#define __VMML__VECTOR__HPP__

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <type_traits>
#include <vector>

namespace vmml
{
template <size_t M, typename T>
class vector
{
public:
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    static const size_t DIMENSION = M;

    // constructors
    vector()
        : array()
    {
    }                            // http://stackoverflow.com/questions/5602030
    explicit vector(const T& a); // sets all components to a;
    vector(const vector& other);
    vector(const T& x, const T& y);
    vector(const T& x, const T& y, const T& z);
    vector(const T& x, const T& y, const T& z, const T& w);

#ifndef SWIG
    // initializes the first M-1 values from vector_, the last from last_
    template <typename TT>
    vector(const vector<M - 1, TT>& vector_, T last_,
           typename std::enable_if<M == 4, TT>::type* = nullptr);
#endif

    explicit vector(const T* values);

#ifdef __OSG_MATH
    template <typename OSGVEC3>
    explicit vector(const OSGVEC3& from,
                    typename std::enable_if<M == 3, OSGVEC3>::type* = nullptr);
#endif

    // vec< M > with homogeneous coordinates <-> vec< M-1 > conversion ctor
    // to-homogenous-coordinates ctor
    template <typename TT>
    vector(const vector<3, TT>& source_,
           typename std::enable_if<M == 4, TT>::type* = nullptr);

    // from-homogenous-coordinates vector
    template <typename TT>
    vector(const vector<4, TT>& source_,
           typename std::enable_if<M == 3, TT>::type* = nullptr);

    template <typename U>
    vector(const vector<M, U>& source_);

    // iterators
    inline iterator begin();
    inline iterator end();
    inline const_iterator begin() const;
    inline const_iterator end() const;
    inline reverse_iterator rbegin();
    inline reverse_iterator rend();
    inline const_reverse_iterator rbegin() const;
    inline const_reverse_iterator rend() const;

    // conversion operators
    inline operator T*();
    inline operator const T*() const;
    // accessors
    inline T& operator()(size_t index);
    inline const T& operator()(size_t index) const;

    inline T& at(size_t index);
    inline const T& at(size_t index) const;

    // element accessors for M <= 4;
    inline T& x();
    inline T& y();
    inline T& z();
    inline T& w();
    inline const T& x() const;
    inline const T& y() const;
    inline const T& z() const;
    inline const T& w() const;

    // pixel color element accessors for M<= 4
    inline T& r();
    inline T& g();
    inline T& b();
    inline T& a();
    inline const T& r() const;
    inline const T& g() const;
    inline const T& b() const;
    inline const T& a() const;

    bool operator==(const vector& other) const;
    bool operator!=(const vector& other) const;
    bool equals(const vector& other,
                T tolerance = std::numeric_limits<T>::epsilon()) const;
    bool operator<(const vector& other) const;

    // remember kids: c_arrays are dangerous and evil!
    T operator=(T filler);

    vector& operator=(const vector& other);

    // returns void to avoid 'silent' loss of precision when chaining
    template <typename U>
    void operator=(const vector<M, U>& other);

    vector operator*(const vector& other) const;
    vector operator/(const vector& other) const;
    vector operator+(const vector& other) const;
    vector operator-(const vector& other) const;

    void operator*=(const vector& other);
    void operator/=(const vector& other);
    void operator+=(const vector& other);
    void operator-=(const vector& other);

    vector operator*(const T other) const;
    vector operator/(const T other) const;
    vector operator+(const T other) const;
    vector operator-(const T other) const;

    void operator*=(const T other);
    void operator/=(const T other);
    void operator+=(const T other);
    void operator-=(const T other);

    vector operator-() const;

    const vector& negate();

    void set(T a); // sets all components to a;
    template <size_t N>
    void set(const vector<N, T>& v);

    // sets the first few components to a certain value
    void set(T x, T y);
    void set(T x, T y, T z);
    void set(T x, T y, T z, T w);

    template <typename input_iterator_t>
    void iter_set(input_iterator_t begin_, input_iterator_t end_);

    // compute the cross product of two vectors
    // note: there's also a free function:
    // vector<> cross( const vector<>, const vector<> )
    template <typename TT>
    vector<M, T>& cross(const vector<M, TT>& b,
                        typename std::enable_if<M == 3, TT>::type* = nullptr);

    // compute the dot product of two vectors
    // note: there's also a free function:
    // T dot( const vector<>, const vector<> );
    inline T dot(const vector& other) const;

    // normalize the vector
    // note: there's also a free function:
    // vector<> normalize( const vector<> );
    inline T normalize();

    // sets all vector components to random values
    // remember to set srand( seed );
    // if seed is set to -1, srand( seed ) was set outside set_random
    // otherwise srand( seed ) will be called with the given seed
    void set_random(int seed = -1);

    inline T length() const;
    inline T squared_length() const;

    inline T distance(const vector& other) const;
    inline T squared_distance(const vector& other) const;

    /** @return the product of all elements of this vector */
    T product() const;

    template <typename TT>
    vector<3, T>& rotate(T theta, vector<M, TT> axis,
                         typename std::enable_if<M == 3, TT>::type* = nullptr);

    /** @return the sub vector of the given length at the given offset. */
    template <size_t N, size_t O>
    vector<N, T> get_sub_vector(
        typename std::enable_if<M >= N + O>::type* = nullptr) const;

    /** Set the sub vector of the given length at the given offset. */
    template <size_t N, size_t O>
    void set_sub_vector(const vector<N, T>& sub,
                        typename std::enable_if<M >= N + O>::type* = nullptr);

    // sphere functions - sphere layout: center xyz, radius w
    template <typename TT>
    inline vector<3, T> project_point_onto_sphere(
        const vector<3, TT>& point,
        typename std::enable_if<M == 4, TT>::type* = nullptr) const;

    // returns a negative distance if the point lies in the sphere
    template <typename TT>
    inline T distance_to_sphere(
        const vector<3, TT>& point,
        typename std::enable_if<M == 4, TT>::type* = nullptr) const;

    // plane functions - plane layout; normal xyz, distance w
    template <typename TT>
    inline T distance_to_plane(
        const vector<3, TT>& point,
        typename std::enable_if<M == 4, TT>::type* = nullptr) const;

    template <typename TT>
    inline vector<3, T> project_point_onto_plane(
        const vector<3, TT>& point,
        typename std::enable_if<M == 4, TT>::type* = nullptr) const;

    // returns the index of the minimal resp. maximal value in the vector
    size_t find_min_index() const;
    size_t find_max_index() const;

    // returns minimal resp. maximal value in the vector
    T& find_min();
    T& find_max();
    const T& find_min() const;
    const T& find_max() const;

    void clamp(const T& min = 0.0, const T& max = 1.0);

    inline static size_t size(); // returns M

    bool is_unit_vector() const;

    // perturbs each component by randomly + or - the perturbation parameter
    void perturb(T perturbation = 0.0001);

    void sqrt_elementwise();
    double norm() const; // l2 norm

    // computes the reciprocal value for each component, x = 1/x;
    // WARNING: might result in nans if division by 0!
    void reciprocal();
    // computes the reciprocal value for each component, x = 1/x;
    // checks every component for 0, sets to max value if zero.
    void reciprocal_safe();

    template <typename TT>
    void cast_from(const vector<M, TT>& other);

    size_t nnz() const;

    friend std::ostream& operator<<(std::ostream& os, const vector& vector_)
    {
        os << "[ ";
        for (size_t index = 0; index < M; ++index)
            os << vector_.at(index) << " ";
        return os << "]";
    }

    /** @name Convenience presets for 3 and 4 component vectors */
    //@{
    static vector<M, T> zero();
    static vector<M, T> forward();
    static vector<M, T> backward();
    static vector<M, T> up();
    static vector<M, T> down();
    static vector<M, T> left();
    static vector<M, T> right();
    static vector<M, T> unitX();
    static vector<M, T> unitY();
    static vector<M, T> unitZ();
    //@}

    T array[M]; //!< storage
};

//
//  some free functions for convenience
//

template <size_t M, typename T>
bool equals(const vector<M, T>& a, const vector<M, T>& b)
{
    return a.equals(b);
}

// allows float * vector, not only vector * float
template <size_t M, typename T>
static vector<M, T> operator*(T factor, const vector<M, T>& vector_)
{
    return vector_ * factor;
}

template <size_t M, typename T>
inline T dot(const vector<M, T>& first, const vector<M, T>& second)
{
    return first.dot(second);
}

template <size_t M, typename T>
inline vector<M, T> cross(vector<M, T> a, const vector<M, T>& b)
{
    return a.cross(b);
}

template <size_t M, typename T>
vector<M, T> compute_normal(const vector<M, T>& a, const vector<M, T>& b,
                            const vector<M, T>& c)
{
    // right hand system, CCW triangle
    const vector<M, T> u = b - a;
    const vector<M, T> v = c - a;
    vector<M, T> w = cross(u, v);
    w.normalize();
    return w;
}

template <typename T>
vector<3, T> rotate(vector<3, T> vec, const T theta, const vector<3, T>& axis)
{
    return vec.rotate(theta, axis);
}

template <size_t M, typename T>
inline vector<M, T> normalize(vector<M, T> vector_)
{
    vector_.normalize();
    return vector_;
}

template <typename T>
inline vector<4, T> compute_plane(const vector<3, T>& a, const vector<3, T>& b,
                                  const vector<3, T>& c)
{
    const vector<3, T> cb = b - c;
    const vector<3, T> ba = a - b;

    vector<4, T> plane = vector<4, T>(cross(cb, ba));
    plane.normalize();
    plane.w() = -plane.x() * a.x() - plane.y() * a.y() - plane.z() * a.z();
    return plane;
}

template <size_t M, typename T>
vector<M, T>::vector(const T& _a)
{
    for (iterator it = begin(), it_end = end(); it != it_end; ++it)
    {
        *it = _a;
    }
}

template <size_t M, typename T>
vector<M, T>::vector(const vector<M, T>& other)
{
    memcpy(array, other.array, M * sizeof(T));
}

template <size_t M, typename T>
vector<M, T>::vector(const T& _x, const T& _y)
{
    array[0] = _x;
    array[1] = _y;
}

template <size_t M, typename T>
vector<M, T>::vector(const T& _x, const T& _y, const T& _z)
{
    array[0] = _x;
    array[1] = _y;
    array[2] = _z;
}

template <size_t M, typename T>
vector<M, T>::vector(const T& _x, const T& _y, const T& _z, const T& _w)
{
    array[0] = _x;
    array[1] = _y;
    array[2] = _z;
    array[3] = _w;
}

template <size_t M, typename T>
vector<M, T>::vector(const T* values)
{
    memcpy(array, values, M * sizeof(T));
}

#ifdef __OSG_MATH
template <size_t M, typename T>
template <typename OSGVEC3>
vector<M, T>::vector(const OSGVEC3& from,
                     typename std::enable_if<M == 3, OSGVEC3>::type*)
{
    array[0] = from.x();
    array[1] = from.y();
    array[2] = from.z();
}
#endif

#ifndef SWIG
template <size_t M, typename T>
template <typename TT>
// initializes the first M-1 values from vector_, the last from last_
vector<M, T>::vector(const vector<M - 1, TT>& vector_, T last_,
                     typename std::enable_if<M == 4, TT>::type*)
{
    array[0] = vector_.array[0];
    array[1] = vector_.array[1];
    array[2] = vector_.array[2];
    array[3] = last_;
}
#endif

// to-homogenous-coordinates ctor
template <size_t M, typename T>
template <typename TT>
vector<M, T>::vector(const vector<3, TT>& source_,
                     typename std::enable_if<M == 4, TT>::type*)
{
    std::copy(source_.begin(), source_.end(), begin());
    at(M - 1) = static_cast<T>(1.0);
}

// from-homogenous-coordinates ctor
template <size_t M, typename T>
template <typename TT>
vector<M, T>::vector(const vector<4, TT>& source_,
                     typename std::enable_if<M == 3, TT>::type*)
{
    const T w_reci = static_cast<T>(1.0) / source_(M);
    iterator it = begin(), it_end = end();
    for (size_t index = 0; it != it_end; ++it, ++index)
        *it = source_(index) * w_reci;
}

template <size_t M, typename T>
template <typename U>
vector<M, T>::vector(const vector<M, U>& source_)
{
    (*this) = source_;
}

namespace
{
template <size_t M, typename T>
vector<M, T> _createVector(const T x, const T y, const T z,
                           typename std::enable_if<M == 3>::type* = nullptr)
{
    return vector<M, T>(x, y, z);
}

template <size_t M, typename T>
vector<M, T> _createVector(const T x, const T y, const T z,
                           typename std::enable_if<M == 4>::type* = nullptr)
{
    return vector<M, T>(x, y, z, 1);
}
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::zero()
{
    return _createVector<M, T>(0, 0, 0);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::forward()
{
    return _createVector<M, T>(0, 0, -1);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::backward()
{
    return _createVector<M, T>(0, 0, 1);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::up()
{
    return _createVector<M, T>(0, 1, 0);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::down()
{
    return _createVector<M, T>(0, -1, 0);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::left()
{
    return _createVector<M, T>(-1, 0, 0);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::right()
{
    return _createVector<M, T>(1, 0, 0);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::unitX()
{
    return _createVector<M, T>(1, 0, 0);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::unitY()
{
    return _createVector<M, T>(0, 1, 0);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::unitZ()
{
    return _createVector<M, T>(0, 0, 1);
}

template <size_t M, typename T>
void vector<M, T>::set(T _a)
{
    for (iterator it = begin(), it_end = end(); it != it_end; ++it)
        *it = _a;
}

template <size_t M, typename T>
template <size_t N>
void vector<M, T>::set(const vector<N, T>& v)
{
    size_t minimum = M;
    if (N < M)
        minimum = N;
    memcpy(array, v.array, sizeof(T) * minimum);
}

template <size_t M, typename T>
void vector<M, T>::set(T _x, T _y)
{
    array[0] = _x;
    array[1] = _y;
}

template <size_t M, typename T>
void vector<M, T>::set(T _x, T _y, T _z)
{
    array[0] = _x;
    array[1] = _y;
    array[2] = _z;
}

template <size_t M, typename T>
void vector<M, T>::set(T _x, T _y, T _z, T _w)
{
    array[0] = _x;
    array[1] = _y;
    array[2] = _z;
    array[3] = _w;
}

template <size_t M, typename T>
inline T& vector<M, T>::operator()(size_t index)
{
    return at(index);
}

template <size_t M, typename T>
inline const T& vector<M, T>::operator()(size_t index) const
{
    return at(index);
}

template <size_t M, typename T>
inline T& vector<M, T>::at(size_t index)
{
    if (index >= M)
        throw std::runtime_error("at() - index out of bounds");
    return array[index];
}

template <size_t M, typename T>
inline const T& vector<M, T>::at(size_t index) const
{
    if (index >= M)
        throw std::runtime_error("at() - index out of bounds");
    return array[index];
}

template <size_t M, typename T>
vector<M, T>::operator T*()
{
    return array;
}

template <size_t M, typename T>
vector<M, T>::operator const T*() const
{
    return array;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator*(const vector<M, T>& other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) * other.at(index);
    return result;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator/(const vector<M, T>& other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) / other.at(index);
    return result;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator+(const vector<M, T>& other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) + other.at(index);
    return result;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator-(const vector<M, T>& other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) - other.at(index);
    return result;
}

template <size_t M, typename T>
void vector<M, T>::operator*=(const vector<M, T>& other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) *= other.at(index);
}

template <size_t M, typename T>
void vector<M, T>::operator/=(const vector<M, T>& other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) /= other.at(index);
}

template <size_t M, typename T>
void vector<M, T>::operator+=(const vector<M, T>& other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) += other.at(index);
}

template <size_t M, typename T>
void vector<M, T>::operator-=(const vector<M, T>& other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) -= other.at(index);
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator*(const T other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) * other;
    return result;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator/(const T other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) / other;
    return result;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator+(const T other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) + other;
    return result;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator-(const T other) const
{
    vector<M, T> result;
    for (size_t index = 0; index < M; ++index)
        result.at(index) = at(index) - other;
    return result;
}

template <size_t M, typename T>
void vector<M, T>::operator*=(const T other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) *= other;
}

template <size_t M, typename T>
void vector<M, T>::operator/=(const T other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) /= other;
}

template <size_t M, typename T>
void vector<M, T>::operator+=(const T other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) += other;
}

template <size_t M, typename T>
void vector<M, T>::operator-=(const T other)
{
    for (size_t index = 0; index < M; ++index)
        at(index) -= other;
}

template <size_t M, typename T>
vector<M, T> vector<M, T>::operator-() const
{
    vector<M, T> v(*this);
    return v.negate();
}

template <size_t M, typename T>
const vector<M, T>& vector<M, T>::negate()
{
    for (size_t index = 0; index < M; ++index)
        array[index] = -array[index];
    return *this;
}

template <size_t M, typename T>
inline T& vector<M, T>::x()
{
    return array[0];
}

template <size_t M, typename T>
inline T& vector<M, T>::y()
{
    return array[1];
}

template <size_t M, typename T>
inline T& vector<M, T>::z()
{
    return array[2];
}

template <size_t M, typename T>
inline T& vector<M, T>::w()
{
    return array[3];
}

template <size_t M, typename T>
inline const T& vector<M, T>::x() const
{
    return array[0];
}

template <size_t M, typename T>
inline const T& vector<M, T>::y() const
{
    return array[1];
}

template <size_t M, typename T>
inline const T& vector<M, T>::z() const
{
    return array[2];
}

template <size_t M, typename T>
inline const T& vector<M, T>::w() const
{
    return array[3];
}

template <size_t M, typename T>
inline T& vector<M, T>::r()
{
    return array[0];
}

template <size_t M, typename T>
inline T& vector<M, T>::g()
{
    return array[1];
}

template <size_t M, typename T>
inline T& vector<M, T>::b()
{
    return array[2];
}

template <size_t M, typename T>
inline T& vector<M, T>::a()
{
    return array[3];
}

template <size_t M, typename T>
inline const T& vector<M, T>::r() const
{
    return array[0];
}

template <size_t M, typename T>
inline const T& vector<M, T>::g() const
{
    return array[1];
}

template <size_t M, typename T>
inline const T& vector<M, T>::b() const
{
    return array[2];
}

template <size_t M, typename T>
inline const T& vector<M, T>::a() const
{
    return array[3];
}

template <size_t M, typename T>
template <typename TT>
vector<M, T>& vector<M, T>::cross(const vector<M, TT>& rhs,
                                  typename std::enable_if<M == 3, TT>::type*)
{
    const T x_ = y() * rhs.z() - z() * rhs.y();
    const T y_ = z() * rhs.x() - x() * rhs.z();
    const T z_ = x() * rhs.y() - y() * rhs.x();
    x() = x_;
    y() = y_;
    z() = z_;
    return *this;
}

template <size_t M, typename T>
inline T vector<M, T>::dot(const vector<M, T>& other) const
{
    T tmp = 0.0;
    for (size_t index = 0; index < M; ++index)
        tmp += at(index) * other.at(index);

    return tmp;
}

template <size_t M, typename T>
inline T vector<M, T>::normalize()
{
    const T len = length();
    if (len <= std::numeric_limits<T>::epsilon())
        return 0;

    const T tmp = 1.0 / len;
    (*this) *= tmp;
    return len;
}

template <size_t M, typename T>
inline T vector<M, T>::length() const
{
    return std::sqrt(squared_length());
}

template <size_t M, typename T>
inline T vector<M, T>::squared_length() const
{
    T _squared_length = 0.0;
    for (const_iterator it = begin(), it_end = end(); it != it_end; ++it)
        _squared_length += (*it) * (*it);

    return _squared_length;
}

template <size_t M, typename T>
inline T vector<M, T>::distance(const vector<M, T>& other) const
{
    return std::sqrt(squared_distance(other));
}

template <size_t M, typename T>
inline T vector<M, T>::squared_distance(const vector<M, T>& other) const
{
    vector<M, T> tmp(*this);
    tmp -= other;
    return tmp.squared_length();
}

template <size_t M, typename T>
inline T vector<M, T>::product() const
{
    T result = at(0);
    for (size_t i = 1; i < M; ++i)
        result *= at(i);
    return result;
}

template <size_t M, typename T>
template <typename TT>
vector<3, T>& vector<M, T>::rotate(const T theta, vector<M, TT> axis,
                                   typename std::enable_if<M == 3, TT>::type*)
{
    const T costheta = std::cos(theta);
    const T sintheta = std::sin(theta);

    axis.normalize();
    return *this = vector<3, T>(
               (costheta + (1 - costheta) * axis.x() * axis.x()) * x() +
                   ((1 - costheta) * axis.x() * axis.y() -
                    axis.z() * sintheta) *
                       y() +
                   ((1 - costheta) * axis.x() * axis.z() +
                    axis.y() * sintheta) *
                       z(),

               ((1 - costheta) * axis.x() * axis.y() + axis.z() * sintheta) *
                       x() +
                   (costheta + (1 - costheta) * axis.y() * axis.y()) * y() +
                   ((1 - costheta) * axis.y() * axis.z() -
                    axis.x() * sintheta) *
                       z(),

               ((1 - costheta) * axis.x() * axis.z() - axis.y() * sintheta) *
                       x() +
                   ((1 - costheta) * axis.y() * axis.z() +
                    axis.x() * sintheta) *
                       y() +
                   (costheta + (1 - costheta) * axis.z() * axis.z()) * z());
}

// sphere layout: center xyz, radius w
template <size_t M, typename T>
template <typename TT>
inline vector<3, T> vector<M, T>::project_point_onto_sphere(
    const vector<3, TT>& point,
    typename std::enable_if<M == 4, TT>::type*) const
{
    const vector<3, T>& center_ = get_sub_vector<3>();

    vector<3, T> projected_point(point);
    projected_point -= center_;
    projected_point.normalize();
    projected_point *= w();
    return center_ + projected_point;
}

// sphere layout: center xyz, radius w
template <size_t M, typename T>
template <typename TT>
inline T vector<M, T>::distance_to_sphere(
    const vector<3, TT>& point,
    typename std::enable_if<M == 4, TT>::type*) const
{
    const vector<3, T>& center_ = get_sub_vector<3>();
    return (point - center_).length() - w();
}

template <size_t M, typename T>
template <size_t N, size_t O>
vector<N, T> vector<M, T>::get_sub_vector(
    typename std::enable_if<M >= N + O>::type*) const
{
    return vector<N, T>(array + O);
}

template <size_t M, typename T>
template <size_t N, size_t O>
void vector<M, T>::set_sub_vector(const vector<N, T>& sub,
                                  typename std::enable_if<M >= N + O>::type*)
{
    ::memcpy(array + O, sub.array, N * sizeof(T));
}

// plane: normal xyz, distance w
template <size_t M, typename T>
template <typename TT>
inline T vector<M, T>::distance_to_plane(
    const vector<3, TT>& point,
    typename std::enable_if<M == 4, TT>::type*) const
{
    const vector<3, T>& normal = get_sub_vector<3>();
    return normal.dot(point) + w();
}

// plane: normal xyz, distance w
template <size_t M, typename T>
template <typename TT>
vector<3, T> vector<M, T>::project_point_onto_plane(
    const vector<3, TT>& point,
    typename std::enable_if<M == 4, TT>::type*) const
{
    const vector<3, T>& normal = get_sub_vector<3>();
    return point - (normal * distance_to_plane(point));
}

template <size_t M, typename T>
bool vector<M, T>::operator==(const vector<M, T>& other) const
{
    return memcmp(array, other.array, sizeof(array)) == 0;
}

template <size_t M, typename T>
bool vector<M, T>::operator!=(const vector<M, T>& other) const
{
    return !this->operator==(other);
}

template <size_t M, typename T>
bool vector<M, T>::equals(const vector<M, T>& other, T tolerance) const
{
    for (size_t index = 0; index < M; ++index)
        if (fabs(at(index) - other(index)) >= tolerance)
            return false;
    return true;
}

template <size_t M, typename T>
bool vector<M, T>::operator<(const vector<M, T>& other) const
{
    for (size_t index = 0; index < M; ++index)
    {
        if (at(index) < other.at(index))
            return true;
        if (other.at(index) < at(index))
            return false;
    }
    return false;
}

template <size_t M, typename T>
T vector<M, T>::operator=(T filler_value)
{
    for (size_t index = 0; index < M; ++index)
        at(index) = filler_value;
    return filler_value;
}

template <size_t M, typename T>
vector<M, T>& vector<M, T>::operator=(const vector<M, T>& other)
{
    if (this != &other)
        memcpy(array, other.array, M * sizeof(T));
    return *this;
}

// returns void to avoid 'silent' loss of precision when chaining
template <size_t M, typename T>
template <typename U>
void vector<M, T>::operator=(const vector<M, U>& source_)
{
    typedef typename vector<M, U>::const_iterator u_c_iter;
    u_c_iter it = source_.begin(), it_end = source_.end();
    for (iterator my_it = begin(); it != it_end; ++it, ++my_it)
        *my_it = static_cast<T>(*it);
}

template <size_t M, typename T>
template <typename input_iterator_t>
void vector<M, T>::iter_set(input_iterator_t begin_, input_iterator_t end_)
{
    input_iterator_t in_it = begin_;
    iterator it = begin(), it_end = end();
    for (; it != it_end && in_it != end_; ++it, ++in_it)
        (*it) = static_cast<T>(*in_it);
}

template <size_t M, typename T>
void vector<M, T>::clamp(const T& min, const T& max)
{
    for (size_t i = 0; i < M; ++i)
    {
        if (array[i] < min)
            array[i] = min;
        if (array[i] > max)
            array[i] = max;
    }
}

template <size_t M, typename T>
inline size_t vector<M, T>::size()
{
    return M;
}

template <size_t M, typename T>
size_t vector<M, T>::find_min_index() const
{
    return std::min_element(begin(), end()) - begin();
}

template <size_t M, typename T>
size_t vector<M, T>::find_max_index() const
{
    return std::max_element(begin(), end()) - begin();
}

template <size_t M, typename T>
T& vector<M, T>::find_min()
{
    return *std::min_element(begin(), end());
}

template <size_t M, typename T>
const T& vector<M, T>::find_min() const
{
    return *std::min_element(begin(), end());
}

template <size_t M, typename T>
T& vector<M, T>::find_max()
{
    return *std::max_element(begin(), end());
}

template <size_t M, typename T>
const T& vector<M, T>::find_max() const
{
    return *std::max_element(begin(), end());
}

template <size_t M, typename T>
inline typename vector<M, T>::iterator vector<M, T>::begin()
{
    return array;
}

template <size_t M, typename T>
inline typename vector<M, T>::iterator vector<M, T>::end()
{
    return array + M;
    ;
}

template <size_t M, typename T>
inline typename vector<M, T>::const_iterator vector<M, T>::begin() const
{
    return array;
}

template <size_t M, typename T>
inline typename vector<M, T>::const_iterator vector<M, T>::end() const
{
    return array + M;
    ;
}

template <size_t M, typename T>
inline typename vector<M, T>::reverse_iterator vector<M, T>::rbegin()
{
    return array + M - 1;
}

template <size_t M, typename T>
inline typename vector<M, T>::reverse_iterator vector<M, T>::rend()
{
    return array - 1;
}

template <size_t M, typename T>
inline typename vector<M, T>::const_reverse_iterator vector<M, T>::rbegin()
    const
{
    return array + M - 1;
}

template <size_t M, typename T>
inline typename vector<M, T>::const_reverse_iterator vector<M, T>::rend() const
{
    return array - 1;
}

template <size_t M, typename T>
bool vector<M, T>::is_unit_vector() const
{
    const_iterator it = begin(), it_end = end();
    bool one = false;
    for (; it != it_end; ++it)
    {
        if (*it == 1.0)
        {
            if (one)
                return false;
            one = true;
        }
        else if (*it != 0.0)
        {
            return false;
        }
    }
    return one;
}

template <size_t M, typename T>
void vector<M, T>::perturb(T perturbation)
{
    for (iterator it = begin(), it_end = end(); it != it_end; ++it)
    {
        (*it) += (rand() & 1u) ? perturbation : -perturbation;
    }
}

template <size_t M, typename T>
void vector<M, T>::sqrt_elementwise()
{
    for (iterator it = begin(), it_end = end(); it != it_end; ++it)
    {
        (*it) = std::sqrt(*it);
    }
}

template <size_t M, typename T>
void vector<M, T>::reciprocal()
{
    for (iterator it = begin(), it_end = end(); it != it_end; ++it)
        (*it) = static_cast<T>(1) / (*it);
}

template <size_t M, typename T>
void vector<M, T>::reciprocal_safe()
{
    for (iterator it = begin(), it_end = end(); it != it_end; ++it)
    {
        T& v = *it;

        if (v == 0)
            v = std::numeric_limits<T>::max();
        else
            v = static_cast<T>(1) / v;
    }
}

template <size_t M, typename T>
template <typename TT>
void vector<M, T>::cast_from(const vector<M, TT>& other)
{
    typedef vmml::vector<M, TT> vector_tt_type;
    typedef typename vector_tt_type::const_iterator tt_const_iterator;

    iterator it = begin(), it_end = end();
    tt_const_iterator other_it = other.begin();
    for (; it != it_end; ++it, ++other_it)
    {
        *it = static_cast<T>(*other_it);
    }
}

template <size_t M, typename T>
size_t vector<M, T>::nnz() const
{
    size_t counter = 0;

    const_iterator it = begin(), it_end = end();
    for (; it != it_end; ++it)
    {
        if (*it != 0)
        {
            ++counter;
        }
    }

    return counter;
}

template <size_t M, typename T>
double vector<M, T>::norm() const
{
    double norm_v = 0.0;

    const_iterator it = begin(), it_end = end();
    for (; it != it_end; ++it)
    {
        norm_v += *it * *it;
    }

    return std::sqrt(norm_v);
}

template <size_t M, typename T>
void vector<M, T>::set_random(int seed)
{
    if (seed >= 0)
        srand(seed);

    for (size_t i = 0; i < M; ++i)
    {
        const double fillValue = double(rand()) / double(RAND_MAX);
        at(i) = -1.0 + 2.0 * fillValue;
    }
}

} // namespace vmml

#endif
