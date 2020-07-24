/*!=============================================================================================
Unless otherwise stated, all files distributed in this package
are licensed under the following terms.
BerkeleyGW, Copyright (c) 2011, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy).
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
(1) Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence
Berkeley National Laboratory, U.S. Dept. of Energy nor the names of
its contributors may be used to endorse or promote products derived
from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
You are under no obligation whatsoever to provide any bug fixes,
patches, or upgrades to the features, functionality or performance of
the source code ("Enhancements") to anyone; however, if you choose to
make your Enhancements available either publicly, or directly to
Lawrence Berkeley National Laboratory, without imposing a separate
written license agreement for such Enhancements, then you hereby grant
the following license: a  non-exclusive, royalty-free perpetual
license to install, use, modify, prepare derivative works, incorporate
into other computer software, distribute, and sublicense such
enhancements or derivative works thereof, in binary and source code
form.

C++ class for complex number arithemetic.
*/
#ifndef __CustomComplex
#define __CustomComplex

#include <iostream>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <sys/time.h>
using namespace std;

template<class type>

class CustomComplex {
    public:
    type x;
    type y;


    explicit CustomComplex () {
        x = 0.00;
        y = 0.00;
    }


     explicit CustomComplex(const double& a, const double& b) {
        x = a;
        y = b;
    }

     CustomComplex(const CustomComplex& src) {
        x = src.x;
        y = src.y;
    }

     CustomComplex& operator =(const CustomComplex& src) {
        x = src.x;
        y = src.y;

        return *this;
    }

     CustomComplex& operator +=(const CustomComplex& src) {
        x = src.x + this->x;
        y = src.y + this->y;

        return *this;
    }

     CustomComplex& operator *=(const CustomComplex& src) {
        x = src.x * this->x;
        y = src.y * this->y;

        return *this;
    }

     CustomComplex& operator *=(const double src) {
        x = src * this->x;
        y = src * this->y;

        return *this;
    }

     CustomComplex& operator -=(const CustomComplex& src) {
        x = src.x - this->x;
        y = src.y - this->y;

        return *this;
    }

     CustomComplex& operator -() {
        x = -this->x;
        y = -this->y;

        return *this;
    }

    CustomComplex& operator ~() {
        return *this;
    }

    void print() const {
        printf("( %f, %f) ", this->x, this->y);
        printf("\n");
    }

    double get_real() const
    {
        return this->x;
    }

    double get_imag() const
    {
        return this->y;
    }

    void set_real(double val)
    {
        this->x = val;
    }

    void set_imag(double val)
    {
        this->y = val;
    }

    template<class T>
     friend inline CustomComplex<T> operator *(const CustomComplex<T> a, const CustomComplex<T> b) {
        T x_this = a.x * b.x - a.y*b.y ;
        T y_this = a.x * b.y + a.y*b.x ;
        CustomComplex<T> result(x_this, y_this);
        return (result);
    }

    template<class T>
     friend inline CustomComplex<T> operator *(const CustomComplex<T>& a, const double &b) {
       CustomComplex<T> result(a.x*b, a.y*b);
       return result;
    }

    template<class T>
     friend inline CustomComplex<T> operator *(const CustomComplex<T> &a, const int &b) {
       CustomComplex<T> result(a.x*b, a.y*b);
       return result;
    }

    template<class T>
     friend inline CustomComplex<T> operator -(CustomComplex<T>& a, CustomComplex<T>& b) {
        CustomComplex<T> result(a.x - b.x, a.y - b.y);
        return result;
    }

    template<class T>
     friend inline CustomComplex<T> operator -(T &a, CustomComplex<T>& src) {
        CustomComplex<T> result(a - src.x, 0 - src.y);
        return result;
    }

    template<class T>
     friend inline CustomComplex<T> operator +(const double &a, CustomComplex<T>& src) {
        CustomComplex<T> result(a + src.x, src.y);
        return result;
    }

    template<class T>
     friend inline CustomComplex<T> operator +(CustomComplex<T> a, CustomComplex<T> b) {
        CustomComplex<T> result(a.x + b.x, a.y+b.y);
        return result;
    }

    template<class T>
     friend inline CustomComplex<T> operator /(CustomComplex<T> a, CustomComplex<T> b) {

        CustomComplex<T> b_conj = CustomComplex_conj(&b);
        CustomComplex<T> numerator = a * b_conj;
        CustomComplex<T> denominator = b * b_conj;

        double re_this = numerator.x / denominator.x;
        double im_this = numerator.y / denominator.x;

        CustomComplex<T> result(re_this, im_this);
        return result;
    }

    template<class T>
     friend inline CustomComplex<T> operator /(CustomComplex<T> a, T b) {
       CustomComplex<T> result(a.x/b, a.y/b);
       return result;
    }

    template<class T>
     friend inline void CustomComplex_equals(const CustomComplex<T>* src, CustomComplex<T>* dest) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_conj(const CustomComplex<T>* src) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_conj(const CustomComplex<T>& src) ;

    template<class T>
     friend inline double CustomComplex_abs(const CustomComplex<T>& src) ;

    template<class T>
     friend inline double CustomComplex_real( const CustomComplex<T>* src) ;

    template<class T>
     friend inline double CustomComplex_imag( const CustomComplex<T>* src) ;

    template<class T>
     friend inline double CustomComplex_real( const CustomComplex<T>& src) ;

    template<class T>
     friend inline double CustomComplex_imag( const CustomComplex<T>& src) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_product(const CustomComplex<T>* src, T* b) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_product(const CustomComplex<T>* src, T b) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_product(const CustomComplex<T>* a, const CustomComplex<T>* b) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_minus(const T* a, const CustomComplex<T>* src) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_minus(const CustomComplex<T>* a, const CustomComplex<T>* b) ;

    template<class T>
     friend inline CustomComplex<T> CustomComplex_plus(const CustomComplex<T>* a, const CustomComplex<T>* b) ;

    template<class T>
     friend inline void CustomComplex_plusEquals(CustomComplex<T>* a, const CustomComplex<T>* b) ;

    template<class T>
     friend inline void CustomComplex_minusEquals(CustomComplex<T>* a, const CustomComplex<T>* b) ;

    template<class T>
    friend inline void CustomComplex_print(const CustomComplex<T>& src);
};

/*Print the complex number*/
template<class T>
inline void CustomComplex_print(const CustomComplex<T>& src)
{
    printf("( %f, %f) ", src.x, src.y);
    printf("\n");
}

/* Return the conjugate of a complex number
flop
*/
template<class T>
inline CustomComplex<T> CustomComplex_conj(const CustomComplex<T>* src) {
    T re_this = src->x;
    T im_this = -1 * src->y;
    CustomComplex<T> result(re_this, im_this);
    return result;
}

template<class T>
inline CustomComplex<T> CustomComplex_conj(const CustomComplex<T>& src) {
    T re_this = src.x;
    T im_this = -1 * src.y;
    CustomComplex<T> result(re_this, im_this);
    return result;
}

/*
 * Return the absolute of a complex number
 */
template<class T>
inline double CustomComplex_abs(const CustomComplex<T>& src) {
    T re_this = src.x * src.x;
    T im_this = src.y * src.y;

    T result = sqrt(re_this+im_this);
    return result;
}

/*
 * Return the real part of a complex number
 */
template<class T>
inline double CustomComplex_real( const CustomComplex<T>* src) {
    return src->x;
}

template<class T>
inline double CustomComplex_real( const CustomComplex<T>& src) {
    return src.x;
}

/*
 * Return the imaginary part of a complex number
 */
template<class T>
inline double CustomComplex_imag( const CustomComplex<T>* src) {
    return src->y;
}

template<class T>
inline double CustomComplex_imag( const CustomComplex<T>& src) {
    return src.y;
}

template<class T>
inline CustomComplex<T> CustomComplex_product(const CustomComplex<T>* src, T* b) {
   T re_this = src->x * (*b);
   T im_this = src->y * (*b);
   return (CustomComplex<T>(re_this, im_this));
}

template<class T>
inline CustomComplex<T> CustomComplex_product(const CustomComplex<T>* src, T b) {
   T re_this = src->x * b;
   T im_this = src->y * b;
   return (CustomComplex<T>(re_this, im_this));
}

template<class T>
inline CustomComplex<T> CustomComplex_product(const CustomComplex<T>* a, const CustomComplex<T>* b){
    T x_this = a->x * b->x - a->y*b->y ;
    T y_this = a->x * b->y + a->y*b->x ;
    CustomComplex<T> result(x_this, y_this);
    return (result);
}

template<class T>
inline CustomComplex<T> CustomComplex_minus(const CustomComplex<T>* a, const CustomComplex<T>* b){
        CustomComplex<T> result(a->x - b->x, a->y - b->y);
        return result;
}

template<class T>
inline CustomComplex<T> CustomComplex_minus(const T* a, const CustomComplex<T>* src) {
        CustomComplex<T> result(*a - src->x, 0 - src->y);
        return result;
}

template<class T>
inline CustomComplex<T> CustomComplex_plus(const CustomComplex<T>* a, const CustomComplex<T>* b){
        CustomComplex<T> result(a->x + b->x, a->y + b->y);
        return result;
}

template<class T>
inline void CustomComplex_equals(const CustomComplex<T>* src, CustomComplex<T>* dest) {
    *dest = CustomComplex<T>(src->x, src->y);
}

template<class T>
inline void CustomComplex_plusEquals(CustomComplex<T>* a, const CustomComplex<T>* b){
        a->x += b->x ;
        a->y += b->y ;
}

template<class T>
inline void CustomComplex_minusEquals(CustomComplex<T>* a, const CustomComplex<T>* b){
        a->x -= b->x ;
        a->y -= b->y ;
}

#endif

