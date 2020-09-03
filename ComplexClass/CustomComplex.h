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

#include <chrono>
#include <cmath>
#include <iostream>
using namespace std;
using namespace chrono;

template <class type>

class CustomComplex {
public:
  type re;
  type im;

  explicit CustomComplex() {
    re = 0.00;
    im = 0.00;
  }

  explicit CustomComplex(const double &a, const double &b) {
    re = a;
    im = b;
  }

  CustomComplex(const CustomComplex &src) {
    re = src.re;
    im = src.im;
  }

  CustomComplex &operator=(const CustomComplex &src) {
    re = src.re;
    im = src.im;

    return *this;
  }

  CustomComplex &operator+=(const CustomComplex &src) {
    re = src.re + this->re;
    im = src.im + this->im;

    return *this;
  }

  CustomComplex &operator*=(const CustomComplex &src) {
    re = src.re * this->re;
    im = src.im * this->im;

    return *this;
  }

  CustomComplex &operator*=(const double src) {
    re = src * this->re;
    im = src * this->im;

    return *this;
  }

  CustomComplex &operator-=(const CustomComplex &src) {
    re = src.re - this->re;
    im = src.im - this->im;

    return *this;
  }

  CustomComplex &operator-() {
    re = -this->re;
    im = -this->im;

    return *this;
  }

  CustomComplex &operator~() { return *this; }

  void print() const {
    printf("( %f, %f) ", this->re, this->im);
    printf("\n");
  }

  double real() const { return this->re; }

  double imag() const { return this->im; }

  void set_real(double val) { this->re = val; }

  void set_imag(double val) { this->im = val; }

  template <class T>
  friend inline CustomComplex<T> operator*(const CustomComplex<T> a,
                                           const CustomComplex<T> b) {
    return (
        CustomComplex<T>(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re));
  }

  template <class T>
  friend inline CustomComplex<T> operator*(const CustomComplex<T> &a,
                                           const double &b) {
    return (CustomComplex<T>(a.re * b, a.im * b));
  }

  template <class T>
  friend inline CustomComplex<T> operator*(const CustomComplex<T> &a,
                                           const int &b) {
    return (CustomComplex<T>(a.re * b, a.im * b));
  }

  template <class T>
  friend inline CustomComplex<T> operator-(CustomComplex<T> &a,
                                           CustomComplex<T> &b) {
    return (CustomComplex<T>(a.re - b.re, a.im - b.im));
  }

  template <class T>
  friend inline CustomComplex<T> operator-(T &a, CustomComplex<T> &src) {
    return (CustomComplex<T>(a - src.re, 0 - src.im));
  }

  template <class T>
  friend inline CustomComplex<T> operator+(const double &a,
                                           CustomComplex<T> &src) {
    return (CustomComplex<T>(a + src.re, src.im));
  }

  template <class T>
  friend inline CustomComplex<T> operator+(CustomComplex<T> a,
                                           CustomComplex<T> b) {
    return (CustomComplex<T>(a.re + b.re, a.im + b.im));
  }

  template <class T> inline CustomComplex<T> conj(const CustomComplex<T> *src) {
    return (CustomComplex<T>(src->re, -src->im));
  }

  template <class T> inline CustomComplex<T> conj(const CustomComplex<T> &src) {
    return (CustomComplex<T>(src.re, -src.im));
  }
};

/*Print the complex number*/
template <class T>
inline void CustomComplex_print(const CustomComplex<T> &src) {
  printf("( %f, %f) ", src.re, src.im);
  printf("\n");
}

template <class T>
inline CustomComplex<T> CustomComplex_conj(const CustomComplex<T> &src) {
  T re_this = src.re;
  T im_this = -1 * src.im;
  CustomComplex<T> result(re_this, im_this);
  return result;
}

#endif
