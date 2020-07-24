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

   Class library implementing simple MD-array structures similar to Fortran,
   except since this is c++ the fast running index are the rightmost
   ( For this "CPU" version the array elements are accessed in row-major format.
   )

   These are the necessary set used for BerkeleyGW mini-apps, and may not be
   completely defined for all combinations of sub-dimensional references
   and copy constructors.

   2D through 4D arrays are supported,
       ArrayxD<type> newarray;

   The '()' operator is overlaaded here to access array elements.

   Copy constructors are provided for openmp private/firstprivate.
   In this case, the data pointer must be allocated with resize or
   assignment as above.

   When using the GPU, the data pointers must be deep-copied to be made available on the device.
*/
#ifndef _ARRAYMDCPU_H
#define _ARRAYMDCPU_H

/* ----------------------------------------------------------------
   Class library implementing simple array structures similar to Fortran,
   except since this is c++ the fast running index are the rightmost
   ( For this "CPU" version. )

   Currently each of the 2D through 4D arrays are implemented in individual classes.
   Using variadic templates to combine all this in a single class has poor performance inside OpenMPTarget region

   The '()' operator is overlaaded here to access array elements.

   Copy constructors are provided for openmp private/firstprivate.
   In this case, the data pointer must be copied specially through deep copy for OpenMP4.5 and OpenACC data mapping.

   Rahulkumar Gayatri, NERSC, LBNL.
   ----------------------------------------------------------------*/

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

template<typename T>
struct Array1D
{
  unsigned n1;
  unsigned size;
  T* dptr;

  inline T& operator()(unsigned i1) { return dptr[i1]; }

  Array1D() = default;

  Array1D(const Array1D& p)
  {
    n1 = p.n1;
    size = 0;
    dptr = p.dptr;
  }

  Array1D(int in1)
  {
    n1 = in1;
    size = n1;
    dptr = new T[size];
  }

  ~Array1D()
  {
    if (size && dptr)
      delete[] dptr;
  }

  unsigned getSize() { return size * sizeof(T); }  // NB: in bytes
};

template<typename T>
struct Array2D
{
  unsigned n1, n2, b1;
  unsigned size;
  T* dptr;

  inline T& operator()(unsigned i1, unsigned i2)
  {
    return dptr[i2 + (n2 * i1)];
  }

  Array2D()
  {
    n1 = n2 = 0, b1 = 0;
    size = 0;
    dptr = NULL;
  }

  Array2D(const Array2D& p)
  {
    n1 = p.n1;
    n2 = p.n2;
    size = 0;
    dptr = p.dptr;
  }

  Array2D(int in1, int in2)
  {
    n1 = in1;
    n2 = in2;
    size = n1 * n2;
    dptr = new T[size];
  }

  ~Array2D()
  {
    if (size && dptr)
      delete[] dptr;
  }

  unsigned getSize() { return size * sizeof(T); }  // NB: in bytes
};

template<typename T>
struct Array3D
{
  unsigned n1, n2, n3;
  unsigned size;
  T* dptr;

  inline T& operator()(unsigned i1, unsigned i2, unsigned i3)
  {
    return dptr[i3 + i2 * f2 + i1 * f1];
  }

  Array3D() = default;

  Array3D(const Array3D& p)
  {
    n1 = p.n1;
    n2 = p.n2;
    n3 = p.n3;
    size = 0;
    dptr = p.dptr;
    f2 = n3;
    f1 = f2 * n2;
  }

  Array3D(unsigned in1, unsigned in2, unsigned in3)
  {
    n1 = in1;
    n2 = in2;
    n3 = in3;
    size = n1 * n2 * n3;
    f2 = n3;
    f1 = f2 * n2;
    dptr = new T[size];
  }

  ~Array3D()
  {
    if (size && dptr)
      delete[] dptr;
  }

  unsigned getSize() { return size * sizeof(T); }  // NB: in bytes

private:
  unsigned f2, f1, b1, b2;
};

template<typename T>
struct Array4D
{
  unsigned n1, n2, n3, n4;
  unsigned size;
  T* dptr;

  inline T& operator()(unsigned i1, unsigned i2, unsigned i3, unsigned i4)
  {
    return dptr[i4 + i3 * f3 + i2 * f2 + i1 * f1];
  }

  Array4D() = default;

  Array4D(const Array4D& p)
  {
    n1 = p.n1;
    n2 = p.n2;
    n3 = p.n3;
    n4 = p.n4;
    size = 0;
    dptr = p.dptr;
    f3 = n4;
    f2 = f3 * n3;
    f1 = f2 * n2;
  }

  Array4D(unsigned in1, unsigned in2, unsigned in3, unsigned in4)
  {
    n1 = in1;
    n2 = in2;
    n3 = in3;
    n4 = in4;
    size = n1 * n2 * n3 * n4;
    f3 = n4;
    f2 = f3 * n3;
    f1 = f2 * n2;
    dptr = new T[size];
  }

  ~Array4D()
  {
    if (size && dptr)
      delete[] dptr;
  }

private:
  unsigned f3, f2, f1, b1, b2, b3;
};

#endif
