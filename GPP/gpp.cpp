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
!
! This Kernel Represents the GPP Self-Energy Summations in BerkeleyGW
!
! Run like:
!
!  gppkernel.x <number_bands> <number_valence_bands> <number_plane_waves>
<nodes_per_mpi_group> <gppsum>
!
! Example run:
!
!  gppkernel.x 8 2 10000 20
!
!==============================================================================================
*/

#include "../common/Defines.h"

#if defined(_OPENACC)
inline void init_structs(size_t number_bands, size_t ngpown, size_t ncouls,
                         ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                         ARRAY2D &I_eps_array, ARRAY2D &wtilde_array,
                         ARRAY1D_DataType &vcoul, ARRAY1D_int &inv_igp_index,
                         ARRAY1D_int &indinv, ARRAY1D_DataType &wx_array) {
  const DataType dw = 1;
  const DataType e_lk = 10;
  const DataType to1 = 1e-6;
  const DataType limittwo = pow(0.5, 2);
  const DataType e_n1kq = 6.0;
  ComplexType expr(0.5, 0.5);

#pragma acc enter data copyin(aqsmtemp, aqsntemp, vcoul, inv_igp_index,        \
                              indinv, I_eps_array, wtilde_array, wx_array)

#pragma acc enter data create(                                                 \
    aqsmtemp.dptr [0:aqsmtemp.size], vcoul.dptr [0:vcoul.size],                \
    inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],    \
    aqsntemp.dptr [0:aqsntemp.size], I_eps_array.dptr [0:I_eps_array.size],    \
    wx_array.dptr [nstart:nend], wtilde_array.dptr [0:wtilde_array.size])

#pragma acc parallel loop present(aqsmtemp, aqsntemp)
  for (int i = 0; i < number_bands; i++)
    for (int j = 0; j < ncouls; j++) {
      aqsmtemp(i, j) = ComplexType(0.5, 0.5);
      aqsntemp(i, j) = ComplexType(0.5, 0.5);
    }

#pragma acc parallel loop copyin(expr) present(I_eps_array, wtilde_array)
  for (int i = 0; i < ngpown; i++)
    for (int j = 0; j < ncouls; j++) {
      I_eps_array(i, j) = expr;
      wtilde_array(i, j) = expr;
    }

#pragma acc parallel loop present(vcoul)
  for (int i = 0; i < ncouls; i++)
    vcoul(i) = 1.0;

#pragma acc parallel loop present(inv_igp_index)
  for (int ig = 0; ig < ngpown; ++ig)
    inv_igp_index(ig) = (ig + 1) * ncouls / ngpown;

#pragma acc parallel loop present(indinv)
  for (int ig = 0; ig < ncouls; ++ig)
    indinv(ig) = ig;
  indinv(ncouls) = ncouls - 1;

#pragma acc parallel loop present(wx_array)
  for (int iw = nstart; iw < nend; ++iw) {
    wx_array(iw) = e_lk - e_n1kq + dw * ((iw + 1) - 2);
    if (wx_array(iw) < to1)
      wx_array(iw) = to1;
  }
}
#endif

inline void correntess(ComplexType result) {
  double re_diff, im_diff;
#if defined(_OPENMP) || defined(_OPENACC)
  re_diff = result.real() - -264241149.849658;
  im_diff = result.imag() - 1321205773.349384;
#else
  re_diff = result.real() - -264241220.914570;
  im_diff = result.imag() - 1321205332.084101;
#endif

  if (re_diff < 0.001 && im_diff < 0.001)
    printf("\n!!!! SUCCESS - !!!! Correctness test passed :-D :-D\n\n");
  else
    printf("\n!!!! FAILURE - Correctness test failed :-( :-(  \n");
}

int main(int argc, char **argv) {

#if defined(OPENMP_TARGET)
  cout << "\n ************OpenMP 4.5**********\n" << endl;
#elif defined(_OPENMP)
  cout << "\n ************OpenMP 3.0**********\n" << endl;
#elif defined(_OPENACC)
  cout << "\n ************OpenACC  **********\n" << endl;
#endif

  int number_bands = 0, nvband = 0, ncouls = 0, nodes_per_group = 0;
  int npes = 1;
  if (argc == 1) {
    number_bands = 512;
    nvband = 2;
    ncouls = 32768;
    nodes_per_group = 20;
  } else if (argc == 5) {
    number_bands = atoi(argv[1]);
    nvband = atoi(argv[2]);
    ncouls = atoi(argv[3]);
    nodes_per_group = atoi(argv[4]);
  } else {
    std::cout << "The correct form of input is : " << endl;
    std::cout << " ./a.out <number_bands> <number_valence_bands> "
                 "<number_plane_waves> <nodes_per_mpi_group> "
              << endl;
    exit(0);
  }
  int ngpown = ncouls / (nodes_per_group * npes);

  // Constants that will be used later
  const DataType e_lk = 10;
  const DataType dw = 1;
  const DataType to1 = 1e-6;
  const DataType limittwo = pow(0.5, 2);
  const DataType e_n1kq = 6.0;

  // Using time point and system_clock
  time_point<system_clock> start, end, k_start, k_end;
  start = system_clock::now();
  double elapsedKernelTimer;

  // OpenMP Printing of threads on Host and Device
#if defined(_OPENMP)
  int numThreads = 1, numTeams = 2;
#if defined(OPENMP_TARGET)
#pragma omp target teams map(tofrom : numTeams, numThreads) \
  shared(numTeams)
  {
    int tid = omp_get_team_num();
    if (tid == 0) {
      numTeams = omp_get_num_teams();
#pragma omp parallel
      {
        int ttid = omp_get_thread_num();
        if (ttid == 0)
          numThreads = omp_get_num_threads();
      }
    }
  }
  std::cout << "Number of OpenMP Teams = " << numTeams << std::endl;
  std::cout << "Number of OpenMP DEVICE Threads = " << numThreads << std::endl;
#else
#pragma omp parallel shared(numThreads)
  {
    int tid = omp_get_thread_num();
    if (tid == 0)
      numThreads = omp_get_num_threads();
  }
  std::cout << "Number of OpenMP Threads = " << numThreads << endl;
#endif
#endif


  // Printing out the params passed.
  std::cout << "Sizeof(ComplexType = "
            << sizeof(ComplexType) << " bytes" << std::endl;
  std::cout << "number_bands = " << number_bands << "\t nvband = " << nvband
            << "\t ncouls = " << ncouls
            << "\t nodes_per_group  = " << nodes_per_group
            << "\t ngpown = " << ngpown << "\t nend = " << nend
            << "\t nstart = " << nstart << endl;

  size_t memFootPrint = 0.00;

  // ALLOCATE statements from fortran gppkernel.
  ARRAY1D achtemp(nend - nstart);
  memFootPrint += (nend - nstart) * sizeof(ComplexType);

  ARRAY2D aqsmtemp(number_bands, ncouls);
  ARRAY2D aqsntemp(number_bands, ncouls);
  memFootPrint += 2 * (number_bands * ncouls) * sizeof(ComplexType);

  ARRAY2D I_eps_array(ngpown, ncouls);
  ARRAY2D wtilde_array(ngpown, ncouls);
  memFootPrint += 2 * (ngpown * ncouls) * sizeof(ComplexType);

  ARRAY1D_DataType vcoul(ncouls);
  memFootPrint += ncouls * sizeof(DataType);

  ARRAY1D_int inv_igp_index(ngpown);
  ARRAY1D_int indinv(ncouls + 1);
  memFootPrint += ngpown * sizeof(int);
  memFootPrint += (ncouls + 1) * sizeof(int);

  ARRAY1D_DataType wx_array(nend - nstart);
  memFootPrint += 3 * (nend - nstart) * sizeof(DataType);

  // Print Memory Foot print
  cout << "Memory Foot Print = " << memFootPrint / pow(1024, 3) << " GBs"
       << endl;

#if !defined(_OPENACC)
  ComplexType expr(.5, .5);
  for (int i = 0; i < number_bands; i++)
    for (int j = 0; j < ncouls; j++) {
      aqsmtemp(i, j) = expr;
      aqsntemp(i, j) = expr;
    }

  for (int i = 0; i < ngpown; i++)
    for (int j = 0; j < ncouls; j++) {
      I_eps_array(i, j) = expr;
      wtilde_array(i, j) = expr;
    }

  for (int i = 0; i < ncouls; i++)
    vcoul(i) = 1.0;

  for (int ig = 0; ig < ngpown; ++ig)
    inv_igp_index(ig) = (ig + 1) * ncouls / ngpown;

  // Do not know yet what this array represents
  for (int ig = 0; ig < ncouls; ++ig)
    indinv(ig) = ig;
  indinv(ncouls) = ncouls - 1;

  for (int iw = nstart; iw < nend; ++iw) {
    wx_array(iw) = e_lk - e_n1kq + dw * ((iw + 1) - 2);
    if (wx_array(iw) < to1)
      wx_array(iw) = to1;
  }

#else
  init_structs(number_bands, ngpown, ncouls, aqsmtemp, aqsntemp, I_eps_array,
               wtilde_array, vcoul, inv_igp_index, indinv, wx_array);
#endif // Initailize OpenACC structures on the device directily

  k_start = system_clock::now();
  noflagOCC_solver(number_bands, ngpown, ncouls, inv_igp_index, indinv,
                   wx_array, wtilde_array, aqsmtemp, aqsntemp, I_eps_array,
                   vcoul, achtemp);

  k_end = system_clock::now();
  duration<double> elapsed = k_end - k_start;
  elapsedKernelTimer = elapsed.count();

  // Check for correctness
  correntess(achtemp(0));
  printf("\n Final achtemp\n");
  ComplexType_print(achtemp(0));

  end = system_clock::now();
  elapsed = end - start;

  cout << "********** Kernel Time Taken **********= " << elapsedKernelTimer
       << " secs" << endl;
  cout << "********** Total Time Taken **********= " << elapsed.count()
       << " secs" << endl;

  return 0;
}

void noflagOCC_solver(size_t number_bands, size_t ngpown, size_t ncouls,
                      ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                      ARRAY1D_DataType &wx_array, ARRAY2D &wtilde_array,
                      ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                      ARRAY2D &I_eps_array, ARRAY1D_DataType &vcoul,
                      ARRAY1D &achtemp) {
  time_point<system_clock> start, end;
  start = system_clock::now();
  // Vars to use for reduction
  DataType ach_re0 = 0.00, ach_re1 = 0.00, ach_re2 = 0.00, ach_im0 = 0.00,
           ach_im1 = 0.00, ach_im2 = 0.00;

#if defined(OPENMP_TARGET)
#pragma omp target enter data map(to                                           \
                                  : aqsmtemp, vcoul, inv_igp_index, indinv,    \
                                    aqsntemp, I_eps_array, wx_array,           \
                                    wtilde_array)

  start = system_clock::now();

#pragma omp target map(                                                        \
    to                                                                         \
    : aqsmtemp.dptr [0:aqsmtemp.size], vcoul.dptr [0:vcoul.size],              \
      inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],  \
      aqsntemp.dptr [0:aqsmtemp.size], I_eps_array.dptr [0:I_eps_array.size],  \
      wx_array.dptr [nstart:nend], wtilde_array.dptr [0:wtilde_array.size])    \
    map(tofrom                                                                 \
        : ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)

#pragma omp teams distribute parallel for collapse(2) \
    reduction(+:ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2) \
  num_teams(ngpown*ncouls) thread_limit(32)
  for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
    for (int ig = 0; ig < ncouls; ++ig) {
      for (int n1 = 0; n1 < number_bands; ++n1) {

#elif defined(_OPENMP)
#pragma omp parallel for collapse(2) \
    reduction(+:ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)
  for (int n1 = 0; n1 < number_bands; ++n1) {
    for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
      for (int ig = 0; ig < ncouls; ++ig) {
#elif defined(_OPENACC)
#pragma acc parallel loop gang vector collapse(3) \
    present(inv_igp_index, indinv, aqsmtemp, aqsntemp, wtilde_array, wx_array, I_eps_array, vcoul) \
    reduction(+:ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)\
    num_gangs(number_bands*ncouls)
  for (int n1 = 0; n1 < number_bands; ++n1) {
    for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
      for (int ig = 0; ig < ncouls; ++ig) {
#else
  for (int n1 = 0; n1 < number_bands; ++n1) {
    for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
      for (int ig = 0; ig < ncouls; ++ig) {
#endif
        int indigp = inv_igp_index(my_igp);
        int igp = indinv(indigp);
        DataType achtemp_re_loc[nend - nstart], achtemp_im_loc[nend - nstart];
        for (int iw = nstart; iw < nend; ++iw) {
          achtemp_re_loc[iw] = 0.00;
          achtemp_im_loc[iw] = 0.00;
        }
        ComplexType sch_store1 =
            ComplexType_conj(aqsmtemp(n1, igp)) * aqsntemp(n1, igp) * 0.5 *
            vcoul(igp) * wtilde_array(my_igp, igp);

        for (int iw = nstart; iw < nend; ++iw) {
          ComplexType wdiff =
              wx_array(iw) - wtilde_array(my_igp, ig);
          ComplexType delw =
              ComplexType_conj(wdiff) *
              (1 / (wdiff * ComplexType_conj(wdiff)).real());
          ComplexType sch_array =
              delw * I_eps_array(my_igp, ig) * sch_store1;

          achtemp_re_loc[iw] += (sch_array).real();
          achtemp_im_loc[iw] += (sch_array).imag();
        }
        ach_re0 += achtemp_re_loc[0];
        ach_re1 += achtemp_re_loc[1];
        ach_re2 += achtemp_re_loc[2];
        ach_im0 += achtemp_im_loc[0];
        ach_im1 += achtemp_im_loc[1];
        ach_im2 += achtemp_im_loc[2];
      }
    } // ngpown
  }   // number_bands

  end = system_clock::now();
  duration<double> elapsed = end - start;
  double elapsedKernelTimer = elapsed.count();

  achtemp(0) = ComplexType(ach_re0, ach_im0);
  achtemp(1) = ComplexType(ach_re1, ach_im1);
  achtemp(2) = ComplexType(ach_re2, ach_im2);
}
