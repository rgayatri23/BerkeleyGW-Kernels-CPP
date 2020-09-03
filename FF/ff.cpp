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
!=============================================================================================
!
! This Kernel Represents the Full-Frequency Self-Energy Summations in BerkeleyGW
!
! Run like:
!
! ffkernel.x <number_bands> <number_valence_bands> <number_plane_waves>
<number_plane_waves_per_mpi_task> <number_frequencies>
<number_evalutation_energies>
!
! Example of bandwidth bound run (large system fewer frequencies):
!
!  ffkernel.x 8 2 10000 600 400 200
!
! Example of compute bound run (small system asking for many frequencies):
!  ffkernel.x 8 2 1000 60 20000 30000
!
!==============================================================================================
*/
#include "../common/Defines.h"

inline double elapsedTime(timeval start_time, timeval end_time) {
  return ((end_time.tv_sec - start_time.tv_sec) +
          1e-6 * (end_time.tv_usec - start_time.tv_usec));
}

inline void correctness(ComplexType &achsDtemp, ComplexType &asxDtemp,
                        ComplexType &achDtemp_cor) {
  double re_diff = 0.0, im_diff = 0.0;
  // achsDtemp correctness
  re_diff = achsDtemp.real() - -2776374339.750000;
  im_diff = achsDtemp.imag() - 2776374339.750000;

  if (re_diff < 0.00001 && im_diff < 0.00001)
    printf("\n!!!! achsDtemp : SUCCESS - !!!! Correctness test passed :-D "
           ":-D\n\n");
  else
    printf("\n!!!! achsDtemp : FAILURE - Correctness test failed :-( :-(  \n");

  // asxDtemp correctness
  re_diff = asxDtemp.real() - 738493767.000000;
  im_diff = asxDtemp.imag() - -730362204.000000;

  if (re_diff < 0.00001 && im_diff < 0.00001)
    printf(
        "\n!!!! axsDtemp : SUCCESS - !!!! Correctness test passed :-D :-D\n\n");
  else
    printf("\n!!!! axsDtemp : FAILURE - Correctness test failed :-( :-(  \n");

  // achDtemp_cor correctness
  re_diff = achDtemp_cor.real() - 0.00;
  im_diff = achDtemp_cor.imag() - -734427985.500000;

  if (re_diff < 0.00001 && im_diff < 0.00001)
    printf("\n!!!! achDtemp_cor : SUCCESS - !!!! Correctness test passed :-D "
           ":-D\n\n");
  else
    printf(
        "\n!!!! achDtemp_cor : FAILURE - Correctness test failed :-( :-(  \n");
}

inline void compute_fact(double wx, int nFreq, ARRAY1D_DataType &dFreqGrid,
                         double &fact1, double &fact2, int &ifreq, int loop,
                         bool flag_occ) {
  if (loop == 1 && wx > 0.00) {
    for (int ijk = 0; ijk < nFreq - 1; ++ijk) {
      if (wx > dFreqGrid(ijk) && wx < dFreqGrid(ijk + 1))
        ifreq = ijk;
    }
    if (ifreq == 0)
      ifreq = nFreq - 2;
    fact1 =
        (dFreqGrid(ifreq + 1) - wx) / (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
    fact2 = (wx - dFreqGrid(ifreq)) / (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
  } else if (loop == 1) {
    for (int ijk = 0; ijk < nFreq - 1; ++ijk) {
      if (-wx > dFreqGrid(ijk) && -wx < dFreqGrid(ijk + 1))
        ifreq = ijk;
    }
    if (ifreq == 0)
      ifreq = nFreq - 2;
    fact1 =
        (dFreqGrid(ifreq + 1) + wx) / (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
    fact2 =
        (-dFreqGrid(ifreq) - wx) / (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
  }
  if (loop == 2 && wx > 0.00) {
    for (int ijk = 0; ijk < nFreq - 1; ++ijk) {
      if (wx > dFreqGrid(ijk) && wx < dFreqGrid(ijk + 1))
        ifreq = ijk;
    }
    if (ifreq == -1)
      ifreq = nFreq - 2;
    fact1 = -0.5 * (dFreqGrid(ifreq + 1) - wx) /
            (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
    fact2 = -0.5 * (wx - dFreqGrid(ifreq)) /
            (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
  } else if (loop == 2 && flag_occ) {
    wx = -wx;
    ifreq = 0;
    for (int ijk = 0; ijk < nFreq - 1; ++ijk) {
      if (wx > dFreqGrid(ijk) && wx < dFreqGrid(ijk + 1))
        ifreq = ijk;
    }
    if (ifreq == 0)
      ifreq = nFreq - 2;
    fact1 =
        (dFreqGrid(ifreq + 1) - wx) / (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
    fact2 = (wx - dFreqGrid(ifreq)) / (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq));
  }
}

inline void ssxDittt_kernel(ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                            ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                            ARRAY1D_DataType &vcoul, ARRAY3D &I_eps_array,
                            ComplexType &ssxDittt, int ngpown, int ncouls,
                            int n1, int ifreq, double fact1, double fact2,
                            int igp, int my_igp) {
  ComplexType ssxDitt(0.00, 0.00);
  for (int ig = 0; ig < ncouls; ++ig) {
    ComplexType ssxDit = I_eps_array(ifreq, my_igp, ig) * fact1 +
                         I_eps_array((ifreq + 1), my_igp, ig) * fact2;

    ssxDitt += aqsntemp(n1, ig) * ComplexType_conj(aqsmtemp(n1, igp)) * ssxDit *
               vcoul(igp);
  }
  ssxDittt = ssxDitt;
}

void achsDtemp_Kernel(int number_bands, int ngpown, int ncouls,
                      ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                      ARRAY2D &aqsntemp, ARRAY2D &aqsmtemp,
                      ARRAY3D &I_epsR_array, ARRAY1D_DataType &vcoul,
                      ComplexType &achsDtemp) {
  double achsDtemp_re = 0.00, achsDtemp_im = 0.00;
#if defined(OPENMP_TARGET)
#pragma omp target teams distribute parallel for collapse(2) \
    reduction(+:achsDtemp_re, achsDtemp_im)

#elif defined(_OPENMP)
#pragma omp parallel for collapse(2) default(shared) \
    reduction(+:achsDtemp_re, achsDtemp_im)

#elif defined(_OPENACC)
#pragma acc parallel loop gang collapse(2) \
    present(inv_igp_index, indinv, aqsntemp, aqsmtemp, \
            I_epsR_array, vcoul) \
    reduction(+:achsDtemp_re, achsDtemp_im)
#endif
  for (int n1 = 0; n1 < number_bands; ++n1) {
    for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
      int indigp = inv_igp_index(my_igp);
      int igp = indinv(indigp);
      ComplexType schsDtemp(0.00, 0.00);

      for (int ig = 0; ig < ncouls; ++ig) {
        ComplexType almost_schsDtemp = aqsntemp(n1, ig) *
                                       ComplexType_conj(aqsmtemp(n1, igp)) *
                                       I_epsR_array(1, my_igp, ig);
        schsDtemp = schsDtemp - almost_schsDtemp;
      }

      achsDtemp_re += (schsDtemp).real() * vcoul(igp) * 0.5;
      achsDtemp_im += (schsDtemp).imag() * vcoul(igp) * 0.5;
    } // my_igp
  }   // n1
  achsDtemp = ComplexType(achsDtemp_re, achsDtemp_im);
}

void asxDtemp_Kernel(int nvband, int nfreqeval, int ncouls, int ngpown,
                     int nFreq, double freqevalmin, double freqevalstep,
                     double occ, ARRAY1D_DataType &ekq,
                     ARRAY1D_DataType &dFreqGrid, ARRAY1D_int &inv_igp_index,
                     ARRAY1D_int &indinv, ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                     ARRAY1D_DataType &vcoul, ARRAY3D &I_epsR_array,
                     ARRAY3D &I_epsA_array, ARRAY1D &asxDtemp) {
  double *asxDtemp_re = new double[nfreqeval];
  double *asxDtemp_im = new double[nfreqeval];
  for (int iw = 0; iw < nfreqeval; ++iw) {
    asxDtemp_re[iw] = 0.00;
    asxDtemp_im[iw] = 0.00;
  }

#if defined(OPENMP_TARGET)

#pragma omp target enter data map(to                                           \
                                  : asxDtemp_re [0:nfreqeval],                 \
                                    asxDtemp_im [0:nfreqeval])

#pragma omp target teams distribute parallel for collapse(3)

#elif defined(_OPENMP)
#pragma omp parallel for collapse(2) default(shared)

#elif defined(_OPENACC)
#pragma acc parallel loop gang vector collapse(3)                              \
    present(ekq, dFreqGrid, inv_igp_index, indinv, aqsmtemp, aqsntemp, vcoul,  \
            I_epsR_array, I_epsA_array)                                        \
        copy(asxDtemp_re [0:nfreqeval], asxDtemp_im [0:nfreqeval])
#endif
  for (int n1 = 0; n1 < nvband; ++n1) {
    for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
      for (int iw = 0; iw < nfreqeval; ++iw) {
        double wx = freqevalmin - ekq(n1) + freqevalstep;
        int indigp = inv_igp_index(my_igp);
        int igp = indinv(indigp);
        double fact1 = 0.00, fact2 = 0.00;
        int ifreq = 0;
        ComplexType ssxDittt(0.00, 0.00);

        compute_fact(wx, nFreq, dFreqGrid, fact1, fact2, ifreq, 1, 0);

        if (wx > 0)
          ssxDittt_kernel(inv_igp_index, indinv, aqsmtemp, aqsntemp, vcoul,
                          I_epsR_array, ssxDittt, ngpown, ncouls, n1, ifreq,
                          fact1, fact2, igp, my_igp);
        else
          ssxDittt_kernel(inv_igp_index, indinv, aqsmtemp, aqsntemp, vcoul,
                          I_epsA_array, ssxDittt, ngpown, ncouls, n1, ifreq,
                          fact1, fact2, igp, my_igp);

        ssxDittt *= occ;
        double ssxDit_re = (ssxDittt).real();
        double ssxDit_im = (ssxDittt).imag();

#if defined(_OPENMP)
#pragma omp atomic
#elif defined(_OPENACC)
#pragma acc atomic update
#endif
        asxDtemp_re[iw] += ssxDit_re;
#if defined(_OPENMP)
#pragma omp atomic
#elif defined(_OPENACC)
#pragma acc atomic update
#endif
        asxDtemp_im[iw] += ssxDit_im;
      } // iw
    }
  }

#if defined(OPENMP_TARGET)
#pragma omp target exit data map(from                                          \
                                 : asxDtemp_re [0:nfreqeval],                  \
                                   asxDtemp_im [0:nfreqeval])
#endif
  for (int iw = 0; iw < nfreqeval; ++iw)
    asxDtemp(iw) = ComplexType(asxDtemp_re[iw], asxDtemp_im[iw]);

  delete[] asxDtemp_re;
  delete[] asxDtemp_im;
}

void achDtemp_cor_Kernel(int number_bands, int nvband, int nfreqeval,
                         int ncouls, int ngpown, int nFreq, double freqevalmin,
                         double freqevalstep, ARRAY1D_DataType &ekq,
                         ARRAY1D_DataType &dFreqGrid,
                         ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                         ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                         ARRAY1D_DataType &vcoul, ARRAY3D &I_epsR_array,
                         ARRAY3D &I_epsA_array, ARRAY1D &achDtemp_cor) {
  double *achDtemp_cor_re = new double[nfreqeval];
  double *achDtemp_cor_im = new double[nfreqeval];

  for (int iw = 0; iw < nfreqeval; ++iw) {
    achDtemp_cor_re[iw] = 0.0;
    achDtemp_cor_im[iw] = 0.0;
  }
  bool flag_occ;
#if defined(OPENMP_TARGET)
#pragma omp target map(tofrom                                                  \
                       : achDtemp_cor_re [0:achDtemp_cor.size],                \
                         achDtemp_cor_im [0:achDtemp_cor.size])

#pragma omp teams distribute parallel for collapse(2)

#elif defined(_OPENMP)
#pragma omp parallel for collapse(2)
#elif defined(_OPENACC)
#pragma acc enter data create(                                                 \
    achDtemp_cor_re[:achDtemp_cor.size],                                       \
                                       achDtemp_cor_im [0:achDtemp_cor.size])
#pragma acc update device(                                                     \
    achDtemp_cor_re[:achDtemp_cor.size],                                       \
                                       achDtemp_cor_im [0:achDtemp_cor.size])

#pragma acc parallel loop gang vector collapse(2)                              \
    present(dFreqGrid, inv_igp_index, indinv, I_epsR_array, I_epsA_array,      \
            vcoul, aqsmtemp, aqsntemp, ekq, achDtemp_cor_re, achDtemp_cor_im)
#endif
  for (int n1 = 0; n1 < number_bands; ++n1) {
    for (int iw = 0; iw < nfreqeval; ++iw) {
      flag_occ = n1 < nvband;
      ComplexType sch2Di(0.00, 0.00);
      ComplexType schDi_cor(0.00, 0.00);
      double wx = freqevalmin - ekq(n1) + freqevalstep;

      double fact1 = 0.00, fact2 = 0.00;
      int ifreq = 0.00;

      compute_fact(wx, nFreq, dFreqGrid, fact1, fact2, ifreq, 2, flag_occ);

      if (wx > 0) {
        if (!flag_occ)
          schDttt_corKernel1(schDi_cor, inv_igp_index, indinv, I_epsR_array,
                             I_epsA_array, aqsmtemp, aqsntemp, vcoul, ncouls,
                             ifreq, ngpown, n1, fact1, fact2);
      } else if (flag_occ)
        schDttt_corKernel2(schDi_cor, inv_igp_index, indinv, I_epsR_array,
                           I_epsA_array, aqsmtemp, aqsntemp, vcoul, ncouls,
                           ifreq, ngpown, n1, fact1, fact2);

      double schDi_re = (schDi_cor).real();
      double schDi_im = (schDi_cor).imag();

// Summing up at the end of iw loop
#if defined(_OPENMP)
#pragma omp atomic
#elif defined(_OPENACC)
#pragma acc atomic update
#endif
      achDtemp_cor_re[iw] += schDi_re;

#if defined(_OPENMP)
#pragma omp atomic
#elif defined(_OPENACC)
#pragma acc atomic update
#endif
      achDtemp_cor_im[iw] += schDi_im;
    } // iw
  }   // n1

  for (int iw = 0; iw < nfreqeval; ++iw)
    achDtemp_cor(iw) = ComplexType(achDtemp_cor_re[iw], achDtemp_cor_im[iw]);

  delete[] achDtemp_cor_re;
  delete[] achDtemp_cor_im;
}

inline void schDttt_corKernel1(ComplexType &schDttt_cor,
                               ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                               ARRAY3D &I_epsR_array, ARRAY3D &I_epsA_array,
                               ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                               ARRAY1D_DataType &vcoul, int ncouls, int ifreq,
                               int ngpown, int n1, double fact1, double fact2) {
  int blkSize = 512;
  double schDttt_cor_re = 0.00, schDttt_cor_im = 0.00, schDttt_re = 0.00,
         schDttt_im = 0.00;
#if defined(_OPENMP) && !defined(OPENMP_TARGET)
#pragma omp parallel for default(shared) collapse(2) reduction(+:schDttt_cor_re, schDttt_cor_im, schDttt_re, schDttt_im)
#endif
  for (int igbeg = 0; igbeg < ncouls; igbeg += blkSize) {
    for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
      for (int ig = igbeg; ig < min(ncouls, igbeg + blkSize); ++ig) {
        int indigp = inv_igp_index(my_igp);
        int igp = indinv(indigp);
        ComplexType sch2Dt = (I_epsR_array(ifreq, my_igp, ig) -
                              I_epsA_array(ifreq, my_igp, ig)) *
                                 fact1 +
                             (I_epsR_array((ifreq + 1), my_igp, ig) -
                              I_epsA_array((ifreq + 1), my_igp, ig)) *
                                 fact2;
        ComplexType sch2Dtt = aqsntemp(n1, ig) *
                              ComplexType_conj(aqsmtemp(n1, igp)) * sch2Dt *
                              vcoul(igp);

        schDttt_re += (sch2Dtt).real();
        schDttt_im += (sch2Dtt).imag();
        schDttt_cor_re += (sch2Dtt).real();
        schDttt_cor_im += (sch2Dtt).imag();
      }
    }
  }
  schDttt_cor = ComplexType(schDttt_cor_re, schDttt_cor_im);
}

inline void schDttt_corKernel2(ComplexType &schDttt_cor,
                               ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                               ARRAY3D &I_epsR_array, ARRAY3D &I_epsA_array,
                               ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                               ARRAY1D_DataType &vcoul, int ncouls, int ifreq,
                               int ngpown, int n1, double fact1, double fact2) {
  int blkSize = 512;
  double schDttt_cor_re = 0.00, schDttt_cor_im = 0.00;
#if defined(_OPENMP) && !defined(OPENMP_TARGET)
#pragma omp parallel for collapse(2) default(shared) reduction(+:schDttt_cor_re, schDttt_cor_im)
#endif
  for (int my_igp = 0; my_igp < ngpown; ++my_igp) {
    for (int igbeg = 0; igbeg < ncouls; igbeg += blkSize) {
      int indigp = inv_igp_index(my_igp);
      int igp = indinv(indigp);
      ComplexType sch2Dtt_store1 =
          ComplexType_conj(aqsmtemp(n1, igp)) * vcoul(igp);

      for (int ig = igbeg; ig < min(ncouls, igbeg + blkSize); ++ig) {
        ComplexType sch2Dt = ((I_epsR_array(ifreq, my_igp, ig) -
                               I_epsA_array(ifreq, my_igp, ig)) *
                                  fact1 +
                              (I_epsR_array((ifreq + 1), my_igp, ig) -
                               I_epsA_array((ifreq + 1), my_igp, ig)) *
                                  fact2) *
                             -0.5;
        ComplexType sch2Dtt = aqsntemp(n1, ig) * sch2Dt * sch2Dtt_store1;
        schDttt_cor_re += (sch2Dtt).real();
        schDttt_cor_im += (sch2Dtt).imag();
      }
    }
  }
  schDttt_cor = ComplexType(schDttt_cor_re, schDttt_cor_im);
}

int main(int argc, char **argv) {
#if defined(OPENMP_TARGET)
  cout << "\n ************OpenMP 4.5**********\n" << endl;
#elif defined(_OPENMP)
  cout << "\n ************OpenMP 3.0**********\n" << endl;
#elif defined(_OPENACC)
  cout << "\n ************OpenACC  **********\n" << endl;
#endif

  int number_bands = 0, nvband = 0, ncouls = 0, ngpown = 0, nFreq = 0,
      nfreqeval = 0;

  const double freqevalmin = 0.00;
  const double freqevalstep = 0.50;
  const double occ = 1.00;
  const double pref_zb = 0.5 / 3.14;
  double dw = -10;

  if (argc == 1) {
    number_bands = 15023;
    nvband = 1998;
    ncouls = 22401;
    ngpown = 66;
    nFreq = 15;
    nfreqeval = 10;
  } else if (argc == 7) {
    number_bands = atoi(argv[1]);
    nvband = atoi(argv[2]);
    ncouls = atoi(argv[3]);
    ngpown = atoi(argv[4]);
    nFreq = atoi(argv[5]);
    nfreqeval = atoi(argv[6]);
  } else {
    cout << "Incorrect Parameters!!! The correct form is " << endl;
    cout << "./a.out number_bands nvband ncouls ngpown nFreq nfreqeval "
         << endl;
    exit(0);
  }

  if (ngpown > ncouls) {
    cout << "Incorrect Parameters!!! ngpown cannot be greater than ncouls. The "
            "correct form is "
         << endl;
    cout << "./a.out number_bands nvband ncouls ngpown nFreq nfreqeval "
         << endl;
    exit(0);
  }

  time_point<system_clock> startTimer_Kernel, endTimer_Kernel,
      start_achsDtemp_Kernel, end_achsDtemp_Kernel, start_asxDtemp_Kernel,
      end_asxDtemp_Kernel, start_achDtemp_Kernel, end_achDtemp_Kernel,
      start_achDtemp_cor_Kernel, end_achDtemp_cor_Kernel, start_preKernel,
      end_preKernel;

  start_preKernel = system_clock::now();

  // OpenMP variables
  int tid = 0, numThreads = 1, numTeams = 1;
#if defined(_OPENMP)
#pragma omp parallel shared(numThreads) private(tid)
  {
    tid = omp_get_thread_num();
    if (tid == 0)
      numThreads = omp_get_num_threads();
  }
  std::cout << "Number of OpenMP Threads = " << numThreads << endl;
#endif

#if defined(OPENMP_TARGET)
#pragma omp target teams map(tofrom                                            \
                             : numTeams, numThreads)                           \
    shared(numTeams) private(tid)
  {
    tid = omp_get_team_num();
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
#endif

  cout << "Number of Threads = " << numThreads
       << "\n number_bands = " << number_bands << "\n nvband = " << nvband
       << "\n ncouls = " << ncouls << "\n ngpown = " << ngpown
       << "\n nFreq = " << nFreq << "\n nfreqeval = " << nfreqeval << endl;

  size_t mem_alloc = 0.00;

  // Start to allocate the data structures;
  ARRAY1D_int inv_igp_index(ngpown);
  mem_alloc += inv_igp_index.getSize();

  ARRAY1D_int indinv(ncouls);
  mem_alloc += indinv.getSize();

  ARRAY1D_DataType vcoul(ncouls);
  mem_alloc += vcoul.getSize();

  ARRAY1D_DataType ekq(number_bands);
  mem_alloc += ekq.getSize();

  ARRAY1D_DataType dFreqGrid(nFreq);
  mem_alloc += dFreqGrid.getSize();

  ARRAY1D_DataType pref(nFreq);
  mem_alloc += pref.getSize();

  ARRAY2D aqsntemp(number_bands, ncouls);
  ARRAY2D aqsmtemp(number_bands, ncouls);
  mem_alloc += 2 * aqsmtemp.getSize();

  ARRAY3D I_epsR_array(nFreq, ngpown, ncouls);
  ARRAY3D I_epsA_array(nFreq, ngpown, ncouls);
  mem_alloc += 2 * I_epsR_array.getSize();

  ARRAY1D sch2Di(nfreqeval);
  ARRAY1D schDi_cor(nfreqeval);
  ARRAY1D achDtemp_cor(nfreqeval);
  ARRAY1D asxDtemp(nfreqeval);
  mem_alloc += 4 * sch2Di.getSize();

  // Variables used :
  ComplexType achsDtemp(0.00, 0.00);

  // Initialize the data structures
  for (int ig = 0; ig < ngpown; ++ig)
    inv_igp_index(ig) = ig;

  for (int ig = 0; ig < ncouls; ++ig)
    indinv(ig) = ig;

  for (int i = 0; i < number_bands; ++i) {
    ekq(i) = dw;
    dw += 1.00;

    for (int j = 0; j < ncouls; ++j) {
      aqsmtemp(i, j) = ComplexType(0.5, 0.5);
      aqsntemp(i, j) = ComplexType(0.5, 0.5);
    }
  }

  for (int i = 0; i < ncouls; ++i)
    vcoul(i) = 1.00;

  for (int i = 0; i < nFreq; ++i) {
    for (int j = 0; j < ngpown; ++j) {
      for (int k = 0; k < ncouls; ++k) {
        I_epsR_array(i, j, k) = ComplexType(0.5, 0.5);
        I_epsA_array(i, j, k) = ComplexType(0.5, -0.5);
      }
    }
  }

  dw = 0.00;
  for (int ijk = 0; ijk < nFreq; ++ijk) {
    dFreqGrid(ijk) = dw;
    dw += 2.00;
  }

  for (int ifreq = 0; ifreq < nFreq; ++ifreq) {
    if (ifreq < nFreq - 1)
      pref(ifreq) = (dFreqGrid(ifreq + 1) - dFreqGrid(ifreq)) / 3.14;
    else
      pref(ifreq) = pref(ifreq - 1);
  }
  pref(0) *= 0.5;
  pref(nFreq - 1) *= 0.5;

  for (int i = 0; i < nfreqeval; ++i) {
    sch2Di(i) = ComplexType(0.0, 0.0);
    schDi_cor(i) = ComplexType(0.0, 0.0);
    asxDtemp(i) = ComplexType(0.0, 0.0);
    achDtemp_cor(i) = ComplexType(0.0, 0.0);
  }

  end_preKernel = system_clock::now();
  duration<double> elapsed = end_preKernel - start_preKernel;

  cout << "Memory Used = " << mem_alloc / (1024 * 1024 * 1024) << " GB" << endl;
  cout << "pre kernel time taken = " << elapsed.count() << " secs" << endl;

  cout << "starting Kernels" << endl;
  startTimer_Kernel = system_clock::now();

#if defined(OPENMP_TARGET)
#pragma omp target enter data map(to                                           \
                                  : inv_igp_index, indinv, aqsmtemp, aqsntemp, \
                                    I_epsR_array, vcoul, I_epsA_array,         \
                                    asxDtemp, dFreqGrid, ekq)

#pragma omp target enter data map(                                             \
    alloc                                                                      \
    : inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],  \
      aqsmtemp.dptr [0:aqsmtemp.size], aqsntemp.dptr [0:aqsntemp.size],        \
      I_epsR_array.dptr [0:I_epsR_array.size], vcoul.dptr [0:vcoul.size],      \
      I_epsA_array.dptr [0:I_epsA_array.size],                                 \
      asxDtemp.dptr [0:asxDtemp.size], dFreqGrid.dptr [0:dFreqGrid.size],      \
      ekq.dptr [0:ekq.size])

#pragma omp target update to(                                                  \
    inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],    \
    aqsmtemp.dptr [0:aqsmtemp.size], aqsntemp.dptr [0:aqsntemp.size],          \
    I_epsR_array.dptr [0:I_epsR_array.size], vcoul.dptr [0:vcoul.size],        \
    I_epsA_array.dptr [0:I_epsA_array.size], asxDtemp.dptr [0:asxDtemp.size],  \
    dFreqGrid.dptr [0:dFreqGrid.size], ekq.dptr [0:ekq.size])

#elif defined(_OPENACC)
#pragma acc enter data copyin(inv_igp_index, indinv, aqsmtemp, aqsntemp,       \
                              I_epsR_array, I_epsA_array, vcoul, asxDtemp,     \
                              dFreqGrid, ekq)
#pragma acc enter data copyin(                                                 \
    inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],    \
    aqsmtemp.dptr [0:aqsmtemp.size], aqsntemp.dptr [0:aqsntemp.size],          \
    I_epsR_array.dptr [0:I_epsR_array.size], vcoul.dptr [0:vcoul.size],        \
    I_epsA_array.dptr [0:I_epsA_array.size], asxDtemp.dptr [0:asxDtemp.size],  \
    dFreqGrid.dptr [0:dFreqGrid.size], ekq.dptr [0:ekq.size])

//#pragma acc enter data create(                                                 \
//    inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],    \
//    aqsmtemp.dptr [0:aqsmtemp.size], aqsntemp.dptr [0:aqsntemp.size],          \
//    I_epsR_array.dptr [0:I_epsR_array.size], vcoul.dptr [0:vcoul.size],        \
//    I_epsA_array.dptr [0:I_epsA_array.size], asxDtemp.dptr [0:asxDtemp.size],  \
//    dFreqGrid.dptr [0:dFreqGrid.size], ekq.dptr [0:ekq.size])
//
//#pragma acc update device( \
//    inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size], \
//    aqsmtemp.dptr [0:aqsmtemp.size], aqsntemp.dptr [0:aqsntemp.size], \
//    I_epsR_array.dptr [0:I_epsR_array.size], vcoul.dptr [0:vcoul.size], \
//    I_epsA_array.dptr [0:I_epsA_array.size], asxDtemp.dptr [0:asxDtemp.size],
//    \ dFreqGrid.dptr [0:dFreqGrid.size], ekq.dptr [0:ekq.size])
#endif

  endTimer_Kernel = system_clock::now();
  elapsed = endTimer_Kernel - startTimer_Kernel;
  cout << "**** Time taken to map the data on the GPU ***** " << elapsed.count()
       << endl;

  /***********achsDtemp Kernel ****************/
  start_achsDtemp_Kernel = system_clock::now();

  achsDtemp_Kernel(number_bands, ngpown, ncouls, inv_igp_index, indinv,
                   aqsntemp, aqsmtemp, I_epsR_array, vcoul, achsDtemp);

  end_achsDtemp_Kernel = system_clock::now();
  elapsed = end_achsDtemp_Kernel - start_achsDtemp_Kernel;

  cout << "********** achsDtemp Time Taken **********= " << elapsed.count()
       << " secs" << endl;

  /***********asxDtemp Kernel ****************/
  start_asxDtemp_Kernel = system_clock::now();

  asxDtemp_Kernel(nvband, nfreqeval, ncouls, ngpown, nFreq, freqevalmin,
                  freqevalstep, occ, ekq, dFreqGrid, inv_igp_index, indinv,
                  aqsmtemp, aqsntemp, vcoul, I_epsR_array, I_epsA_array,
                  asxDtemp);

  end_asxDtemp_Kernel = system_clock::now();
  elapsed = end_asxDtemp_Kernel - start_asxDtemp_Kernel;
  cout << "********** asxDtemp Time Taken **********= " << elapsed.count()
       << " secs" << endl;

  /***********achDtemp_cor Kernel ****************/
  start_achDtemp_cor_Kernel = system_clock::now();

  achDtemp_cor_Kernel(number_bands, nvband, nfreqeval, ncouls, ngpown, nFreq,
                      freqevalmin, freqevalstep, ekq, dFreqGrid, inv_igp_index,
                      indinv, aqsmtemp, aqsntemp, vcoul, I_epsR_array,
                      I_epsA_array, achDtemp_cor);

  end_achDtemp_cor_Kernel = system_clock::now();
  elapsed = end_achDtemp_cor_Kernel - start_achDtemp_cor_Kernel;
  cout << "********** achDtemp_cor Time Taken **********= " << elapsed.count()
       << " secs" << endl;

#if defined(OPENMP_TARGET)
#pragma omp target exit data map(                                              \
    delete                                                                     \
    : inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],  \
      aqsmtemp.dptr [0:aqsmtemp.size], aqsntemp.dptr [0:aqsntemp.size],        \
      I_epsR_array.dptr [0:I_epsR_array.size], vcoul.dptr [0:vcoul.size],      \
      I_epsA_array.dptr [0:I_epsA_array.size],                                 \
      asxDtemp.dptr [0:asxDtemp.size], dFreqGrid.dptr [0:dFreqGrid.size],      \
      ekq.dptr [0:ekq.size])
#endif

  endTimer_Kernel = system_clock::now();
  elapsed = endTimer_Kernel - startTimer_Kernel;

  /*******Correctness Checks*****/
  correctness(achsDtemp, asxDtemp(0), achDtemp_cor(0));

  cout << "achsDtemp = ";
  ComplexType_print(achsDtemp);
  cout << "asxDtemp = ";
  ComplexType_print(asxDtemp(0));
  cout << "achDtemp_cor = ";
  ComplexType_print(achDtemp_cor(0));

  cout
      << "********** Kernel  Time Taken (Includes data allocation) **********= "
      << elapsed.count() << " secs" << endl;

  return 0;
}
