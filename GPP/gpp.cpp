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
#include <Kokkos_Core.hpp>

inline void correntess(KokkosComplex result) {
  double re_diff, im_diff;
  re_diff = result.real() - -264241156.144946;
  im_diff = result.imag() - 1321205754.194689;

  if (re_diff < 0.001 && im_diff < 0.001)
    printf("\n!!!! SUCCESS - !!!! Correctness test passed :-D :-D\n\n");
  else
    printf("\n!!!! FAILURE - Correctness test failed :-( :-(.\n You are away "
           "by (%f,%f) \n",
           re_diff, im_diff);
}

int main(int argc, char **argv) {

  Kokkos::initialize(argc, argv);
  {
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
      std::cout << "The correct form of input is : "
                << "\n";
      std::cout << " ./a.out <number_bands> <number_valence_bands> "
                   "<number_plane_waves> <nodes_per_mpi_group> "
                << "\n";
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

    // Printing out the params passed.
    std::cout << "Sizeof(KokkosComplex) = " << sizeof(KokkosComplex) << " bytes"
              << "\n";
    std::cout << "number_bands = " << number_bands << "\t nvband = " << nvband
              << "\t ncouls = " << ncouls
              << "\t nodes_per_group  = " << nodes_per_group
              << "\t ngpown = " << ngpown << "\t nend = " << nend
              << "\t nstart = " << nstart << "\n";

    size_t memFootPrint = 0.00;

    // ALLOCATE statements from fortran gppkernel.
    ComplexStruct achtemp;
    memFootPrint += (nend - nstart) * sizeof(KokkosComplex);

    ARRAY2D aqsmtemp("aqsmtemp", number_bands, ncouls);
    ARRAY2D aqsntemp("aqsntemp", number_bands, ncouls);
    memFootPrint += 2 * (number_bands * ncouls) * sizeof(KokkosComplex);

    ARRAY2D I_eps_array("I_eps_array", ngpown, ncouls);
    ARRAY2D wtilde_array("wtilde_array", ngpown, ncouls);
    memFootPrint += 2 * (ngpown * ncouls) * sizeof(KokkosComplex);

    ARRAY1D_DataType vcoul("vcoul", ncouls);
    memFootPrint += ncouls * sizeof(DataType);

    ARRAY1D_int inv_igp_index("inv_igp_index", ngpown);
    ARRAY1D_int indinv("indinv", ncouls + 1);

    memFootPrint += ngpown * sizeof(int);
    memFootPrint += (ncouls + 1) * sizeof(int);

    ARRAY1D_DataType wx_array("wx_array", nend - nstart);
    memFootPrint += 3 * (nend - nstart) * sizeof(DataType);

    // Print Memory Foot print
    std::cout << "Memory Foot Print = " << memFootPrint / pow(1024, 3) << " GBs"
              << "\n";

    KokkosComplex expr(.5, .5);
    for (int i = 0; i < number_bands; i++)
      for (int j = 0; j < ncouls; j++) {
        aqsmtemp.h_view(i, j) = expr;
        aqsntemp.h_view(i, j) = expr;
      }

    for (int i = 0; i < ngpown; i++)
      for (int j = 0; j < ncouls; j++) {
        I_eps_array.h_view(i, j) = expr;
        wtilde_array.h_view(i, j) = expr;
      }

    for (int i = 0; i < ncouls; i++)
      vcoul.h_view(i) = 1.0;

    for (int ig = 0; ig < ngpown; ++ig)
      inv_igp_index.h_view(ig) = (ig + 1) * ncouls / ngpown;

    // Do not know yet what this array represents
    for (int ig = 0; ig < ncouls; ++ig)
      indinv.h_view(ig) = ig;
    indinv.h_view(ncouls) = ncouls - 1;

    for (int iw = nstart; iw < nend; ++iw) {
      wx_array.h_view(iw) = e_lk - e_n1kq + dw * ((iw + 1) - 2);
      if (wx_array.h_view(iw) < to1)
        wx_array.h_view(iw) = to1;
    }

    k_start = system_clock::now();
    noflagOCC_solver(number_bands, ngpown, ncouls, inv_igp_index, indinv,
                     wx_array, wtilde_array, aqsmtemp, aqsntemp, I_eps_array,
                     vcoul, achtemp);

    k_end = system_clock::now();
    duration<double> elapsed = k_end - k_start;
    elapsedKernelTimer = elapsed.count();

    // Check for correctness
    correntess(achtemp.val[0]);
    printf("\n Final achtemp\n");
    printf("(%f,%f) \n", achtemp.val[0].real(), achtemp.val[0].imag());

    end = system_clock::now();
    elapsed = end - start;

    std::cout << "********** Kernel Time Taken **********= "
              << elapsedKernelTimer << " secs"
              << "\n";
    std::cout << "********** Total Time Taken **********= " << elapsed.count()
              << " secs"
              << "\n";
  }
  Kokkos::finalize();
  return 0;
}

void noflagOCC_solver(size_t number_bands, size_t ngpown, size_t ncouls,
                      ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                      ARRAY1D_DataType &wx_array, ARRAY2D &wtilde_array,
                      ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                      ARRAY2D &I_eps_array, ARRAY1D_DataType &vcoul,
                      ComplexStruct &achtemp) {
  time_point<system_clock> start, end;
  start = system_clock::now();

  // Deep copies between host views to device views of the DualViews.
  Kokkos::deep_copy(inv_igp_index.d_view, inv_igp_index.h_view);
  Kokkos::deep_copy(indinv.d_view, indinv.h_view);
  Kokkos::deep_copy(vcoul.d_view, vcoul.h_view);
  Kokkos::deep_copy(wx_array.d_view, wx_array.h_view);
  Kokkos::deep_copy(aqsmtemp.d_view, aqsmtemp.h_view);
  Kokkos::deep_copy(aqsntemp.d_view, aqsntemp.h_view);
  Kokkos::deep_copy(I_eps_array.d_view, I_eps_array.h_view);
  Kokkos::deep_copy(wtilde_array.d_view, wtilde_array.h_view);

  // An MDRangePolicy for 3 nested loops.
  using md_policy =
      Kokkos::Experimental::MDRangePolicy<ExecSpace,
                                          Kokkos::Experimental::Rank<3>, int>;
#if defined(KOKKOS_ENABLE_CUDA)
  md_policy iteration({0, 0, 0}, {ncouls, number_bands, ngpown}, {8, 8, 8});
#else
  md_policy iteration({0, 0, 0}, {number_bands, ngpown, ncouls}, {0, 0, 0});
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  Kokkos::parallel_reduce(
      "noflagOCC_solver", iteration,
      KOKKOS_LAMBDA(const int ig, const int n1, const int my_igp,
                    ComplexStruct &update) {
#else
  Kokkos::parallel_reduce(
      "noflagOCC_solver", iteration,
      KOKKOS_LAMBDA(const int n1, const int my_igp, const int ig,
                    ComplexStruct &update) {
#endif
        int indigp = inv_igp_index.d_view(my_igp);
        int igp = indinv.d_view(indigp);
        KokkosComplex aqsmtemp_conj = KokkosComplex(
            aqsmtemp.d_view(n1, igp).real(), -aqsmtemp.d_view(n1, igp).imag());
        KokkosComplex sch_store1 = aqsmtemp_conj * aqsntemp.d_view(n1, igp) *
                                   0.5 * vcoul.d_view(igp) *
                                   wtilde_array.d_view(my_igp, igp);

        for (int iw = nstart; iw < nend; ++iw) {
          KokkosComplex wdiff =
              wx_array.d_view(iw) - wtilde_array.d_view(my_igp, ig);
          KokkosComplex wdiff_conj = KokkosComplex(wdiff.real(), -wdiff.imag());
          KokkosComplex delw = wdiff_conj * (1 / (wdiff * wdiff_conj).real());
          KokkosComplex sch_array =
              delw * I_eps_array.d_view(my_igp, ig) * sch_store1;

          update.val[iw].real() += sch_array.real();
          update.val[iw].imag() += sch_array.imag();
        }
      },
      achtemp); // number_bands

  end = system_clock::now();
  duration<double> elapsed = end - start;
  double elapsedKernelTimer = elapsed.count();
}
