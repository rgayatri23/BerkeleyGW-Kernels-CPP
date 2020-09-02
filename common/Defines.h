/*
 * Defines a common set of headers and typdefs used in both GPP and FF mini-apps
 * */
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <chrono>
#include <thread>

#if _OPENMP
#include <omp.h>
#endif
#if _OPENACC
#include <openacc.h>
#endif

#define nstart 0
#define nend 3

using DataType = double;
using namespace std::chrono;

using KokkosComplex = Kokkos::complex<DataType>;
using ExecSpace = Kokkos::DefaultExecutionSpace;
using HostExecSpace = Kokkos::DefaultHostExecutionSpace;
using Layout = ExecSpace::array_layout;
using MemSpace = ExecSpace::memory_space;

// ArrayMD definitions
#define ARRAY3D Kokkos::View<KokkosComplex ***, Layout, MemSpace>
#define ARRAY2D Kokkos::DualView<KokkosComplex **, Layout, MemSpace>
#define ARRAY1D Kokkos::DualView<KokkosComplex *, Layout, MemSpace>
#define ARRAY1D_DataType Kokkos::DualView<DataType *, Layout, MemSpace>
#define ARRAY1D_int Kokkos::DualView<int *, Layout, MemSpace>

// Structure to store nend-nstart number of complex-numbers and define a
// reduction on those elements.
struct ComplexStruct {
  KokkosComplex val[nend - nstart];

  KOKKOS_INLINE_FUNCTION
  void operator+=(ComplexStruct const &other) {
    for (int i = 0; i < 3; ++i)
      val[i] += other.val[i];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(ComplexStruct const volatile &other) volatile {
    for (int i = 0; i < 3; ++i)
      val[i] += other.val[i];
  }
};

// Function Definitions

void schDttt_corKernel1(KokkosComplex &schDttt_cor, ARRAY1D_int &inv_igp_index,
                        ARRAY1D_int &indinv, ARRAY3D &I_epsR_array,
                        ARRAY3D &I_epsA_array, ARRAY2D &aqsmtemp,
                        ARRAY2D &aqsntemp, ARRAY1D_DataType &vcoul, int ncouls,
                        int ifreq, int ngpown, int n1, double fact1,
                        double fact2);

void schDttt_corKernel2(KokkosComplex &schDttt_cor, ARRAY1D_int &inv_igp_index,
                        ARRAY1D_int &indinv, ARRAY3D &I_epsR_array,
                        ARRAY3D &I_epsA_array, ARRAY2D &aqsmtemp,
                        ARRAY2D &aqsntemp, ARRAY1D_DataType &vcoul, int ncouls,
                        int ifreq, int ngpown, int n1, double fact1,
                        double fact2);

void noflagOCC_solver(size_t number_bands, size_t ngpown, size_t ncouls,
                      ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                      ARRAY1D_DataType &wx_array, ARRAY2D &wtilde_array,
                      ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                      ARRAY2D &I_eps_array, ARRAY1D_DataType &vcoul,
                      ComplexStruct &achtemp);
