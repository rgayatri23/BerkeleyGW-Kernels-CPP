/*
 * Defines a common set of headers and typdefs used in both GPP and FF mini-apps
 * */
#include "../ComplexClass/CustomComplex.h"
#include "../arrayMD/arrayMDcpu.h"
#include "../arrayMD/arrayMD.h"
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

// ArrayMD definitions
#define ARRAY3D ArrayMD<CustomComplex<DataType>, 3, Device::cpu>
#define ARRAY2D ArrayMD<CustomComplex<DataType>, 2, Device::cpu>
#define ARRAY1D ArrayMD<CustomComplex<DataType>, 1, Device::cpu>
#define ARRAY1D_int ArrayMD<int, 1, Device::cpu>
#define ARRAY1D_DataType ArrayMD<DataType, 1, Device::cpu>

// Function Definitions

void schDttt_corKernel1(CustomComplex<double> &schDttt_cor,
                        ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                        ARRAY3D &I_epsR_array, ARRAY3D &I_epsA_array,
                        ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                        ARRAY1D_DataType &vcoul, int ncouls, int ifreq,
                        int ngpown, int n1, double fact1, double fact2);

void schDttt_corKernel2(CustomComplex<double> &schDttt_cor,
                        ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
                        ARRAY3D &I_epsR_array, ARRAY3D &I_epsA_array,
                        ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
                        ARRAY1D_DataType &vcoul, int ncouls, int ifreq,
                        int ngpown, int n1, double fact1, double fact2);

void
noflagOCC_solver(size_t number_bands,
                 size_t ngpown,
                 size_t ncouls,
                 ARRAY1D_int& inv_igp_index,
                 ARRAY1D_int& indinv,
                 ARRAY1D_DataType& wx_array,
                 ARRAY2D& wtilde_array,
                 ARRAY2D& aqsmtemp,
                 ARRAY2D& aqsntemp,
                 ARRAY2D& I_eps_array,
                 ARRAY1D_DataType& vcoul,
                 ARRAY1D& achtemp);
