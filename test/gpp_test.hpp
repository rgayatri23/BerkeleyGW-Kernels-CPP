
#include "Defines.h"

#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <utility>

enum GPP_Backend
{
    OPENMP_BACKEND = 0,
    OPENACC_BACKEND,
    OPENMP_TARGET_BACKEND,
};

//--------------------------------------------------------------------------------------//

void
noflagOCC_solver(size_t, size_t, size_t, ARRAY1D_int&, ARRAY1D_int&,
                 ARRAY1D_DataType& wxarray, ARRAY2D&, ARRAY2D&, ARRAY2D&, ARRAY2D&,
                 ARRAY1D_DataType&, ARRAY1D&);

//--------------------------------------------------------------------------------------//

template <size_t Idx>
inline void
init_structs(size_t number_bands, size_t ngpown, size_t ncouls, ARRAY2D& aqsmtemp,
             ARRAY2D& aqsntemp, ARRAY2D& I_eps_array, ARRAY2D& wtilde_array,
             ARRAY1D_DataType& vcoul, ARRAY1D_int& inv_igp_index, ARRAY1D_int& indinv,
             ARRAY1D_DataType& wx_array)
{
    // Constants that will be used later
    const DataType e_lk   = 10;
    const DataType dw     = 1;
    const DataType to1    = 1e-6;
    const DataType e_n1kq = 6.0;

    ComplexType expr(.5, .5);
    for(int i = 0; i < number_bands; i++)
        for(int j = 0; j < ncouls; j++)
        {
            aqsmtemp(i, j) = expr;
            aqsntemp(i, j) = expr;
        }

    for(int i = 0; i < ngpown; i++)
        for(int j = 0; j < ncouls; j++)
        {
            I_eps_array(i, j)  = expr;
            wtilde_array(i, j) = expr;
        }

    for(int i = 0; i < ncouls; i++)
        vcoul(i) = 1.0;

    for(int ig = 0; ig < ngpown; ++ig)
        inv_igp_index(ig) = (ig + 1) * ncouls / ngpown;

    // Do not know yet what this array represents
    for(int ig = 0; ig < ncouls; ++ig)
        indinv(ig) = ig;
    indinv(ncouls) = ncouls - 1;

    for(int iw = nstart; iw < nend; ++iw)
    {
        wx_array(iw) = e_lk - e_n1kq + dw * ((iw + 1) - 2);
        if(wx_array(iw) < to1)
            wx_array(iw) = to1;
    }
}

//--------------------------------------------------------------------------------------//

template <>
inline void
init_structs<OPENACC_BACKEND>(size_t number_bands, size_t ngpown, size_t ncouls,
                              ARRAY2D& aqsmtemp, ARRAY2D& aqsntemp, ARRAY2D& I_eps_array,
                              ARRAY2D& wtilde_array, ARRAY1D_DataType& vcoul,
                              ARRAY1D_int& inv_igp_index, ARRAY1D_int& indinv,
                              ARRAY1D_DataType& wx_array)
{
    const DataType dw   = 1;
    const DataType e_lk = 10;
    const DataType to1  = 1e-6;
    // const DataType limittwo = pow(0.5, 2);
    const DataType e_n1kq = 6.0;
    ComplexType    expr(0.5, 0.5);

#pragma acc enter data copyin(aqsmtemp, aqsntemp, vcoul, inv_igp_index, indinv,          \
                              I_eps_array, wtilde_array, wx_array)

#pragma acc enter data create(                                                           \
    aqsmtemp.dptr [0:aqsmtemp.size], vcoul.dptr [0:vcoul.size],                          \
    inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],              \
    aqsntemp.dptr [0:aqsntemp.size], I_eps_array.dptr [0:I_eps_array.size],              \
    wx_array.dptr [nstart:nend], wtilde_array.dptr [0:wtilde_array.size])

#pragma acc parallel loop present(aqsmtemp, aqsntemp)
    for(int i = 0; i < number_bands; i++)
        for(int j = 0; j < ncouls; j++)
        {
            aqsmtemp(i, j) = ComplexType(0.5, 0.5);
            aqsntemp(i, j) = ComplexType(0.5, 0.5);
        }

#pragma acc parallel loop copyin(expr) present(I_eps_array, wtilde_array)
    for(int i = 0; i < ngpown; i++)
        for(int j = 0; j < ncouls; j++)
        {
            I_eps_array(i, j)  = expr;
            wtilde_array(i, j) = expr;
        }

#pragma acc parallel loop present(vcoul)
    for(int i = 0; i < ncouls; i++)
        vcoul(i) = 1.0;

#pragma acc parallel loop present(inv_igp_index)
    for(int ig = 0; ig < ngpown; ++ig)
        inv_igp_index(ig) = (ig + 1) * ncouls / ngpown;

#pragma acc parallel loop present(indinv)
    for(int ig = 0; ig < ncouls; ++ig)
        indinv(ig) = ig;
    indinv(ncouls) = ncouls - 1;

#pragma acc parallel loop present(wx_array)
    for(int iw = nstart; iw < nend; ++iw)
    {
        wx_array(iw) = e_lk - e_n1kq + dw * ((iw + 1) - 2);
        if(wx_array(iw) < to1)
            wx_array(iw) = to1;
    }
}

//--------------------------------------------------------------------------------------//

template <size_t Idx>
auto
get_num_threads_and_teams()
{
    return std::tuple<int, int>{ 1, 2 };
}

//--------------------------------------------------------------------------------------//

template <>
auto
get_num_threads_and_teams<OPENMP_BACKEND>()
{
    int numThreads = 1;
    int numTeams   = 2;
#if defined(_OPENMP)
#    pragma omp parallel shared(numThreads)
    {
        int tid = omp_get_thread_num();
        if(tid == 0)
            numThreads = omp_get_num_threads();
    }
#endif
    return std::make_tuple(numThreads, numTeams);
}

//--------------------------------------------------------------------------------------//

template <>
auto
get_num_threads_and_teams<OPENMP_TARGET_BACKEND>()
{
    int numThreads = 1;
    int numTeams   = 2;
#if defined(GPP_OPENMP_TARGET)
#    pragma omp target teams map(tofrom : numTeams, numThreads) shared(numTeams)
    {
        int tid = omp_get_team_num();
        if(tid == 0)
        {
            numTeams = omp_get_num_teams();
#    pragma omp parallel
            {
                int ttid = omp_get_thread_num();
                if(ttid == 0)
                    numThreads = omp_get_num_threads();
            }
        }
    }
#endif
    return std::make_tuple(numThreads, numTeams);
}

//--------------------------------------------------------------------------------------//

template <size_t Idx>
auto
run(int number_bands = 512, int nvband = 2, int ncouls = 32768, int nodes_per_group = 20,
    int npes = 1)
{
    int ngpown = ncouls / (nodes_per_group * npes);

    // Using time point and system_clock
    // time_point<system_clock> start, end, k_start, k_end;
    // start = system_clock::now();
    // double elapsedKernelTimer;

    int numThreads                 = 1;
    int numTeams                   = 2;
    std::tie(numThreads, numTeams) = get_num_threads_and_teams<Idx>();

    // Printing out the params passed.
    std::cout << "Sizeof(ComplexType = " << sizeof(ComplexType) << " bytes" << std::endl;
    std::cout << "number_bands = " << number_bands << "\t nvband = " << nvband
              << "\t ncouls = " << ncouls << "\t nodes_per_group  = " << nodes_per_group
              << "\t ngpown = " << ngpown << "\t nend = " << nend
              << "\t nstart = " << nstart << endl;

    // size_t memFootPrint = 0.00;

    // ALLOCATE statements from fortran gppkernel.
    ARRAY1D achtemp(nend - nstart);
    // memFootPrint += (nend - nstart) * sizeof(ComplexType);

    ARRAY2D aqsmtemp(number_bands, ncouls);
    ARRAY2D aqsntemp(number_bands, ncouls);
    // memFootPrint += 2 * (number_bands * ncouls) * sizeof(ComplexType);

    ARRAY2D I_eps_array(ngpown, ncouls);
    ARRAY2D wtilde_array(ngpown, ncouls);
    // memFootPrint += 2 * (ngpown * ncouls) * sizeof(ComplexType);

    ARRAY1D_DataType vcoul(ncouls);
    // memFootPrint += ncouls * sizeof(DataType);

    ARRAY1D_int inv_igp_index(ngpown);
    ARRAY1D_int indinv(ncouls + 1);
    // memFootPrint += ngpown * sizeof(int);
    // memFootPrint += (ncouls + 1) * sizeof(int);

    ARRAY1D_DataType wx_array(nend - nstart);
    // memFootPrint += 3 * (nend - nstart) * sizeof(DataType);

    // Print Memory Foot print
    // cout << "Memory Foot Print = " << memFootPrint / pow(1024, 3) << " GBs"
    //     << endl;

    init_structs<Idx>(number_bands, ngpown, ncouls, aqsmtemp, aqsntemp, I_eps_array,
                      wtilde_array, vcoul, inv_igp_index, indinv, wx_array);

    // k_start = system_clock::now();
    noflagOCC_solver(number_bands, ngpown, ncouls, inv_igp_index, indinv, wx_array,
                     wtilde_array, aqsmtemp, aqsntemp, I_eps_array, vcoul, achtemp);

    // k_end = system_clock::now();
    // duration<double> elapsed = k_end - k_start;
    // elapsedKernelTimer = elapsed.count();

    // end = system_clock::now();
    // elapsed = end - start;

    // cout << "********** Kernel Time Taken **********= " << elapsedKernelTimer
    //     << " secs" << endl;
    // cout << "********** Total Time Taken **********= " << elapsed.count()
    //     << " secs" << endl;

    return achtemp(0);
}

//--------------------------------------------------------------------------------------//

inline auto
occ_solver_impl(size_t n1, size_t my_igp, size_t ig, ARRAY1D_int& inv_igp_index,
                ARRAY1D_int& indinv, ARRAY1D_DataType& wx_array, ARRAY2D& wtilde_array,
                ARRAY2D& aqsmtemp, ARRAY2D& aqsntemp, ARRAY2D& I_eps_array,
                ARRAY1D_DataType& vcoul)
{
    int      indigp = inv_igp_index(my_igp);
    int      igp    = indinv(indigp);
    DataType achtemp_re_loc[nend - nstart];
    DataType achtemp_im_loc[nend - nstart];
    for(int iw = nstart; iw < nend; ++iw)
    {
        achtemp_re_loc[iw] = 0.00;
        achtemp_im_loc[iw] = 0.00;
    }
    ComplexType sch_store1 = ComplexType_conj(aqsmtemp(n1, igp)) * aqsntemp(n1, igp) *
                             0.5 * vcoul(igp) * wtilde_array(my_igp, igp);

    for(int iw = nstart; iw < nend; ++iw)
    {
        ComplexType wdiff = wx_array(iw) - wtilde_array(my_igp, ig);
        ComplexType delw =
            ComplexType_conj(wdiff) * (1 / (wdiff * ComplexType_conj(wdiff)).real());
        ComplexType sch_array = delw * I_eps_array(my_igp, ig) * sch_store1;

        achtemp_re_loc[iw] += (sch_array).real();
        achtemp_im_loc[iw] += (sch_array).imag();
    }

    return std::make_tuple(achtemp_re_loc[0], achtemp_re_loc[1], achtemp_re_loc[2],
                           achtemp_im_loc[0], achtemp_im_loc[1], achtemp_im_loc[2]);
}

//--------------------------------------------------------------------------------------//

#if !defined(FOLD_EXPRESSION)
#    define FOLD_EXPRESSION(...)                                                         \
        (void) ::std::initializer_list<int> { (__VA_ARGS__, 0)... }
#endif

//--------------------------------------------------------------------------------------//

template <typename LhsT, typename RhsT, size_t... Idx>
decltype(auto)
plus_eq(LhsT&& _lhs, RhsT&& _rhs, index_sequence<Idx...>)
{
    FOLD_EXPRESSION(std::get<Idx>(_lhs) += std::get<Idx>(_rhs));
    return std::forward<LhsT>(_lhs);
}

//--------------------------------------------------------------------------------------//

template <typename LhsT, typename RhsT>
decltype(auto)
plus_eq(LhsT&& _lhs, RhsT&& _rhs)
{
    constexpr auto N = std::tuple_size<LhsT>::value;
    return plus_eq(std::forward<LhsT>(_lhs), std::forward<RhsT>(_rhs),
                   std::make_index_sequence<N>{});
}

//--------------------------------------------------------------------------------------//

void
noflagOCC_solver(size_t number_bands, size_t ngpown, size_t ncouls,
                 ARRAY1D_int& inv_igp_index, ARRAY1D_int& indinv,
                 ARRAY1D_DataType& wx_array, ARRAY2D& wtilde_array, ARRAY2D& aqsmtemp,
                 ARRAY2D& aqsntemp, ARRAY2D& I_eps_array, ARRAY1D_DataType& vcoul,
                 ARRAY1D& achtemp)
{
    // Vars to use for reduction
    DataType ach_re0 = 0.00, ach_re1 = 0.00, ach_re2 = 0.00, ach_im0 = 0.00,
             ach_im1 = 0.00, ach_im2 = 0.00;

#if defined(OPENMP_TARGET)
#    pragma omp target enter data map(to                                                 \
                                      : aqsmtemp, vcoul, inv_igp_index, indinv,          \
                                        aqsntemp, I_eps_array, wx_array, wtilde_array)

    // start = system_clock::now();

#    pragma omp target map(                                                              \
        to                                                                               \
        : aqsmtemp.dptr [0:aqsmtemp.size], vcoul.dptr [0:vcoul.size],                    \
          inv_igp_index.dptr [0:inv_igp_index.size], indinv.dptr [0:indinv.size],        \
          aqsntemp.dptr [0:aqsmtemp.size], I_eps_array.dptr [0:I_eps_array.size],        \
          wx_array.dptr [nstart:nend], wtilde_array.dptr [0:wtilde_array.size])          \
        map(tofrom                                                                       \
            : ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)

#pragma omp teams distribute parallel for collapse(2) \
    reduction(+:ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2) \
  num_teams(ngpown*ncouls) thread_limit(32)
    for(int my_igp = 0; my_igp < ngpown; ++my_igp)
    {
        for(int ig = 0; ig < ncouls; ++ig)
        {
            for(int n1 = 0; n1 < number_bands; ++n1)
            {
#elif defined(_OPENMP)
#pragma omp parallel for collapse(2) \
    reduction(+:ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)
    for(int n1 = 0; n1 < number_bands; ++n1)
    {
        for(int my_igp = 0; my_igp < ngpown; ++my_igp)
        {
            for(int ig = 0; ig < ncouls; ++ig)
            {
#elif defined(_OPENACC)
#pragma acc parallel loop gang vector collapse(3) \
    present(inv_igp_index, indinv, aqsmtemp, aqsntemp, wtilde_array, wx_array, I_eps_array, vcoul) \
    reduction(+:ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)\
    num_gangs(number_bands*ncouls)
    for(int n1 = 0; n1 < number_bands; ++n1)
    {
        for(int my_igp = 0; my_igp < ngpown; ++my_igp)
        {
            for(int ig = 0; ig < ncouls; ++ig)
            {
#else
    for(int n1 = 0; n1 < number_bands; ++n1)
    {
        for(int my_igp = 0; my_igp < ngpown; ++my_igp)
        {
            for(int ig = 0; ig < ncouls; ++ig)
            {
#endif
                plus_eq(std::forward_as_tuple(ach_re0, ach_re1, ach_im2, ach_im0, ach_im1,
                                              ach_im2),
                        occ_solver_impl(n1, my_igp, ig, inv_igp_index, indinv, wx_array,
                                        wtilde_array, aqsmtemp, aqsmtemp, I_eps_array,
                                        vcoul));
            }
        }  // ngpown
    }      // number_bands

    // end                                 = system_clock::now();
    // duration<double> elapsed            = end - start;
    // double           elapsedKernelTimer = elapsed.count();

    achtemp(0) = ComplexType(ach_re0, ach_im0);
    achtemp(1) = ComplexType(ach_re1, ach_im1);
    achtemp(2) = ComplexType(ach_re2, ach_im2);
}
