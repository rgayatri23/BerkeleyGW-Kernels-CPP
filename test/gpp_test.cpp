// MIT License
//
// Copyright (c) 2020, The Regents of the University of California,
// through Lawrence Berkeley National Laboratory (subject to receipt of any
// required approvals from the U.S. Dept. of Energy).  All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "gtest/gtest.h"

#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include "gpp_test.hpp"

static int    _argc = 0;
static char** _argv = nullptr;

using mutex_t = std::mutex;
using lock_t  = std::unique_lock<mutex_t>;

//--------------------------------------------------------------------------------------//

namespace details
{
//  Get the current tests name
inline std::string
get_test_name()
{
    return ::testing::UnitTest::GetInstance()->current_test_info()->name();
}
}  // namespace details

//--------------------------------------------------------------------------------------//

class GPP_FIXTURE_NAME : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

//--------------------------------------------------------------------------------------//

#if defined(GPP_USE_OPENMP_TARGET)
static constexpr auto BACKEND_TYPE = OPENMP_TARGET_BACKEND;
#elif defined(GPP_USE_OPENMP)
static constexpr auto BACKEND_TYPE = OPENMP_BACKEND;
#else
static constexpr auto BACKEND_TYPE = OPENACC_BACKEND;
#endif

//--------------------------------------------------------------------------------------//

TEST_F(GPP_FIXTURE_NAME, standard)
{
    auto result = run<BACKEND_TYPE>();

    printf("\n %s :: Final achtemp\n", details::get_test_name().c_str());
    ComplexType_print(result);

#if defined(_OPENMP) || defined(_OPENACC)
    EXPECT_NEAR(result.real(), -264241149.849658, 0.001);
    EXPECT_NEAR(result.imag(), 1321205773.349384, 0.001);
#else
    EXPECT_NEAR(result.real(), -264241220.914570, 0.001);
    EXPECT_NEAR(result.imag(), 1321205332.084101, 0.001);
#endif
}

//--------------------------------------------------------------------------------------//

int
main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    _argc = argc;
    _argv = argv;

    auto ret = RUN_ALL_TESTS();

    // tim::timemory_finalize();
    // tim::dmp::finalize();
    return ret;
}

//--------------------------------------------------------------------------------------//
