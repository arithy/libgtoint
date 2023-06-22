/*
 * LibGtoint: an analytical GTO integral library for C and Fortran.
 *
 * Copyright (c) 2020-2023 Arihiro Yoshida. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#ifndef NDEBUG
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif
#endif

#include "gtoint-private.h"

#include <stdio.h>
#include <math.h>

int main(int argc, char **argv) {
#ifdef _MSC_VER
#ifndef NDEBUG
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
#endif
#endif
    int r = 0;
    gtoint_error_t e = GTOINT_ERROR_OK;
    gtoint_integrator_t itg = 0;
    {
        if ((e = gtoint_integrator_create(&itg)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_integrator_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
        int c = 0;
        {
            const gtoint_double3_t p0 = { -0.94448511705458826, 2.8224949396413224, -0.30046645380488102 };
            const gtoint_int3_t m0[] = { { 0, 0, 0 } };
            const double e0[] = { 130.70932139999999, 23.808866049999999, 6.4436083130000004 };
            const double c0[] = { 4.2519432778903488, 4.1122944243038289, 1.2816225514654456 };
            const gtoint_double3_t p1 = { -0.94448511705458826, 2.8224949396413224, -0.30046645380488102 };
            const gtoint_int3_t m1[] = { { 0, 0, 0 } };
            const double e1[] = { 130.70932139999999, 23.808866049999999, 6.4436083130000004 };
            const double c1[] = { 4.2519432778903488, 4.1122944243038289, 1.2816225514654456 };
            const gtoint_double3_t pc = { 0.17687836525872241, -0.00037794522491180004, -0.013417055484368902 };
            const double gc[] = { 11.793279999999999 };
            const double cc[] = { 6.2076000000000002 };
            const gtoint_int3_t d[] = { { 0, 0, 0 } };
            if ((e = gtoint__compute_scalar_ecp_type2_integrals(
                itg,
                &p0, 1, m0, 3, e0, c0,
                &p1, 1, m1, 3, e1, c1,
                &pc, 0, 0, 1, gc, cc,
                1, d, d, d
            )) != GTOINT_ERROR_OK) {
                printf("ERROR: %d: gtoint__compute_scalar_ecp_type2_integrals() -> %d\n", __LINE__, e);
                goto EXIT;
            }
            for (size_t i1 = 0; i1 < 3; i1++) {
            for (size_t i0 = 0; i0 < 3; i0++) {
                const double v = itg->v.p[i0 + 3 * i1];
                if (isnan(v) || isinf(v)) r = 1;
            }
            }
            printf("Case %d, Result: %s\n", c, r ? "FAIL" : "PASS");
        }
        c++;
    }
EXIT:;
    gtoint_integrator_destroy(itg);
    return (r == 0 && e == GTOINT_ERROR_OK) ? 0 : 1;
}
