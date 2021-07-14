/*
 * LibGtoint: an analytical GTO integral library for C and Fortran.
 *
 * Copyright (c) 2020-2021 Arihiro Yoshida. All rights reserved.
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

#include "function.h"

#include <math.h>

#define D_RPID2 0.88622692545275801364908374167057 /* sqrt(PI)/2 */

#define BOYS_NM 100
#define BOYS_NX 300
#define BOYS_DX 0.125
#define BOYS_AX 8.0 /* 1.0 / BOYS_DX */

static double g_boys_table_[BOYS_NM][BOYS_NX];

static void compute_boys_function_table_(void) {
    {
        const int m = BOYS_NM - 1;
        for (int i = 0; i < BOYS_NX; i++) {
            const double x = BOYS_DX * i;
            g_boys_table_[m][i] = exp(-x) / (2 * m + 1);
        }
    }
    for (int m = BOYS_NM - 2; m >= 0; m--) {
        for (int i = 0; i < BOYS_NX; i++) {
            const double x = BOYS_DX * i;
            g_boys_table_[m][i] = (2.0 * x * g_boys_table_[m + 1][i] + exp(-x)) / (2 * m + 1);
        }
    }
}

double gtoint__compute_boys_function(int m, double x, double tol) {
    if (x < 0.0) x = 0.0; /* fail-safe */
    if (x < BOYS_DX * (BOYS_NX - 0.5)) {
        const int ix = (int)floor(x * BOYS_AX + 0.5);
        const double x0 = BOYS_DX * ix;
        const double dx = x - x0;
        double f = 1.0;
        double y = g_boys_table_[m][ix];
        for (int i = 1;; i++) {
            f *= -dx / i;
            const double d = f * ((m + i < BOYS_NM) ? g_boys_table_[m + i][ix] : 0.0);
            y += d;
            if (fabs(d) <= fabs(y) * tol) break;
        }
        return y;
    }
    else {
        const double r = 0.5 / x;
        double y = D_RPID2 / sqrt(x);
        for (int i = 0; i < m; i++) {
            y *= (2 * i + 1) * r;
        }
        return y;
    }
}

#undef BOYS_NM
#undef BOYS_NX
#undef BOYS_DX
#undef BOYS_AX

#ifdef _MSC_VER
#pragma section(".CRT$XCU",read)
#else
__attribute__((constructor))
#endif
static void construct(void) {
    compute_boys_function_table_();
}

#ifdef _MSC_VER
__declspec(allocate(".CRT$XCU")) void (*construct_function_boys_)(void) = construct;
#ifdef _WIN64
__pragma(comment(linker,"/include:construct_function_boys_"))
#else
__pragma(comment(linker,"/include:_construct_function_boys_"))
#endif
#endif
