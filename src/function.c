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

#define D_2DPI 0.63661977236758134307553505349006 /* 2/PI */

double gtoint__compute_cartesian_normalization_constant(int3_t a, double g) {
    double v = D_2DPI * g;
    v *= sqrt(v);
    for (int i = 0; i < a.x; i++) v *= g / (0.5 * i + 0.25);
    for (int i = 0; i < a.y; i++) v *= g / (0.5 * i + 0.25);
    for (int i = 0; i < a.z; i++) v *= g / (0.5 * i + 0.25);
    return sqrt(v);
}

double gtoint__compute_spherical_harmonics_normalization_constant(int a, double g) {
    double v = D_2DPI * g;
    v *= sqrt(v);
    for (int i = 0; i < a; i++) v *= g / (0.5 * i + 0.25);
    return sqrt(v);
}
