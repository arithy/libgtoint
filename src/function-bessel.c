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

inline static double power_factorial_half_shifted_(double a, int n) {
    double v = 0.5;
    for (int i = 1; i <= n; i++) v *= a * (i + 0.5);
    return v;
}

double gtoint__compute_modified_spherical_bessel_function_first_kind(int n, double x, double e) { /* The relative error is about 1e-6. */
    if (n >= 0) {
        if (x == 0.0) return (n > 0) ? 0.0 : (e == 0.0) ? 1.0 : exp(e);
        if (n > 0 && x <= 16.0) {
            const int k = 18;
            const double s = x * 0.5;
            const double t = s * s;
            double v = 1.0;
            for (int i = k; i > 0; i--) {
                v = 1.0 + v * t / (i * (i + n + 0.5));
            }
            v *= 0.5 / power_factorial_half_shifted_(1.0 / s, n);
            return (e == 0.0) ? v : v * exp(e);
        }
        const double t = 1.0 / x;
        const double p = exp(e + x);
        const double q = exp(e - x);
        double y1 = t * 0.5 * (p - q);
        if (n > 0) {
            double y0 = t * 0.5 * (p + q);
            for (int i = 1; i <= n; i++) {
                const double y = y0 - (2 * i - 1) * t * y1; y0 = y1; y1 = y;
            }
        }
        return y1;
    }
    else {
        const double t = 1.0 / x;
        const double p = exp(e + x);
        const double q = exp(e - x);
        double y0 = t * 0.5 * (p + q);
        if (n < -1) {
            double y1 = t * 0.5 * (p - q);
            for (int i = -2; i >= n; i--) {
                const double y = y1 + (2 * i + 3) * t * y0; y1 = y0; y0 = y;
            }
        }
        return y0;
    }
}
