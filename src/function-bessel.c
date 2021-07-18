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
#include <assert.h>

#ifndef VL
#define VL 256
#endif

#define K 18

inline static double power_factorial_half_shifted_(double a, int n) {
    double v = 0.5;
    for (int i = 1; i <= n; i++) v *= a * (i + 0.5);
    return v;
}

double gtoint__compute_modified_spherical_bessel_function_first_kind(int n, double x, double e) { /* The relative error is about 1e-6. */
    if (n >= 0) {
        if (x == 0.0) return (n > 0) ? 0.0 : (e == 0.0) ? 1.0 : exp(e);
        if (n > 0 && x <= 16.0) {
            const double s = x * 0.5;
            const double t = s * s;
            double v = 1.0;
            for (int i = K; i > 0; i--) {
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

inline static void compute_function_batch_helper_0_(int n, size_t m, const double *restrict x, const double *restrict e, double *restrict v) {
    assert(m <= VL);
    double t[VL];
    double u[VL];
    double w[VL];
    for (size_t i = 0; i < m; i++) {
        const double s = x[i] * 0.5;
        t[i] = s * s;
        u[i] = 1.0 / s;
        v[i] = 1.0;
        w[i] = 1.0;
    }
    for (int k = K; k > 0; k--) {
        for (size_t i = 0; i < m; i++) {
            v[i] = 1.0 + v[i] * t[i] / (k * (k + n + 0.5));
        }
    }
    for (int k = 1; k <= n; k++) {
        for (size_t i = 0; i < m; i++) {
            w[i] *= u[i] * (k + 0.5);
        }
    }
    for (size_t i = 0; i < m; i++) {
        v[i] *= exp(e[i]) / w[i];
    }
}

inline static void compute_function_batch_helper_1_(int n, size_t m, const double *restrict x, const double *restrict e, double *restrict v) {
    assert(m <= VL);
    double t[VL];
    double u[VL];
    for (size_t i = 0; i < m; i++) {
        t[i] = 1.0 / x[i];
        const double p = exp(e[i] + x[i]);
        const double q = exp(e[i] - x[i]);
        v[i] = t[i] * 0.5 * (p - q);
        u[i] = t[i] * 0.5 * (p + q);
    }
    for (int k = 1; k <= n; k++) {
        for (size_t i = 0; i < m; i++) {
            const double w = u[i] - (2 * k - 1) * t[i] * v[i]; u[i] = v[i]; v[i] = w;
        }
    }
}

void gtoint__compute_modified_spherical_bessel_function_first_kind_batch(int n, size_t m, const double *restrict x, const double *restrict e, double *restrict v) {
    if (n > 0) {
        for (size_t i = 0; i < m; i++) {
            if (x[i] == 0.0) v[i] = 0.0;
        }
        {
            size_t w_i[VL];
            double w_x[VL];
            double w_e[VL];
            double w_v[VL];
            size_t l = 0;
            for (size_t i = 0; i < m; i++) {
                if (x[i] <= 16.0) {
                    w_i[l] = i;
                    w_x[l] = x[i];
                    w_e[l] = e[i];
                    l++;
                    if (l == VL) {
                        compute_function_batch_helper_0_(n, l, w_x, w_e, w_v);
                        for (size_t j = 0; j < l; j++) v[w_i[j]] = w_v[j];
                        l = 0;
                    }
                }
            }
            if (l > 0) {
                compute_function_batch_helper_0_(n, l, w_x, w_e, w_v);
                for (size_t j = 0; j < l; j++) v[w_i[j]] = w_v[j];
            }
        }
        {
            size_t w_i[VL];
            double w_x[VL];
            double w_e[VL];
            double w_v[VL];
            size_t l = 0;
            for (size_t i = 0; i < m; i++) {
                if (x[i] > 16.0) {
                    w_i[l] = i;
                    w_x[l] = x[i];
                    w_e[l] = e[i];
                    l++;
                    if (l == VL) {
                        compute_function_batch_helper_1_(n, l, w_x, w_e, w_v);
                        for (size_t j = 0; j < l; j++) v[w_i[j]] = w_v[j];
                        l = 0;
                    }
                }
            }
            if (l > 0) {
                compute_function_batch_helper_1_(n, l, w_x, w_e, w_v);
                for (size_t j = 0; j < l; j++) v[w_i[j]] = w_v[j];
            }
        }
    }
    else if (n == 0) {
        for (size_t i = 0; i < m; i++) {
            v[i] = (x[i] > 0.0) ? 0.5 * (exp(e[i] + x[i]) - exp(e[i] - x[i])) / x[i] : exp(e[i]);
        }
    }
    else if (n == -1) {
        for (size_t i = 0; i < m; i++) {
            v[i] = 0.5 * (exp(e[i] + x[i]) + exp(e[i] - x[i])) / x[i];
        }
    }
    else {
        double t[VL];
        double y0[VL];
        double y1[VL];
        for (size_t i = 0; i < m; i += VL) {
            const size_t l = (m - i < VL) ? m - i : VL;
            for (size_t j = 0; j < l; j++) {
                const size_t k = i + j;
                t[j] = 1.0 / x[k];
                const double p = exp(e[k] + x[k]);
                const double q = exp(e[k] - x[k]);
                y0[j] = t[j] * 0.5 * (p + q);
                y1[j] = t[j] * 0.5 * (p - q);
            }
            for (int h = -2; h >= n; h--) {
                for (size_t j = 0; j < l; j++) {
                    const double y = y1[j] + (2 * h + 3) * t[j] * y0[j]; y1[j] = y0[j]; y0[j] = y;
                }
            }
            for (size_t j = 0; j < l; j++) {
                v[i + j] = y0[j];
            }
        }
    }
}
