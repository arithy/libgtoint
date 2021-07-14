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

/*
 * References:
 * - Michael E. Mura, and Peter J. Knowles;
 *   "Improved radial grids for quadrature in molecular density-functional calculations";
 *   Journal of Chemical Physics, Volume 104, Issue 24, Pages 9848-9858; June 1996.
 * - Chris-Kriton Skylaris, Laura Gagliardi, Nicholas C. Handy,
 *   Andrew G. Ioannou, Steven Spencer, Andrew Willetts, and Adrian M. Simper;
 *   "An efficient method for calculating effective core potential integrals which involve projection operators";
 *   Chemical Physics Letters, Volume 296, Issues 5-6, Pages 445-451; November 1998.
 */

#include "integral-ecp-2.h"

#include "function.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef ARRAY_INIT_ALLOC
#define ARRAY_INIT_ALLOC 16
#endif

#define A 7.0
#define L_MIN 5
#define L_MAX 12

#define D_RPI    1.7724538509055160272981674833411  /* sqrt(PI) */
#define D_RPID2  1.2533141373155002512078826424055  /* sqrt(PI/2) */
#define D_1D2RPI 0.28209479177387814347403972578039 /* 1/(2 sqrt(PI)) */

static double g_x_[(1 << 12) - 1];

static double g_r_[(1 << 12) - 1];

static double g_w05_[(1 <<  5) - 1];
static double g_w06_[(1 <<  6) - 1];
static double g_w07_[(1 <<  7) - 1];
static double g_w08_[(1 <<  8) - 1];
static double g_w09_[(1 <<  9) - 1];
static double g_w10_[(1 << 10) - 1];
static double g_w11_[(1 << 11) - 1];
static double g_w12_[(1 << 12) - 1];
static double *const g_w_[] = { g_w05_, g_w06_, g_w07_, g_w08_, g_w09_, g_w10_, g_w11_, g_w12_ };

static void compute_quadrature_table_(void) {
    for (int i = 0, l = 0; l < L_MAX; l++) {
        const int m = 1 << l;
        for (int k = 0; k < m; k++) {
            g_x_[i++] = (double)(2 * k + 1) / (2 * m);
        }
    }
    for (int i = 0; i < (1 << L_MAX) - 1; i++) {
        const double x = g_x_[i];
        g_r_[i] = -A * log(1.0 - x * x * x);
    }
    for (int l = L_MIN; l <= L_MAX; l++) {
        const int j = l - L_MIN;
        const int n = (1 << l) - 1;
        const double c = (3.0 * A * A * A) / (n + 1);
        for (int i = 0; i < n; i++) {
            const double x = g_x_[i];
            const double y = x * x;
            const double z = 1.0 - x * y;
            const double w = log(z);
            g_w_[j][i] = c * y * w * w / z;
        }
    }
}

inline static double factorial2_(int n) {
#define N 20
    static const double t[N] = {
        1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135,
        645120, 2027025, 10321920, 34459425, 185794560, 654729075
    };
    if (n < 0) return 1.0; /* For -1 and safety. */
    if (n < N) return t[n];
    double v = 1.0;
    for (int i = 2 - (n & 1); i <= n; i += 2) v *= i;
    return v;
#undef N
}

inline static double power_(double a, int b) {
    if (b < 0) { b = -b; a = 1.0 / a; }
    double v = 1.0;
    while (b != 0) {
        if (b & 1) v *= a;
        b >>= 1; a *= a;
    }
    return v;
}

inline static double power_scaled_(double a, int b, double s) {
    if (b < 0) { b = -b; a = 1.0 / a; }
    double v = s;
    while (b != 0) {
        if (b & 1) v *= a;
        b >>= 1; a *= a;
    }
    return v;
}

static void ecp_type2_radial_integral_entry__compute_(
    ecp_type2_radial_integral_entry_t *obj, const ecp_type2_radial_integral_index_t *index,
    double gr0, double gr1, double g0, double g1, double gc, double tol, double *wrk
) {
    obj->i = *index;
    if (gr0 > 0.0 && gr1 > 0.0) {
        double v0 = 0.0;
        for (int n0 = 0, l = L_MIN; l <= L_MAX; l++) {
            const int j = l - L_MIN;
            const int n = (1 << l) - 1;
            for (int i = n0; i < n; i++) {
                const double r = g_r_[i];
                const double r2 = r * r;
                wrk[i] = power_scaled_(r, index->k, exp(-gc * r2)) *
                    gtoint__compute_modified_spherical_bessel_function_first_kind(index->l[0], gr0 * r, -g0 * r2) *
                    gtoint__compute_modified_spherical_bessel_function_first_kind(index->l[1], gr1 * r, -g1 * r2);
            }
            double v = 0.0;
            for (int i = 0; i < n; i++) v += g_w_[j][i] * wrk[i];
            if (n0 > 0) {
                if (isnan(v) || isinf(v) || fabs(v - v0) <= fabs(v0) * tol) { v0 = v; break; }
            }
            n0 = n;
            v0 = v;
        }
        obj->v = v0;
    }
    else if (gr0 > 0.0) {
        double v0 = 0.0;
        for (int n0 = 0, l = L_MIN; l <= L_MAX; l++) {
            const int j = l - L_MIN;
            const int n = (1 << l) - 1;
            for (int i = n0; i < n; i++) {
                const double r = g_r_[i];
                const double r2 = r * r;
                wrk[i] = power_scaled_(r, index->k, exp(-(g1 + gc) * r2)) *
                    gtoint__compute_modified_spherical_bessel_function_first_kind(index->l[0], gr0 * r, -g0 * r2);
            }
            double v = 0.0;
            for (int i = 0; i < n; i++) v += g_w_[j][i] * wrk[i];
            if (n0 > 0) {
                if (isnan(v) || isinf(v) || fabs(v - v0) <= fabs(v0) * tol) { v0 = v; break; }
            }
            n0 = n;
            v0 = v;
        }
        obj->v = v0;
    }
    else if (gr1 > 0.0) {
        double v0 = 0.0;
        for (int n0 = 0, l = L_MIN; l <= L_MAX; l++) {
            const int j = l - L_MIN;
            const int n = (1 << l) - 1;
            for (int i = n0; i < n; i++) {
                const double r = g_r_[i];
                const double r2 = r * r;
                wrk[i] = power_scaled_(r, index->k, exp(-(g0 + gc) * r2)) *
                    gtoint__compute_modified_spherical_bessel_function_first_kind(index->l[1], gr1 * r, -g1 * r2);
            }
            double v = 0.0;
            for (int i = 0; i < n; i++) v += g_w_[j][i] * wrk[i];
            if (n0 > 0) {
                if (isnan(v) || isinf(v) || fabs(v - v0) <= fabs(v0) * tol) { v0 = v; break; }
            }
            n0 = n;
            v0 = v;
        }
        obj->v = v0;
    }
    else {
        const int k = index->k;
        obj->v = factorial2_(k + 1) * sqrt(power_(0.5 / (g0 + g1 + gc), k + 3)) * ((k & 1) ? 1.0 : D_RPID2);
    }
}

/*== ecp_type2_radial_integral_array_t ==*/

inline static void ecp_type2_radial_integral_array__clear_(ecp_type2_radial_integral_array_t *obj) {
    obj->p = NULL;
}

inline static bool ecp_type2_radial_integral_array__realloc_(ecp_type2_radial_integral_array_t *obj, size_t num) {
    ecp_type2_radial_integral_entry_t *const p = (ecp_type2_radial_integral_entry_t *)realloc(obj->p, sizeof(ecp_type2_radial_integral_entry_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void ecp_type2_radial_integral_array__free_(ecp_type2_radial_integral_array_t *obj) {
    free(obj->p);
}

void gtoint__ecp_type2_radial_integral_array__initialize(ecp_type2_radial_integral_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    ecp_type2_radial_integral_array__clear_(obj);
}

void gtoint__ecp_type2_radial_integral_array__finalize(ecp_type2_radial_integral_array_t *obj) {
    ecp_type2_radial_integral_array__free_(obj);
}

bool gtoint__ecp_type2_radial_integral_array__resize(ecp_type2_radial_integral_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m <= 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!ecp_type2_radial_integral_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__ecp_type2_radial_integral_array__copy(ecp_type2_radial_integral_array_t *obj, const ecp_type2_radial_integral_array_t *src) {
    if (!gtoint__ecp_type2_radial_integral_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(ecp_type2_radial_integral_entry_t) * src->n);
    return true;
}

void gtoint__ecp_type2_radial_integral_array__move(ecp_type2_radial_integral_array_t *obj, ecp_type2_radial_integral_array_t *src) {
    gtoint__ecp_type2_radial_integral_array__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_radial_integral_array__initialize(src);
}

void gtoint__ecp_type2_radial_integral_array__compact(ecp_type2_radial_integral_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n) m <<= 1;
        if (m >= obj->m) return;
        ecp_type2_radial_integral_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        ecp_type2_radial_integral_array__free_(obj);
        obj->m = 0;
        ecp_type2_radial_integral_array__clear_(obj);
    }
}

/*== ecp_type2_radial_integral_database_t ==*/

inline static int compare_indices_(const ecp_type2_radial_integral_index_t *a, const ecp_type2_radial_integral_index_t *b) {
    if (a->k > b->k) return 1;
    if (a->k < b->k) return -1;
    if (a->l[0] > b->l[0]) return 1;
    if (a->l[0] < b->l[0]) return -1;
    if (a->l[1] > b->l[1]) return 1;
    if (a->l[1] < b->l[1]) return -1;
    return 0;
}

inline static bool ecp_type2_radial_integral_database__find_entry_(const ecp_type2_radial_integral_database_t *obj, const ecp_type2_radial_integral_index_t *index, size_t *out) {
    if (obj->n <= 0) { *out = 0; return false; }
    size_t min = 0;
    size_t max = obj->n - 1;
    for (;;) {
        const size_t mid = min + ((max - min) >> 1);
        const int c = compare_indices_(index, &(obj->a.p[obj->o.p[mid]].i));
        if (c > 0) {
            if (min == mid) { *out = min + 1; return false; }
            min = mid + 1;
        }
        else if (c < 0) {
            if (min == mid) { *out = min; return false; }
            max = mid - 1;
        }
        else {
            *out = mid;
            return true;
        }
    }
}

void gtoint__ecp_type2_radial_integral_database__initialize(ecp_type2_radial_integral_database_t *obj) {
    obj->gr0 = 0.0;
    obj->gr1 = 0.0;
    obj->g0 = 0.0;
    obj->g1 = 0.0;
    obj->gc = 0.0;
    obj->tol = 1e-10;
    obj->n = 0;
    gtoint__size_t_array__initialize(&(obj->o));
    gtoint__ecp_type2_radial_integral_array__initialize(&(obj->a));
    gtoint__double_array__initialize(&(obj->w));
}

void gtoint__ecp_type2_radial_integral_database__finalize(ecp_type2_radial_integral_database_t *obj) {
    gtoint__size_t_array__finalize(&(obj->o));
    gtoint__ecp_type2_radial_integral_array__finalize(&(obj->a));
    gtoint__double_array__finalize(&(obj->w));
}

void gtoint__ecp_type2_radial_integral_database__clear(ecp_type2_radial_integral_database_t *obj) {
    obj->n = 0;
    gtoint__size_t_array__resize(&(obj->o), 0);
    gtoint__ecp_type2_radial_integral_array__resize(&(obj->a), 0);
    gtoint__double_array__resize(&(obj->w), 0);
}

void gtoint__ecp_type2_radial_integral_database__reset(ecp_type2_radial_integral_database_t *obj, double3_t r0c, double3_t r1c, double g0, double g1, double gc) {
    const double gr0 = 2.0 * g0 * sqrt(r0c.x * r0c.x + r0c.y * r0c.y + r0c.z * r0c.z);
    const double gr1 = 2.0 * g1 * sqrt(r1c.x * r1c.x + r1c.y * r1c.y + r1c.z * r1c.z);
    if (obj->gr0 == gr0 && obj->gr1 == gr1 && obj->g0 == g0 && obj->g1 == g1 && obj->gc == gc) return;
    obj->gr0 = gr0;
    obj->gr1 = gr1;
    obj->g0 = g0;
    obj->g1 = g1;
    obj->gc = gc;
    obj->n = 0;
    gtoint__size_t_array__resize(&(obj->o), 0);
    gtoint__ecp_type2_radial_integral_array__resize(&(obj->a), 0);
    gtoint__double_array__resize(&(obj->w), 0);
}

bool gtoint__ecp_type2_radial_integral_database__fetch(ecp_type2_radial_integral_database_t *obj, const ecp_type2_radial_integral_index_t *index, double *out) {
    size_t m, k;
    if (ecp_type2_radial_integral_database__find_entry_(obj, index, &m)) {
        k = obj->o.p[m];
    }
    else {
        k = obj->n;
        if (!gtoint__size_t_array__resize(&(obj->o), obj->n + 1)) return false;
        if (!gtoint__ecp_type2_radial_integral_array__resize(&(obj->a), obj->n + 1)) return false;
        if (!gtoint__double_array__resize(&(obj->w), (1 << L_MAX) - 1)) return false;
        ecp_type2_radial_integral_entry__compute_(
            &(obj->a.p[k]), index, obj->gr0, obj->gr1, obj->g0, obj->g1, obj->gc, obj->tol, obj->w.p
        );
        obj->n++;
        memmove(obj->o.p + m + 1, obj->o.p + m, sizeof(size_t) * (obj->o.n - m - 1));
        obj->o.p[m] = k;
        obj->a.p[k].i = *index;
    }
    if (out) *out = obj->a.p[k].v;
    return true;
}

bool gtoint__ecp_type2_radial_integral_database__copy(ecp_type2_radial_integral_database_t *obj, const ecp_type2_radial_integral_database_t *src) {
    if (!gtoint__size_t_array__copy(&(obj->o), &(src->o))) return false;
    if (!gtoint__ecp_type2_radial_integral_array__copy(&(obj->a), &(src->a))) return false;
    obj->n = src->n;
    obj->gr0 = src->gr0;
    obj->gr1 = src->gr1;
    obj->g0 = src->g0;
    obj->g1 = src->g1;
    obj->gc = src->gc;
    obj->tol = src->tol;
    return true;
}

void gtoint__ecp_type2_radial_integral_database__move(ecp_type2_radial_integral_database_t *obj, ecp_type2_radial_integral_database_t *src) {
    gtoint__ecp_type2_radial_integral_database__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_radial_integral_database__initialize(src);
}

void gtoint__ecp_type2_radial_integral_database__compact(ecp_type2_radial_integral_database_t *obj) {
    gtoint__size_t_array__compact(&(obj->o));
    gtoint__ecp_type2_radial_integral_array__compact(&(obj->a));
    gtoint__double_array__compact(&(obj->w));
}

#ifdef _MSC_VER
#pragma section(".CRT$XCU",read)
#else
__attribute__((constructor))
#endif
static void construct(void) {
    compute_quadrature_table_();
}

#ifdef _MSC_VER
__declspec(allocate(".CRT$XCU")) void (*construct_integral_ecp_2_radial_)(void) = construct;
#ifdef _WIN64
__pragma(comment(linker,"/include:construct_integral_ecp_2_radial_"))
#else
__pragma(comment(linker,"/include:_construct_integral_ecp_2_radial_"))
#endif
#endif
