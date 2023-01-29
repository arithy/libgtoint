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

/*
 * References:
 * - Larry E. McMurchie, and Ernest R. Davidson;
 *   "Calculation of integrals over ab initio pseudopotentials";
 *   Journal of Computational Physics, Volume 44, Issue 2, Pages 289-301; December 1981.
 * - Simon C. McKenzie, Evgeny Epifanovsky, Giuseppe M. J. Barca, Andrew T. B. Gilbert, and Peter M. W. Gill;
 *   "Efficient Method for Calculating Effective Core Potential Integrals";
 *   The Journal of Physical Chemistry A, Volume 122, Issue 11, Pages 3066-3075; February 2018.
 */

#include "gtoint-private.h"

#include <math.h>
#include <assert.h>

#define D_4PI   12.566370614359172953850573533118 /* 4 PI */
#define D_8RPI3 44.546623974653662762278543856951 /* 8 sqrt(PI)^3 */
#define D_16PI2 157.91367041742973790135185599802 /* 16 PI^2 */
#define D_RPI3  5.5683279968317078452848179821188 /* PI^(3/2) */
#define D_2E    5.4365636569180904707205749427053 /* 2 E */
#define D_SINH1 1.1752011936438014568823818505956 /* sinh(1) */
#define D_PI3D  1.1483806177888822287213450024852 /* (PI/3)^3 */

inline static double power_(double a, int b) {
    if (b < 0) { b = -b; a = 1.0 / a; }
    double v = 1.0;
    while (b != 0) {
        if (b & 1) v *= a;
        b >>= 1; a *= a;
    }
    return v;
}

inline static double binomial_(int n, int k) {
#define N 10
    static const double t[(N * (N + 1)) >> 1] = {
        1,
        1, 1,
        1, 2, 1,
        1, 3, 3, 1,
        1, 4, 6, 4, 1,
        1, 5, 10, 10, 5, 1,
        1, 6, 15, 20, 15, 6, 1,
        1, 7, 21, 35, 35, 21, 7, 1,
        1, 8, 28, 56, 70, 56, 28, 8, 1,
        1, 9, 36, 84, 126, 126, 84, 36, 9, 1
    };
    if (n < 0) return 1.0; /* For safety. */
    if (n < N) return t[((n * (n + 1)) >> 1) + k];
    double v = 1.0;
    for (int i = k + 1; i <= n; i++) v *= i;
    double w = 1.0;
    for (int i = 1; i <= n - k; i++) w *= i;
    return v / w;
#undef N
}

gtoint_error_t gtoint__compute_scalar_ecp_type2_integrals(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    const double3_t *pc, int ac, int rc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
) {
#define NVAR 1
#define EXPAND_VRR_0(xyz) /* reference variables: itg, s, is, ng0, ng1, ngc, ng01c, g0, r0c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double r = 2.0 * r0c.xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = r * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2.k[0].xyz++; \
            t.i.ecp2.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = 2.0 * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2.d[0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.d[0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = -2.0 * s.i.ecp2.d[0].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2.a[0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.a[0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = -s.i.ecp2.a[0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = 2.0 * s.i.ecp2.d[2].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
    }
#define EXPAND_VRR_1(xyz) /* reference variables: itg, s, is, ng0, ng1, ngc, ng01c, g1, r1c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double r = 2.0 * r1c.xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = r * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2.k[1].xyz++; \
            t.i.ecp2.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = 2.0 * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2.d[1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.d[1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = -2.0 * s.i.ecp2.d[1].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2.a[1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.a[1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = -s.i.ecp2.a[1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = 2.0 * s.i.ecp2.d[2].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
    }
#define EXPAND_VRR_C(xyz) /* reference variables: itg, s, is, ng0, ng1, ngc, ng01c, g0, g1, r0c, r1c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double r0 = -2.0 * r0c.xyz; \
            const double r1 = -2.0 * r1c.xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = r0 * g0[i0] + r1 * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2.k[0].xyz++; \
            t.i.ecp2.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = -2.0 * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2.k[1].xyz++; \
            t.i.ecp2.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = -2.0 * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2.a[0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.a[0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = s.i.ecp2.a[0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2.a[1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.a[1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = s.i.ecp2.a[1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = -2.0 * s.i.ecp2.d[2].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * (g0[i0] + g1[i1]); \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
    }
    if (!gtoint__double_array__resize(&(itg->v), na0 * na1 * nd)) return GTOINT_ERROR_MEMORY;
    const double3_t r0c = { pc->x - p0->x, pc->y - p0->y, pc->z - p0->z };
    const double3_t r1c = { pc->x - p1->x, pc->y - p1->y, pc->z - p1->z };
    const double r0cw = r0c.x * r0c.x + r0c.y * r0c.y + r0c.z * r0c.z;
    const double r1cw = r1c.x * r1c.x + r1c.y * r1c.y + r1c.z * r1c.z;
    if (!gtoint__ecp_type2_spherical_factor_database_array__resize(&(itg->ecp.s), 2)) return GTOINT_ERROR_MEMORY;
    if (r0cw > 0.0) {
        const double r = 1.0 / sqrt(r0cw);
        gtoint__ecp_type2_spherical_factor_database__reset(&(itg->ecp.s.p[0]), double3__new(r0c.x * r, r0c.y * r, r0c.z * r));
    }
    if (r1cw > 0.0) {
        const double r = 1.0 / sqrt(r1cw);
        gtoint__ecp_type2_spherical_factor_database__reset(&(itg->ecp.s.p[1]), double3__new(r1c.x * r, r1c.y * r, r1c.z * r));
    }
    const size_t ng01 = ng0 * ng1;
    const size_t ng01c = ng01 * ngc;
    if (!gtoint__cache__reset(&(itg->c), ng01c)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_array__resize(&(itg->w), ng01 * NVAR)) return GTOINT_ERROR_MEMORY;
    size_t ivar = 0;
    double *const e01 = itg->w.p + ng01 * ivar++;
    assert(ivar == NVAR);
    for (size_t i0 = 0; i0 < ng0; i0++) {
        const double e0 = exp(-g0[i0] * r0cw);
    for (size_t i1 = 0; i1 < ng1; i1++) {
        e01[i0 + ng0 * i1] = e0;
    }
    }
    for (size_t i1 = 0; i1 < ng1; i1++) {
        const double e1 = exp(-g1[i1] * r1cw);
    for (size_t i0 = 0; i0 < ng0; i0++) {
        e01[i0 + ng0 * i1] *= e1;
    }
    }
    size_t io = 0;
    for (size_t id = 0; id < nd; id++) {
    for (size_t ia1 = 0; ia1 < na1; ia1++) {
    for (size_t ia0 = 0; ia0 < na0; ia0++) {
        itg->v.p[io] = 0.0;
        gtoint__stack__reset(&(itg->s), ng01c);
        {
            stack_index_t t = { 0 };
            t.i.ecp2.a[0] = a0[ia0];
            t.i.ecp2.a[1] = a1[ia1];
            t.i.ecp2.d[0] = d0[id];
            t.i.ecp2.d[1] = d1[id];
            t.i.ecp2.d[2] = dc[id];
            t.i.ecp2.r = rc - 2;
            t.o = STACK_VOID_INDEX;
            t.b = false;
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY;
            size_t it;
            gtoint__stack__write(&(itg->s), &t, &it);
            double *const c = gtoint__stack__coefficients(&(itg->s), it);
            double *const v = gtoint__stack__integrals(&(itg->s), it);
            size_t i = 0;
            for (size_t i1 = 0; i1 < ng1; i1++) {
            for (size_t i0 = 0; i0 < ng0; i0++) {
            for (size_t ic = 0; ic < ngc; ic++, i++) {
                c[i] = c0[i0 + ng0 * ia0] * c1[i1 + ng1 * ia1] * cc[ic];
                v[i] = 0.0;
            }
            }
            }
        }
        while (!gtoint__stack__is_empty(&(itg->s))) {
            stack_index_t s;
            size_t is;
            gtoint__stack__read(&(itg->s), &s, &is);
            if (s.b) {
                const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                const double *const v = gtoint__stack__integrals(&(itg->s), is);
                double *w;
                if (!gtoint__cache__reference_to_store_ecp2(&(itg->c), &s.i, &w)) return GTOINT_ERROR_MEMORY;
                for (size_t i = 0; i < ng01c; i++) w[i] = v[i];
                if (s.o == STACK_VOID_INDEX) {
                    double o = 0.0;
                    for (size_t i = 0; i < ng01c; i++) o += c[i] * v[i];
                    itg->v.p[io] += o;
                }
                else {
                    double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                    for (size_t i = 0; i < ng01c; i++) o[i] += c[i] * v[i];
                }
                gtoint__stack__pop(&(itg->s));
            }
            else {
                const double *w;
                if (gtoint__cache__reference_to_fetch_ecp2(&(itg->c), &s.i, &w)) {
                    const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                    if (s.o == STACK_VOID_INDEX) {
                        double o = 0.0;
                        for (size_t i = 0; i < ng01c; i++) o += c[i] * w[i];
                        itg->v.p[io] += o;
                    }
                    else {
                        double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                        for (size_t i = 0; i < ng01c; i++) o[i] += c[i] * w[i];
                    }
                    gtoint__stack__pop(&(itg->s));
                }
                else {
                    {
                        s.b = true;
                        gtoint__stack__write(&(itg->s), &s, &is);
                        double *const v = gtoint__stack__integrals(&(itg->s), is);
                        for (size_t i = 0; i < ng01c; i++) v[i] = 0.0;
                    }
                    if (s.i.ecp2.d[0].x > 0) {
                        s.i.ecp2.d[0].x--;
                        EXPAND_VRR_0(x);
                    }
                    else if (s.i.ecp2.d[0].y > 0) {
                        s.i.ecp2.d[0].y--;
                        EXPAND_VRR_0(y);
                    }
                    else if (s.i.ecp2.d[0].z > 0) {
                        s.i.ecp2.d[0].z--;
                        EXPAND_VRR_0(z);
                    }
                    else if (s.i.ecp2.d[1].x > 0) {
                        s.i.ecp2.d[1].x--;
                        EXPAND_VRR_1(x);
                    }
                    else if (s.i.ecp2.d[1].y > 0) {
                        s.i.ecp2.d[1].y--;
                        EXPAND_VRR_1(y);
                    }
                    else if (s.i.ecp2.d[1].z > 0) {
                        s.i.ecp2.d[1].z--;
                        EXPAND_VRR_1(z);
                    }
                    else if (s.i.ecp2.d[2].x > 0) {
                        s.i.ecp2.d[2].x--;
                        EXPAND_VRR_C(x);
                    }
                    else if (s.i.ecp2.d[2].y > 0) {
                        s.i.ecp2.d[2].y--;
                        EXPAND_VRR_C(y);
                    }
                    else if (s.i.ecp2.d[2].z > 0) {
                        s.i.ecp2.d[2].z--;
                        EXPAND_VRR_C(z);
                    }
                    else {
                        double *v;
                        if (!gtoint__cache__reference_to_store_ecp2(&(itg->c), &s.i, &v)) return GTOINT_ERROR_MEMORY;
                        if (r0cw > 0.0 && r1cw > 0.0) {
                            const int h0 = s.i.ecp2.k[0].x + s.i.ecp2.k[0].y + s.i.ecp2.k[0].z;
                            const int h1 = s.i.ecp2.k[1].x + s.i.ecp2.k[1].y + s.i.ecp2.k[1].z;
                            size_t i = 0;
                            for (size_t i1 = 0; i1 < ng1; i1++) {
                            for (size_t i0 = 0; i0 < ng0; i0++) {
                            for (size_t ic = 0; ic < ngc; ic++, i++) {
                                { /* screening by approximated estimation */
                                    /* Reference:
                                     * "Efficient Method for Calculating Effective Core Potential Integrals"
                                     *  Simon C. McKenzie, Evgeny Epifanovsky, Giuseppe M. J. Barca, Andrew T. B. Gilbert, and Peter M. W. Gill
                                     *  J. Phys. Chem. A 2018, 122, 3066-3075
                                     */
                                    const int a0s = s.i.ecp2.a[0].x + s.i.ecp2.a[0].y + s.i.ecp2.a[0].z;
                                    const int a1s = s.i.ecp2.a[1].x + s.i.ecp2.a[1].y + s.i.ecp2.a[1].z;
                                    const double g0c = g0[i0] + gc[ic];
                                    const double g1c = g1[i1] + gc[ic];
                                    const double cg0 = a0s * g0c * g0c / (2.0 * g0[i0] * (gc[ic] * gc[ic] * r0cw + a0s * g0c));
                                    const double cg1 = a1s * g1c * g1c / (2.0 * g1[i1] * (gc[ic] * gc[ic] * r1cw + a1s * g1c));
                                    const double ge0 = (1.0 - cg0) * g0[i0];
                                    const double ge1 = (1.0 - cg1) * g1[i1];
                                    const double ges = 1.0 / (ge0 + ge1 + gc[ic]);
                                    const double t = 2.0 * ge0 * ge1 * sqrt(r0cw * r1cw) * ges;
                                    if (
                                        (2 * ac + 1) * (2 * ac + 1) *
                                        sqrt(power_(a0s / (D_2E * g0[i0] * cg0), a0s)) *
                                        sqrt(power_(a1s / (D_2E * g1[i1] * cg1), a1s)) *
                                        (D_RPI3 / (gc[ic] * sqrt(gc[ic]))) *
                                        exp(ge0 * r0cw * (ge0 * ges - 1.0) + ge1 * r1cw * (ge1 * ges - 1.0)) *
                                        ((t <= 1.0) ? D_SINH1 : exp(t) / (2.0 * t)) <= itg->cut
                                    ) { v[i] = 0.0; continue; }
                                }
                                gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), r0c, r1c, g0[i0], g1[i1], gc[ic]);
                                double v0x = 0.0;
                                for (int k0x = 0; k0x <= s.i.ecp2.a[0].x; k0x++) {
                                double v0y = 0.0;
                                for (int k0y = 0; k0y <= s.i.ecp2.a[0].y; k0y++) {
                                double v0z = 0.0;
                                for (int k0z = 0; k0z <= s.i.ecp2.a[0].z; k0z++) {
                                double v1x = 0.0;
                                for (int k1x = 0; k1x <= s.i.ecp2.a[1].x; k1x++) {
                                double v1y = 0.0;
                                for (int k1y = 0; k1y <= s.i.ecp2.a[1].y; k1y++) {
                                double v1z = 0.0;
                                for (int k1z = 0; k1z <= s.i.ecp2.a[1].z; k1z++) {
                                const int k0 = k0x + k0y + k0z;
                                const int k1 = k1x + k1y + k1z;
                                double v01 = 0.0;
                                for (int l0 = 0; l0 <= ac + k0 + h0; l0++) {
                                    size_t sh0, si0;
                                    if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[0]), &(itg->ecp.h), l0, &sh0)) return GTOINT_ERROR_MEMORY;
                                    {
                                        const ecp_type2_angular_integral_index_t t = {
                                            { l0, ac }, { k0x + s.i.ecp2.k[0].x, k0y + s.i.ecp2.k[0].y, k0z + s.i.ecp2.k[0].z }
                                        };
                                        if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                                    }
                                for (int l1 = 0; l1 <= ac + k1 + h1; l1++) {
                                    size_t sh1, si1;
                                    if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[1]), &(itg->ecp.h), l1, &sh1)) return GTOINT_ERROR_MEMORY;
                                    {
                                        const ecp_type2_angular_integral_index_t t = {
                                            { l1, ac }, { k1x + s.i.ecp2.k[1].x, k1y + s.i.ecp2.k[1].y, k1z + s.i.ecp2.k[1].z }
                                        };
                                        if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                                    }
                                    const double *const th0 = itg->ecp.s.p[0].a.p[sh0].v.p;
                                    const double *const th1 = itg->ecp.s.p[1].a.p[sh1].v.p;
                                    const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                    const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                    double vm = 0.0;
                                    for (int m = 0; m <= 2 * ac; m++) {
                                        double vj0 = 0.0;
                                        for (int j0 = 0; j0 <= 2 * l0; j0++) {
                                            vj0 += th0[j0] * ti0[j0 + (2 * l0 + 1) * m];
                                        }
                                        double vj1 = 0.0;
                                        for (int j1 = 0; j1 <= 2 * l1; j1++) {
                                            vj1 += th1[j1] * ti1[j1 + (2 * l1 + 1) * m];
                                        }
                                        vm += vj0 * vj1;
                                    }
                                    if (vm != 0.0) {
                                        const ecp_type2_radial_integral_index_t t = { { l0, l1 }, s.i.ecp2.r + k0 + k1 };
                                        double q;
                                        if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                        v01 += vm * (((l0 + l1) & 1) ? -q : q);
                                    }
                                }
                                }
                                v1z += v01 * binomial_(s.i.ecp2.a[1].z, k1z) * power_(r1c.z, s.i.ecp2.a[1].z - k1z);
                                }
                                v1y += v1z * binomial_(s.i.ecp2.a[1].y, k1y) * power_(r1c.y, s.i.ecp2.a[1].y - k1y);
                                }
                                v1x += v1y * binomial_(s.i.ecp2.a[1].x, k1x) * power_(r1c.x, s.i.ecp2.a[1].x - k1x);
                                }
                                v0z += v1x * binomial_(s.i.ecp2.a[0].z, k0z) * power_(r0c.z, s.i.ecp2.a[0].z - k0z);
                                }
                                v0y += v0z * binomial_(s.i.ecp2.a[0].y, k0y) * power_(r0c.y, s.i.ecp2.a[0].y - k0y);
                                }
                                v0x += v0y * binomial_(s.i.ecp2.a[0].x, k0x) * power_(r0c.x, s.i.ecp2.a[0].x - k0x);
                                }
                                v[i] = v0x * e01[i0 + ng0 * i1] * D_16PI2;
                            }
                            }
                            }
                        }
                        else if (r0cw > 0.0) {
                            const int h0 = s.i.ecp2.k[0].x + s.i.ecp2.k[0].y + s.i.ecp2.k[0].z;
                            const int k1 = s.i.ecp2.a[1].x + s.i.ecp2.a[1].y + s.i.ecp2.a[1].z;
                            size_t si1;
                            {
                                const ecp_type2_angular_integral_index_t t = {
                                    { 0, ac }, { s.i.ecp2.a[1].x + s.i.ecp2.k[1].x, s.i.ecp2.a[1].y + s.i.ecp2.k[1].y, s.i.ecp2.a[1].z + s.i.ecp2.k[1].z }
                                };
                                if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                            }
                            size_t i = 0;
                            for (size_t i1 = 0; i1 < ng1; i1++) {
                            for (size_t i0 = 0; i0 < ng0; i0++) {
                            for (size_t ic = 0; ic < ngc; ic++, i++) {
                                { /* screening by approximated estimation */
                                    const int a0s = s.i.ecp2.a[0].x + s.i.ecp2.a[0].y + s.i.ecp2.a[0].z;
                                    const double g0c = g0[i0] + gc[ic];
                                    const double cg0 = a0s * g0c * g0c / (2.0 * g0[i0] * (gc[ic] * gc[ic] * r0cw + a0s * g0c));
                                    const double ge0 = (1.0 - cg0) * g0[i0];
                                    const double t = (a0s + 3) / g0[i0];
                                    if (
                                        (2 * ac + 1) * (2 * ac + 1) *
                                        sqrt(power_(a0s * t / ((D_2E * D_2E) * g0[i0] * cg0), a0s) * t * D_PI3D) * t *
                                        exp(-ge0 * gc[ic] * r0cw / (ge0 + gc[ic])) <= itg->cut
                                    ) { v[i] = 0.0; continue; }
                                }
                                gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), r0c, r1c, g0[i0], g1[i1], gc[ic]);
                                double v0x = 0.0;
                                for (int k0x = 0; k0x <= s.i.ecp2.a[0].x; k0x++) {
                                double v0y = 0.0;
                                for (int k0y = 0; k0y <= s.i.ecp2.a[0].y; k0y++) {
                                double v0z = 0.0;
                                for (int k0z = 0; k0z <= s.i.ecp2.a[0].z; k0z++) {
                                const int k0 = k0x + k0y + k0z;
                                double v01 = 0.0;
                                for (int l0 = 0; l0 <= ac + k0 + h0; l0++) {
                                    size_t sh0, si0;
                                    if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[0]), &(itg->ecp.h), l0, &sh0)) return GTOINT_ERROR_MEMORY;
                                    {
                                        const ecp_type2_angular_integral_index_t t = {
                                            { l0, ac }, { k0x + s.i.ecp2.k[0].x, k0y + s.i.ecp2.k[0].y, k0z + s.i.ecp2.k[0].z }
                                        };
                                        if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                                    }
                                    const double *const th0 = itg->ecp.s.p[0].a.p[sh0].v.p;
                                    const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                    const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                    double vm = 0.0;
                                    for (int m = 0; m <= 2 * ac; m++) {
                                        double vj0 = 0.0;
                                        for (int j0 = 0; j0 <= 2 * l0; j0++) {
                                            vj0 += th0[j0] * ti0[j0 + (2 * l0 + 1) * m];
                                        }
                                        vm += vj0 * ti1[m];
                                    }
                                    if (vm != 0.0) {
                                        const ecp_type2_radial_integral_index_t t = { { l0, 0 }, s.i.ecp2.r + k0 + k1 };
                                        double q;
                                        if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                        v01 += vm * ((l0 & 1) ? -q : q);
                                    }
                                }
                                v0z += v01 * binomial_(s.i.ecp2.a[0].z, k0z) * power_(r0c.z, s.i.ecp2.a[0].z - k0z);
                                }
                                v0y += v0z * binomial_(s.i.ecp2.a[0].y, k0y) * power_(r0c.y, s.i.ecp2.a[0].y - k0y);
                                }
                                v0x += v0y * binomial_(s.i.ecp2.a[0].x, k0x) * power_(r0c.x, s.i.ecp2.a[0].x - k0x);
                                }
                                v[i] = v0x * e01[i0 + ng0 * i1] * D_8RPI3;
                            }
                            }
                            }
                        }
                        else if (r1cw > 0.0) {
                            const int h1 = s.i.ecp2.k[1].x + s.i.ecp2.k[1].y + s.i.ecp2.k[1].z;
                            const int k0 = s.i.ecp2.a[0].x + s.i.ecp2.a[0].y + s.i.ecp2.a[0].z;
                            size_t si0;
                            {
                                const ecp_type2_angular_integral_index_t t = {
                                    { 0, ac }, { s.i.ecp2.a[0].x + s.i.ecp2.k[0].x, s.i.ecp2.a[0].y + s.i.ecp2.k[0].y, s.i.ecp2.a[0].z + s.i.ecp2.k[0].z }
                                };
                                if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                            }
                            size_t i = 0;
                            for (size_t i1 = 0; i1 < ng1; i1++) {
                            for (size_t i0 = 0; i0 < ng0; i0++) {
                            for (size_t ic = 0; ic < ngc; ic++, i++) {
                                { /* screening by approximated estimation */
                                    const int a1s = s.i.ecp2.a[1].x + s.i.ecp2.a[1].y + s.i.ecp2.a[1].z;
                                    const double g1c = g1[i1] + gc[ic];
                                    const double cg1 = a1s * g1c * g1c / (2.0 * g1[i1] * (gc[ic] * gc[ic] * r1cw + a1s * g1c));
                                    const double ge1 = (1.0 - cg1) * g1[i1];
                                    const double t = (a1s + 3) / g1[i1];
                                    if (
                                        (2 * ac + 1) * (2 * ac + 1) *
                                        sqrt(power_(a1s * t / ((D_2E * D_2E) * g1[i1] * cg1), a1s) * t * D_PI3D) * t *
                                        exp(-ge1 * gc[ic] * r1cw / (ge1 + gc[ic])) <= itg->cut
                                    ) { v[i] = 0.0; continue; }
                                }
                                gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), r0c, r1c, g0[i0], g1[i1], gc[ic]);
                                double v1x = 0.0;
                                for (int k1x = 0; k1x <= s.i.ecp2.a[1].x; k1x++) {
                                double v1y = 0.0;
                                for (int k1y = 0; k1y <= s.i.ecp2.a[1].y; k1y++) {
                                double v1z = 0.0;
                                for (int k1z = 0; k1z <= s.i.ecp2.a[1].z; k1z++) {
                                const int k1 = k1x + k1y + k1z;
                                double v01 = 0.0;
                                for (int l1 = 0; l1 <= ac + k1 + h1; l1++) {
                                    size_t sh1, si1;
                                    if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[1]), &(itg->ecp.h), l1, &sh1)) return GTOINT_ERROR_MEMORY;
                                    {
                                        const ecp_type2_angular_integral_index_t t = {
                                            { l1, ac }, { k1x + s.i.ecp2.k[1].x, k1y + s.i.ecp2.k[1].y, k1z + s.i.ecp2.k[1].z }
                                        };
                                        if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                                    }
                                    const double *const th1 = itg->ecp.s.p[1].a.p[sh1].v.p;
                                    const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                    const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                    double vm = 0.0;
                                    for (int m = 0; m <= 2 * ac; m++) {
                                        double vj1 = 0.0;
                                        for (int j1 = 0; j1 <= 2 * l1; j1++) {
                                            vj1 += th1[j1] * ti1[j1 + (2 * l1 + 1) * m];
                                        }
                                        vm += vj1 * ti0[m];
                                    }
                                    if (vm != 0.0) {
                                        const ecp_type2_radial_integral_index_t t = { { 0, l1 }, s.i.ecp2.r + k0 + k1 };
                                        double q;
                                        if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                        v01 += vm * ((l1 & 1) ? -q : q);
                                    }
                                }
                                v1z += v01 * binomial_(s.i.ecp2.a[1].z, k1z) * power_(r1c.z, s.i.ecp2.a[1].z - k1z);
                                }
                                v1y += v1z * binomial_(s.i.ecp2.a[1].y, k1y) * power_(r1c.y, s.i.ecp2.a[1].y - k1y);
                                }
                                v1x += v1y * binomial_(s.i.ecp2.a[1].x, k1x) * power_(r1c.x, s.i.ecp2.a[1].x - k1x);
                                }
                                v[i] = v1x * e01[i0 + ng0 * i1] * D_8RPI3;
                            }
                            }
                            }
                        }
                        else {
                            const int k0 = s.i.ecp2.a[0].x + s.i.ecp2.a[0].y + s.i.ecp2.a[0].z;
                            const int k1 = s.i.ecp2.a[1].x + s.i.ecp2.a[1].y + s.i.ecp2.a[1].z;
                            size_t si0;
                            {
                                const ecp_type2_angular_integral_index_t t = {
                                    { 0, ac }, { s.i.ecp2.a[0].x + s.i.ecp2.k[0].x, s.i.ecp2.a[0].y + s.i.ecp2.k[0].y, s.i.ecp2.a[0].z + s.i.ecp2.k[0].z }
                                };
                                if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                            }
                            size_t si1;
                            {
                                const ecp_type2_angular_integral_index_t t = {
                                    { 0, ac }, { s.i.ecp2.a[1].x + s.i.ecp2.k[1].x, s.i.ecp2.a[1].y + s.i.ecp2.k[1].y, s.i.ecp2.a[1].z + s.i.ecp2.k[1].z }
                                };
                                if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                            }
                            size_t i = 0;
                            for (size_t i1 = 0; i1 < ng1; i1++) {
                            for (size_t i0 = 0; i0 < ng0; i0++) {
                            for (size_t ic = 0; ic < ngc; ic++, i++) {
                                gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), r0c, r1c, g0[i0], g1[i1], gc[ic]);
                                const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                double vm = 0.0;
                                for (int m = 0; m <= 2 * ac; m++) {
                                    vm += ti0[m] * ti1[m];
                                }
                                if (vm != 0.0) {
                                    const ecp_type2_radial_integral_index_t t = { { 0, 0 }, s.i.ecp2.r + k0 + k1 };
                                    double q;
                                    if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                    v[i] = vm * q * D_4PI;
                                }
                                else {
                                    v[i] = 0.0;
                                }
                            }
                            }
                            }
                        }
                        const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                        if (s.o == STACK_VOID_INDEX) {
                            double o = 0.0;
                            for (size_t i = 0; i < ng01c; i++) o += c[i] * v[i];
                            itg->v.p[io] += o;
                        }
                        else {
                            double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                            for (size_t i = 0; i < ng01c; i++) o[i] += c[i] * v[i];
                        }
                        gtoint__stack__pop(&(itg->s));
                    }
                }
            }
        }
        io++;
    }
    }
    }
    return GTOINT_ERROR_OK;
#undef NVAR
#undef EXPAND_VRR_0
#undef EXPAND_VRR_1
#undef EXPAND_VRR_C
}

gtoint_error_t gtoint__compute_weighted_scalar_ecp_type2_integrals(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0, const double3_t *pw0, double gw0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1, const double3_t *pw1, double gw1,
    const double3_t *pc, int ac, int rc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
) {
#define NVAR 1
#define EXPAND_VRR_0(xyz) /* reference variables: itg, s, is, ng0, ng1, ngc, ng01c, g0, r0c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double r = 2.0 * r0c.xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = r * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2_w.k[0].xyz++; \
            t.i.ecp2_w.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = 2.0 * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2_w.d[0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.d[0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = -2.0 * s.i.ecp2_w.d[0].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2_w.a[0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.a[0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = -s.i.ecp2_w.a[0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = 2.0 * s.i.ecp2_w.d[2].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g0[i0]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
    }
#define EXPAND_VRR_1(xyz) /* reference variables: itg, s, is, ng0, ng1, ngc, ng01c, g1, r1c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double r = 2.0 * r1c.xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = r * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2_w.k[1].xyz++; \
            t.i.ecp2_w.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = 2.0 * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2_w.d[1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.d[1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = -2.0 * s.i.ecp2_w.d[1].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2_w.a[1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.a[1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = -s.i.ecp2_w.a[1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = 2.0 * s.i.ecp2_w.d[2].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * g1[i1]; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
    }
#define EXPAND_VRR_W0(xyz) /* reference variables: itg, s, is, ng01c, gw0, rw0c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = 2.0 * gw0 * rw0c.xyz; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2_w.k[0].xyz++; \
            t.i.ecp2_w.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = 2.0 * gw0; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.dw[0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.dw[0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = -2.0 * gw0 * s.i.ecp2_w.dw[0].xyz; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = 2.0 * gw0 * s.i.ecp2_w.d[2].xyz; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
    }
#define EXPAND_VRR_W1(xyz) /* reference variables: itg, s, is, ng01c, gw1, rw1c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = 2.0 * gw1 * rw1c.xyz; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2_w.k[1].xyz++; \
            t.i.ecp2_w.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = 2.0 * gw1; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.dw[1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.dw[1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = -2.0 * gw1 * s.i.ecp2_w.dw[1].xyz; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = 2.0 * gw1 * s.i.ecp2_w.d[2].xyz; \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = u; \
                v[i] = 0.0; \
            } \
        } \
    }
#define EXPAND_VRR_C(xyz) /* reference variables: itg, s, is, ng0, ng1, ngc, ng01c, g0, g1, r0c, r1c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double r0 = -2.0 * r0c.xyz; \
            const double r1 = -2.0 * r1c.xyz; \
            const double u = -2.0 * (gw0 * rw0c.xyz + gw1 * rw1c.xyz); \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = r0 * g0[i0] + r1 * g1[i1] + u; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2_w.k[0].xyz++; \
            t.i.ecp2_w.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = -2.0 * gw0; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = -2.0 * g0[i0] + u; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.ecp2_w.k[1].xyz++; \
            t.i.ecp2_w.r++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double u = -2.0 * gw1; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = -2.0 * g1[i1] + u; \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
        if (s.i.ecp2_w.a[0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.a[0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = s.i.ecp2_w.a[0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.a[1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.a[1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++)  { \
                c[i] = s.i.ecp2_w.a[1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp2_w.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp2_w.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            const double d = -2.0 * s.i.ecp2_w.d[2].xyz; \
            size_t i = 0; \
            for (size_t i1 = 0; i1 < ng1; i1++) { \
            for (size_t i0 = 0; i0 < ng0; i0++) { \
            for (size_t ic = 0; ic < ngc; ic++, i++) { \
                c[i] = d * (g0[i0] + g1[i1] + gw0 + gw1); \
                v[i] = 0.0; \
            } \
            } \
            } \
        } \
    }
    if (!gtoint__double_array__resize(&(itg->v), na0 * na1 * nd)) return GTOINT_ERROR_MEMORY;
    const double3_t r0c = { pc->x - p0->x, pc->y - p0->y, pc->z - p0->z };
    const double3_t r1c = { pc->x - p1->x, pc->y - p1->y, pc->z - p1->z };
    const double3_t rw0c = { pc->x - pw0->x, pc->y - pw0->y, pc->z - pw0->z };
    const double3_t rw1c = { pc->x - pw1->x, pc->y - pw1->y, pc->z - pw1->z };
    const double3_t rw00 = { p0->x - pw0->x, p0->y - pw0->y, p0->z - pw0->z };
    const double3_t rw11 = { p1->x - pw1->x, p1->y - pw1->y, p1->z - pw1->z };
    const double rw00w = rw00.x * rw00.x + rw00.y * rw00.y + rw00.z * rw00.z;
    const double rw11w = rw11.x * rw11.x + rw11.y * rw11.y + rw11.z * rw11.z;
    if (!gtoint__ecp_type2_spherical_factor_database_array__resize(&(itg->ecp.s), ng0 + ng1)) return GTOINT_ERROR_MEMORY;
    const size_t ng01 = ng0 * ng1;
    const size_t ng01c = ng01 * ngc;
    if (!gtoint__cache__reset(&(itg->c), ng01c)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_array__resize(&(itg->w), ng01 * NVAR)) return GTOINT_ERROR_MEMORY;
    size_t ivar = 0;
    double *const e01 = itg->w.p + ng01 * ivar++;
    assert(ivar == NVAR);
    for (size_t i0 = 0; i0 < ng0; i0++) {
        const double g0cw = g0[i0] + gw0;
        const double gr0 = 1.0 / g0cw;
        const double3_t q0 = { (g0[i0] * p0->x + gw0 * pw0->x) * gr0, (g0[i0] * p0->y + gw0 * pw0->y) * gr0, (g0[i0] * p0->z + gw0 * pw0->z) * gr0 };
        const double3_t s0c = { pc->x - q0.x, pc->y - q0.y, pc->z - q0.z };
        const double s0cw = s0c.x * s0c.x + s0c.y * s0c.y + s0c.z * s0c.z;
        if (s0cw > 0.0) {
            const double r = 1.0 / sqrt(s0cw);
            gtoint__ecp_type2_spherical_factor_database__reset(&(itg->ecp.s.p[i0]), double3__new(s0c.x * r, s0c.y * r, s0c.z * r));
        }
        const double e0 = exp(-g0cw * s0cw);
        for (size_t i1 = 0; i1 < ng1; i1++) {
            e01[i0 + ng0 * i1] = e0;
        }
    }
    for (size_t i1 = 0; i1 < ng1; i1++) {
        const double g1cw = g1[i1] + gw1;
        const double gr1 = 1.0 / g1cw;
        const double3_t q1 = { (g1[i1] * p1->x + gw1 * pw1->x) * gr1, (g1[i1] * p1->y + gw1 * pw1->y) * gr1, (g1[i1] * p1->z + gw1 * pw1->z) * gr1 };
        const double3_t s1c = { pc->x - q1.x, pc->y - q1.y, pc->z - q1.z };
        const double s1cw = s1c.x * s1c.x + s1c.y * s1c.y + s1c.z * s1c.z;
        if (s1cw > 0.0) {
            const double r = 1.0 / sqrt(s1cw);
            gtoint__ecp_type2_spherical_factor_database__reset(&(itg->ecp.s.p[ng0 + i1]), double3__new(s1c.x * r, s1c.y * r, s1c.z * r));
        }
        const double e1 = exp(-g1cw * s1cw);
        for (size_t i0 = 0; i0 < ng0; i0++) {
            e01[i0 + ng0 * i1] *= e1;
        }
    }
    size_t io = 0;
    for (size_t id = 0; id < nd; id++) {
    for (size_t ia1 = 0; ia1 < na1; ia1++) {
    for (size_t ia0 = 0; ia0 < na0; ia0++) {
        itg->v.p[io] = 0.0;
        gtoint__stack__reset(&(itg->s), ng01c);
        {
            stack_index_t t = { 0 };
            t.i.ecp2_w.a[0] = a0[ia0];
            t.i.ecp2_w.a[1] = a1[ia1];
            t.i.ecp2_w.d[0] = d0[id];
            t.i.ecp2_w.d[1] = d1[id];
            t.i.ecp2_w.d[2] = dc[id];
            t.i.ecp2_w.r = rc - 2;
            t.o = STACK_VOID_INDEX;
            t.b = false;
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY;
            size_t it;
            gtoint__stack__write(&(itg->s), &t, &it);
            double *const c = gtoint__stack__coefficients(&(itg->s), it);
            double *const v = gtoint__stack__integrals(&(itg->s), it);
            size_t i = 0;
            for (size_t i1 = 0; i1 < ng1; i1++) {
            for (size_t i0 = 0; i0 < ng0; i0++) {
            for (size_t ic = 0; ic < ngc; ic++, i++) {
                c[i] = c0[i0 + ng0 * ia0] * c1[i1 + ng1 * ia1] * cc[ic];
                v[i] = 0.0;
            }
            }
            }
        }
        while (!gtoint__stack__is_empty(&(itg->s))) {
            stack_index_t s;
            size_t is;
            gtoint__stack__read(&(itg->s), &s, &is);
            if (s.b) {
                const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                const double *const v = gtoint__stack__integrals(&(itg->s), is);
                double *w;
                if (!gtoint__cache__reference_to_store_ecp2_w(&(itg->c), &s.i, &w)) return GTOINT_ERROR_MEMORY;
                for (size_t i = 0; i < ng01c; i++) w[i] = v[i];
                if (s.o == STACK_VOID_INDEX) {
                    double o = 0.0;
                    for (size_t i = 0; i < ng01c; i++) o += c[i] * v[i];
                    itg->v.p[io] += o;
                }
                else {
                    double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                    for (size_t i = 0; i < ng01c; i++) o[i] += c[i] * v[i];
                }
                gtoint__stack__pop(&(itg->s));
            }
            else {
                const double *w;
                if (gtoint__cache__reference_to_fetch_ecp2_w(&(itg->c), &s.i, &w)) {
                    const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                    if (s.o == STACK_VOID_INDEX) {
                        double o = 0.0;
                        for (size_t i = 0; i < ng01c; i++) o += c[i] * w[i];
                        itg->v.p[io] += o;
                    }
                    else {
                        double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                        for (size_t i = 0; i < ng01c; i++) o[i] += c[i] * w[i];
                    }
                    gtoint__stack__pop(&(itg->s));
                }
                else {
                    {
                        s.b = true;
                        gtoint__stack__write(&(itg->s), &s, &is);
                        double *const v = gtoint__stack__integrals(&(itg->s), is);
                        for (size_t i = 0; i < ng01c; i++) v[i] = 0.0;
                    }
                    if (s.i.ecp2_w.d[0].x > 0) {
                        s.i.ecp2_w.d[0].x--;
                        EXPAND_VRR_0(x);
                    }
                    else if (s.i.ecp2_w.d[0].y > 0) {
                        s.i.ecp2_w.d[0].y--;
                        EXPAND_VRR_0(y);
                    }
                    else if (s.i.ecp2_w.d[0].z > 0) {
                        s.i.ecp2_w.d[0].z--;
                        EXPAND_VRR_0(z);
                    }
                    else if (s.i.ecp2_w.d[1].x > 0) {
                        s.i.ecp2_w.d[1].x--;
                        EXPAND_VRR_1(x);
                    }
                    else if (s.i.ecp2_w.d[1].y > 0) {
                        s.i.ecp2_w.d[1].y--;
                        EXPAND_VRR_1(y);
                    }
                    else if (s.i.ecp2_w.d[1].z > 0) {
                        s.i.ecp2_w.d[1].z--;
                        EXPAND_VRR_1(z);
                    }
                    else if (s.i.ecp2_w.dw[0].x > 0) {
                        s.i.ecp2_w.dw[0].x--;
                        EXPAND_VRR_W0(x);
                    }
                    else if (s.i.ecp2_w.dw[0].y > 0) {
                        s.i.ecp2_w.dw[0].y--;
                        EXPAND_VRR_W0(y);
                    }
                    else if (s.i.ecp2_w.dw[0].z > 0) {
                        s.i.ecp2_w.dw[0].z--;
                        EXPAND_VRR_W0(z);
                    }
                    else if (s.i.ecp2_w.dw[1].x > 0) {
                        s.i.ecp2_w.dw[1].x--;
                        EXPAND_VRR_W1(x);
                    }
                    else if (s.i.ecp2_w.dw[1].y > 0) {
                        s.i.ecp2_w.dw[1].y--;
                        EXPAND_VRR_W1(y);
                    }
                    else if (s.i.ecp2_w.dw[1].z > 0) {
                        s.i.ecp2_w.dw[1].z--;
                        EXPAND_VRR_W1(z);
                    }
                    else if (s.i.ecp2_w.d[2].x > 0) {
                        s.i.ecp2_w.d[2].x--;
                        EXPAND_VRR_C(x);
                    }
                    else if (s.i.ecp2_w.d[2].y > 0) {
                        s.i.ecp2_w.d[2].y--;
                        EXPAND_VRR_C(y);
                    }
                    else if (s.i.ecp2_w.d[2].z > 0) {
                        s.i.ecp2_w.d[2].z--;
                        EXPAND_VRR_C(z);
                    }
                    else {
                        double *v;
                        if (!gtoint__cache__reference_to_store_ecp2_w(&(itg->c), &s.i, &v)) return GTOINT_ERROR_MEMORY;
                        size_t i01 = 0;
                        for (size_t i1 = 0; i1 < ng1; i1++) {
                            const double g1cw = g1[i1] + gw1;
                            const double gr1 = 1.0 / g1cw;
                            const double3_t q1 = { (g1[i1] * p1->x + gw1 * pw1->x) * gr1, (g1[i1] * p1->y + gw1 * pw1->y) * gr1, (g1[i1] * p1->z + gw1 * pw1->z) * gr1 };
                            const double3_t s1c = { pc->x - q1.x, pc->y - q1.y, pc->z - q1.z };
                            const double s1cw = s1c.x * s1c.x + s1c.y * s1c.y + s1c.z * s1c.z;
                        for (size_t i0 = 0; i0 < ng0; i0++) {
                            const double g0cw = g0[i0] + gw0;
                            const double gr0 = 1.0 / g0cw;
                            const double3_t q0 = { (g0[i0] * p0->x + gw0 * pw0->x) * gr0, (g0[i0] * p0->y + gw0 * pw0->y) * gr0, (g0[i0] * p0->z + gw0 * pw0->z) * gr0 };
                            const double3_t s0c = { pc->x - q0.x, pc->y - q0.y, pc->z - q0.z };
                            const double s0cw = s0c.x * s0c.x + s0c.y * s0c.y + s0c.z * s0c.z;
                            const double cw = exp(-(g0[i0] * gw0 * gr0 * rw00w + g1[i1] * gw1 * gr1 * rw11w));
                            if (s0cw > 0.0 && s1cw > 0.0) {
                                const int h0 = s.i.ecp2_w.k[0].x + s.i.ecp2_w.k[0].y + s.i.ecp2_w.k[0].z;
                                const int h1 = s.i.ecp2_w.k[1].x + s.i.ecp2_w.k[1].y + s.i.ecp2_w.k[1].z;
                                for (size_t ic = 0; ic < ngc; ic++) {
                                    { /* screening by approximated estimation */
                                        /* Reference:
                                         * "Efficient Method for Calculating Effective Core Potential Integrals"
                                         *  Simon C. McKenzie, Evgeny Epifanovsky, Giuseppe M. J. Barca, Andrew T. B. Gilbert, and Peter M. W. Gill
                                         *  J. Phys. Chem. A 2018, 122, 3066-3075
                                         */
                                        const int a0s = s.i.ecp2_w.a[0].x + s.i.ecp2_w.a[0].y + s.i.ecp2_w.a[0].z;
                                        const int a1s = s.i.ecp2_w.a[1].x + s.i.ecp2_w.a[1].y + s.i.ecp2_w.a[1].z;
                                        const double g0c = g0cw + gc[ic];
                                        const double g1c = g1cw + gc[ic];
                                        const double cg0 = a0s * g0c * g0c / (2.0 * g0cw * (gc[ic] * gc[ic] * s0cw + a0s * g0c));
                                        const double cg1 = a1s * g1c * g1c / (2.0 * g1cw * (gc[ic] * gc[ic] * s1cw + a1s * g1c));
                                        const double ge0 = (1.0 - cg0) * g0cw;
                                        const double ge1 = (1.0 - cg1) * g1cw;
                                        const double ges = 1.0 / (ge0 + ge1 + gc[ic]);
                                        const double t = 2.0 * ge0 * ge1 * sqrt(s0cw * s1cw) * ges;
                                        if (
                                            (2 * ac + 1) * (2 * ac + 1) *
                                            sqrt(power_(a0s / (D_2E * g0cw * cg0), a0s)) *
                                            sqrt(power_(a1s / (D_2E * g1cw * cg1), a1s)) *
                                            (D_RPI3 / (gc[ic] * sqrt(gc[ic]))) *
                                            exp(ge0 * s0cw * (ge0 * ges - 1.0) + ge1 * s1cw * (ge1 * ges - 1.0)) *
                                            ((t <= 1.0) ? D_SINH1 : exp(t) / (2.0 * t)) * cw <= itg->cut
                                        ) { v[i01 + ic] = 0.0; continue; }
                                    }
                                    gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), s0c, s1c, g0cw, g1cw, gc[ic]);
                                    double v0x = 0.0;
                                    for (int k0x = 0; k0x <= s.i.ecp2_w.a[0].x; k0x++) {
                                    double v0y = 0.0;
                                    for (int k0y = 0; k0y <= s.i.ecp2_w.a[0].y; k0y++) {
                                    double v0z = 0.0;
                                    for (int k0z = 0; k0z <= s.i.ecp2_w.a[0].z; k0z++) {
                                    double v1x = 0.0;
                                    for (int k1x = 0; k1x <= s.i.ecp2_w.a[1].x; k1x++) {
                                    double v1y = 0.0;
                                    for (int k1y = 0; k1y <= s.i.ecp2_w.a[1].y; k1y++) {
                                    double v1z = 0.0;
                                    for (int k1z = 0; k1z <= s.i.ecp2_w.a[1].z; k1z++) {
                                    const int k0 = k0x + k0y + k0z;
                                    const int k1 = k1x + k1y + k1z;
                                    double v01 = 0.0;
                                    for (int l0 = 0; l0 <= ac + k0 + h0; l0++) {
                                        size_t sh0, si0;
                                        if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[i0]), &(itg->ecp.h), l0, &sh0)) return GTOINT_ERROR_MEMORY;
                                        {
                                            const ecp_type2_angular_integral_index_t t = {
                                                { l0, ac }, { k0x + s.i.ecp2_w.k[0].x, k0y + s.i.ecp2_w.k[0].y, k0z + s.i.ecp2_w.k[0].z }
                                            };
                                            if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                                        }
                                    for (int l1 = 0; l1 <= ac + k1 + h1; l1++) {
                                        size_t sh1, si1;
                                        if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[ng0 + i1]), &(itg->ecp.h), l1, &sh1)) return GTOINT_ERROR_MEMORY;
                                        {
                                            const ecp_type2_angular_integral_index_t t = {
                                                { l1, ac }, { k1x + s.i.ecp2_w.k[1].x, k1y + s.i.ecp2_w.k[1].y, k1z + s.i.ecp2_w.k[1].z }
                                            };
                                            if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                                        }
                                        const double *const th0 = itg->ecp.s.p[i0].a.p[sh0].v.p;
                                        const double *const th1 = itg->ecp.s.p[ng0 + i1].a.p[sh1].v.p;
                                        const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                        const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                        double vm = 0.0;
                                        for (int m = 0; m <= 2 * ac; m++) {
                                            double vj0 = 0.0;
                                            for (int j0 = 0; j0 <= 2 * l0; j0++) {
                                                vj0 += th0[j0] * ti0[j0 + (2 * l0 + 1) * m];
                                            }
                                            double vj1 = 0.0;
                                            for (int j1 = 0; j1 <= 2 * l1; j1++) {
                                                vj1 += th1[j1] * ti1[j1 + (2 * l1 + 1) * m];
                                            }
                                            vm += vj0 * vj1;
                                        }
                                        if (vm != 0.0) {
                                            const ecp_type2_radial_integral_index_t t = { { l0, l1 }, s.i.ecp2_w.r + k0 + k1 };
                                            double q;
                                            if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                            v01 += vm * (((l0 + l1) & 1) ? -q : q);
                                        }
                                    }
                                    }
                                    v1z += v01 * binomial_(s.i.ecp2_w.a[1].z, k1z) * power_(r1c.z, s.i.ecp2_w.a[1].z - k1z);
                                    }
                                    v1y += v1z * binomial_(s.i.ecp2_w.a[1].y, k1y) * power_(r1c.y, s.i.ecp2_w.a[1].y - k1y);
                                    }
                                    v1x += v1y * binomial_(s.i.ecp2_w.a[1].x, k1x) * power_(r1c.x, s.i.ecp2_w.a[1].x - k1x);
                                    }
                                    v0z += v1x * binomial_(s.i.ecp2_w.a[0].z, k0z) * power_(r0c.z, s.i.ecp2_w.a[0].z - k0z);
                                    }
                                    v0y += v0z * binomial_(s.i.ecp2_w.a[0].y, k0y) * power_(r0c.y, s.i.ecp2_w.a[0].y - k0y);
                                    }
                                    v0x += v0y * binomial_(s.i.ecp2_w.a[0].x, k0x) * power_(r0c.x, s.i.ecp2_w.a[0].x - k0x);
                                    }
                                    v[i01 + ic] = v0x * e01[i0 + ng0 * i1] * D_16PI2 * cw;
                                }
                            }
                            else if (s0cw > 0.0) {
                                const int h0 = s.i.ecp2_w.k[0].x + s.i.ecp2_w.k[0].y + s.i.ecp2_w.k[0].z;
                                const int k1 = s.i.ecp2_w.a[1].x + s.i.ecp2_w.a[1].y + s.i.ecp2_w.a[1].z;
                                size_t si1;
                                {
                                    const ecp_type2_angular_integral_index_t t = {
                                        { 0, ac }, { s.i.ecp2_w.a[1].x + s.i.ecp2_w.k[1].x, s.i.ecp2_w.a[1].y + s.i.ecp2_w.k[1].y, s.i.ecp2_w.a[1].z + s.i.ecp2_w.k[1].z }
                                    };
                                    if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                                }
                                for (size_t ic = 0; ic < ngc; ic++) {
                                    { /* screening by approximated estimation */
                                        const int a0s = s.i.ecp2_w.a[0].x + s.i.ecp2_w.a[0].y + s.i.ecp2_w.a[0].z;
                                        const double g0c = g0cw + gc[ic];
                                        const double cg0 = a0s * g0c * g0c / (2.0 * g0cw * (gc[ic] * gc[ic] * s0cw + a0s * g0c));
                                        const double ge0 = (1.0 - cg0) * g0cw;
                                        const double t = (a0s + 3) / g0cw;
                                        if (
                                            (2 * ac + 1) * (2 * ac + 1) *
                                            sqrt(power_(a0s * t / ((D_2E * D_2E) * g0cw* cg0), a0s) * t * D_PI3D) * t *
                                            exp(-ge0 * gc[ic] * s0cw / (ge0 + gc[ic])) * cw <= itg->cut
                                        ) { v[i01 + ic] = 0.0; continue; }
                                    }
                                    gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), s0c, s1c, g0cw, g1cw, gc[ic]);
                                    double v0x = 0.0;
                                    for (int k0x = 0; k0x <= s.i.ecp2_w.a[0].x; k0x++) {
                                    double v0y = 0.0;
                                    for (int k0y = 0; k0y <= s.i.ecp2_w.a[0].y; k0y++) {
                                    double v0z = 0.0;
                                    for (int k0z = 0; k0z <= s.i.ecp2_w.a[0].z; k0z++) {
                                    double v1x = 0.0;
                                    for (int k1x = 0; k1x <= s.i.ecp2_w.a[1].x; k1x++) {
                                    double v1y = 0.0;
                                    for (int k1y = 0; k1y <= s.i.ecp2_w.a[1].y; k1y++) {
                                    double v1z = 0.0;
                                    for (int k1z = 0; k1z <= s.i.ecp2_w.a[1].z; k1z++) {
                                    const int k0 = k0x + k0y + k0z;
                                    double v01 = 0.0;
                                    for (int l0 = 0; l0 <= ac + k0 + h0; l0++) {
                                        size_t sh0, si0;
                                        if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[i0]), &(itg->ecp.h), l0, &sh0)) return GTOINT_ERROR_MEMORY;
                                        {
                                            const ecp_type2_angular_integral_index_t t = {
                                                { l0, ac }, { k0x + s.i.ecp2_w.k[0].x, k0y + s.i.ecp2_w.k[0].y, k0z + s.i.ecp2_w.k[0].z }
                                            };
                                            if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                                        }
                                        const double *const th0 = itg->ecp.s.p[i0].a.p[sh0].v.p;
                                        const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                        const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                        double vm = 0.0;
                                        for (int m = 0; m <= 2 * ac; m++) {
                                            double vj0 = 0.0;
                                            for (int j0 = 0; j0 <= 2 * l0; j0++) {
                                                vj0 += th0[j0] * ti0[j0 + (2 * l0 + 1) * m];
                                            }
                                            vm += vj0 * ti1[m];
                                        }
                                        if (vm != 0.0) {
                                            const ecp_type2_radial_integral_index_t t = { { l0, 0 }, s.i.ecp2_w.r + k0 + k1 };
                                            double q;
                                            if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                            v01 += vm * ((l0 & 1) ? -q : q);
                                        }
                                    }
                                    v1z += v01 * binomial_(s.i.ecp2_w.a[1].z, k1z) * power_(r1c.z, s.i.ecp2_w.a[1].z - k1z);
                                    }
                                    v1y += v1z * binomial_(s.i.ecp2_w.a[1].y, k1y) * power_(r1c.y, s.i.ecp2_w.a[1].y - k1y);
                                    }
                                    v1x += v1y * binomial_(s.i.ecp2_w.a[1].x, k1x) * power_(r1c.x, s.i.ecp2_w.a[1].x - k1x);
                                    }
                                    v0z += v1x * binomial_(s.i.ecp2_w.a[0].z, k0z) * power_(r0c.z, s.i.ecp2_w.a[0].z - k0z);
                                    }
                                    v0y += v0z * binomial_(s.i.ecp2_w.a[0].y, k0y) * power_(r0c.y, s.i.ecp2_w.a[0].y - k0y);
                                    }
                                    v0x += v0y * binomial_(s.i.ecp2_w.a[0].x, k0x) * power_(r0c.x, s.i.ecp2_w.a[0].x - k0x);
                                    }
                                    v[i01 + ic] = v0x * e01[i0 + ng0 * i1] * D_8RPI3 * cw;
                                }
                            }
                            else if (s1cw > 0.0) {
                                const int h1 = s.i.ecp2_w.k[1].x + s.i.ecp2_w.k[1].y + s.i.ecp2_w.k[1].z;
                                const int k0 = s.i.ecp2_w.a[0].x + s.i.ecp2_w.a[0].y + s.i.ecp2_w.a[0].z;
                                size_t si0;
                                {
                                    const ecp_type2_angular_integral_index_t t = {
                                        { 0, ac }, { s.i.ecp2_w.a[0].x + s.i.ecp2_w.k[0].x, s.i.ecp2_w.a[0].y + s.i.ecp2_w.k[0].y, s.i.ecp2_w.a[0].z + s.i.ecp2_w.k[0].z }
                                    };
                                    if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                                }
                                for (size_t ic = 0; ic < ngc; ic++) {
                                    { /* screening by approximated estimation */
                                        const int a1s = s.i.ecp2_w.a[1].x + s.i.ecp2_w.a[1].y + s.i.ecp2_w.a[1].z;
                                        const double g1c = g1cw + gc[ic];
                                        const double cg1 = a1s * g1c * g1c / (2.0 * g1cw * (gc[ic] * gc[ic] * s1cw + a1s * g1c));
                                        const double ge1 = (1.0 - cg1) * g1cw;
                                        const double t = (a1s + 3) / g1cw;
                                        if (
                                            (2 * ac + 1) * (2 * ac + 1) *
                                            sqrt(power_(a1s * t / ((D_2E * D_2E) * g1cw * cg1), a1s) * t * D_PI3D) * t *
                                            exp(-ge1 * gc[ic] * s1cw / (ge1 + gc[ic])) * cw <= itg->cut
                                        ) { v[i01 + ic] = 0.0; continue; }
                                    }
                                    gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), s0c, s1c, g0cw, g1cw, gc[ic]);
                                    double v0x = 0.0;
                                    for (int k0x = 0; k0x <= s.i.ecp2_w.a[0].x; k0x++) {
                                    double v0y = 0.0;
                                    for (int k0y = 0; k0y <= s.i.ecp2_w.a[0].y; k0y++) {
                                    double v0z = 0.0;
                                    for (int k0z = 0; k0z <= s.i.ecp2_w.a[0].z; k0z++) {
                                    double v1x = 0.0;
                                    for (int k1x = 0; k1x <= s.i.ecp2_w.a[1].x; k1x++) {
                                    double v1y = 0.0;
                                    for (int k1y = 0; k1y <= s.i.ecp2_w.a[1].y; k1y++) {
                                    double v1z = 0.0;
                                    for (int k1z = 0; k1z <= s.i.ecp2_w.a[1].z; k1z++) {
                                    const int k1 = k1x + k1y + k1z;
                                    double v01 = 0.0;
                                    for (int l1 = 0; l1 <= ac + k1 + h1; l1++) {
                                        size_t sh1, si1;
                                        if (!gtoint__ecp_type2_spherical_factor_database__index(&(itg->ecp.s.p[ng0 + i1]), &(itg->ecp.h), l1, &sh1)) return GTOINT_ERROR_MEMORY;
                                        {
                                            const ecp_type2_angular_integral_index_t t = {
                                                { l1, ac }, { k1x + s.i.ecp2_w.k[1].x, k1y + s.i.ecp2_w.k[1].y, k1z + s.i.ecp2_w.k[1].z }
                                            };
                                            if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                                        }
                                        const double *const th1 = itg->ecp.s.p[ng0 + i1].a.p[sh1].v.p;
                                        const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                        const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                        double vm = 0.0;
                                        for (int m = 0; m <= 2 * ac; m++) {
                                            double vj1 = 0.0;
                                            for (int j1 = 0; j1 <= 2 * l1; j1++) {
                                                vj1 += th1[j1] * ti1[j1 + (2 * l1 + 1) * m];
                                            }
                                            vm += vj1 * ti0[m];
                                        }
                                        if (vm != 0.0) {
                                            const ecp_type2_radial_integral_index_t t = { { 0, l1 }, s.i.ecp2_w.r + k0 + k1 };
                                            double q;
                                            if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                            v01 += vm * ((l1 & 1) ? -q : q);
                                        }
                                    }
                                    v1z += v01 * binomial_(s.i.ecp2_w.a[1].z, k1z) * power_(r1c.z, s.i.ecp2_w.a[1].z - k1z);
                                    }
                                    v1y += v1z * binomial_(s.i.ecp2_w.a[1].y, k1y) * power_(r1c.y, s.i.ecp2_w.a[1].y - k1y);
                                    }
                                    v1x += v1y * binomial_(s.i.ecp2_w.a[1].x, k1x) * power_(r1c.x, s.i.ecp2_w.a[1].x - k1x);
                                    }
                                    v0z += v1x * binomial_(s.i.ecp2_w.a[0].z, k0z) * power_(r0c.z, s.i.ecp2_w.a[0].z - k0z);
                                    }
                                    v0y += v0z * binomial_(s.i.ecp2_w.a[0].y, k0y) * power_(r0c.y, s.i.ecp2_w.a[0].y - k0y);
                                    }
                                    v0x += v0y * binomial_(s.i.ecp2_w.a[0].x, k0x) * power_(r0c.x, s.i.ecp2_w.a[0].x - k0x);
                                    }
                                    v[i01 + ic] = v0x * e01[i0 + ng0 * i1] * D_8RPI3 * cw;
                                }
                            }
                            else {
                                const int k0 = s.i.ecp2_w.a[0].x + s.i.ecp2_w.a[0].y + s.i.ecp2_w.a[0].z;
                                const int k1 = s.i.ecp2_w.a[1].x + s.i.ecp2_w.a[1].y + s.i.ecp2_w.a[1].z;
                                size_t si0;
                                {
                                    const ecp_type2_angular_integral_index_t t = {
                                        { 0, ac }, { s.i.ecp2_w.a[0].x + s.i.ecp2_w.k[0].x, s.i.ecp2_w.a[0].y + s.i.ecp2_w.k[0].y, s.i.ecp2_w.a[0].z + s.i.ecp2_w.k[0].z }
                                    };
                                    if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si0)) return GTOINT_ERROR_MEMORY;
                                }
                                size_t si1;
                                {
                                    const ecp_type2_angular_integral_index_t t = {
                                        { 0, ac }, { s.i.ecp2_w.a[1].x + s.i.ecp2_w.k[1].x, s.i.ecp2_w.a[1].y + s.i.ecp2_w.k[1].y, s.i.ecp2_w.a[1].z + s.i.ecp2_w.k[1].z }
                                    };
                                    if (!gtoint__ecp_type2_angular_integral_database__index(&(itg->ecp.a), &(itg->ecp.h), &t, &si1)) return GTOINT_ERROR_MEMORY;
                                }
                                for (size_t ic = 0; ic < ngc; ic++) {
                                    gtoint__ecp_type2_radial_integral_database__reset(&(itg->ecp.r), s0c, s1c, g0cw, g1cw, gc[ic]);
                                    const double *const ti0 = itg->ecp.a.a.p[si0].v.p;
                                    const double *const ti1 = itg->ecp.a.a.p[si1].v.p;
                                    double v0x = 0.0;
                                    for (int k0x = 0; k0x <= s.i.ecp2_w.a[0].x; k0x++) {
                                    double v0y = 0.0;
                                    for (int k0y = 0; k0y <= s.i.ecp2_w.a[0].y; k0y++) {
                                    double v0z = 0.0;
                                    for (int k0z = 0; k0z <= s.i.ecp2_w.a[0].z; k0z++) {
                                    double v1x = 0.0;
                                    for (int k1x = 0; k1x <= s.i.ecp2_w.a[1].x; k1x++) {
                                    double v1y = 0.0;
                                    for (int k1y = 0; k1y <= s.i.ecp2_w.a[1].y; k1y++) {
                                    double v1z = 0.0;
                                    for (int k1z = 0; k1z <= s.i.ecp2_w.a[1].z; k1z++) {
                                    double v01 = 0.0;
                                    {
                                        double vm = 0.0;
                                        for (int m = 0; m <= 2 * ac; m++) {
                                            vm += ti0[m] * ti1[m];
                                        }
                                        if (vm != 0.0) {
                                            const ecp_type2_radial_integral_index_t t = { { 0, 0 }, s.i.ecp2_w.r + k0 + k1 };
                                            double q;
                                            if (!gtoint__ecp_type2_radial_integral_database__fetch(&(itg->ecp.r), &t, &q)) return GTOINT_ERROR_MEMORY;
                                            v01 += vm * q;
                                        }
                                    }
                                    v1z += v01 * binomial_(s.i.ecp2_w.a[1].z, k1z) * power_(r1c.z, s.i.ecp2_w.a[1].z - k1z);
                                    }
                                    v1y += v1z * binomial_(s.i.ecp2_w.a[1].y, k1y) * power_(r1c.y, s.i.ecp2_w.a[1].y - k1y);
                                    }
                                    v1x += v1y * binomial_(s.i.ecp2_w.a[1].x, k1x) * power_(r1c.x, s.i.ecp2_w.a[1].x - k1x);
                                    }
                                    v0z += v1x * binomial_(s.i.ecp2_w.a[0].z, k0z) * power_(r0c.z, s.i.ecp2_w.a[0].z - k0z);
                                    }
                                    v0y += v0z * binomial_(s.i.ecp2_w.a[0].y, k0y) * power_(r0c.y, s.i.ecp2_w.a[0].y - k0y);
                                    }
                                    v0x += v0y * binomial_(s.i.ecp2_w.a[0].x, k0x) * power_(r0c.x, s.i.ecp2_w.a[0].x - k0x);
                                    }
                                    v[i01 + ic] = v0x * e01[i0 + ng0 * i1] * D_4PI * cw;
                                }
                            }
                            i01 += ngc;
                        }
                        }
                        const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                        if (s.o == STACK_VOID_INDEX) {
                            double o = 0.0;
                            for (size_t i = 0; i < ng01c; i++) o += c[i] * v[i];
                            itg->v.p[io] += o;
                        }
                        else {
                            double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                            for (size_t i = 0; i < ng01c; i++) o[i] += c[i] * v[i];
                        }
                        gtoint__stack__pop(&(itg->s));
                    }
                }
            }
        }
        io++;
    }
    }
    }
    return GTOINT_ERROR_OK;
#undef NVAR
#undef EXPAND_VRR_0
#undef EXPAND_VRR_1
#undef EXPAND_VRR_W0
#undef EXPAND_VRR_W1
#undef EXPAND_VRR_C
}
