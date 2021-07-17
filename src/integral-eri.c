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

#include "gtoint-private.h"

#include "function.h"

#include <math.h>
#include <assert.h>

#define D_PI    3.1415926535897932384626433832795
#define D_2DRPI 1.1283791670955125738961589031215 /* 2/sqrt(PI) */

static gtoint_error_t compute_electron_repulsion_integrals_(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    const double3_t *p2, size_t na2, const int3_t *a2, size_t ng2, const double *g2, const double *c2,
    const double3_t *p3, size_t na3, const int3_t *a3, size_t ng3, const double *g3, const double *c3,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *d2, const int3_t *d3
) {
#define NVAR (6 * 2 + 2 + 14 * 3 * 8)
#define INIT_VRR_COEFFS_A(out, g0_, g2_, g3_, g01_, p0_, p01_, p0123_, g0123_, h0123_) /* referrence variable: i */ \
    { \
        const double vg0 = (g0_); \
        const double vg2 = (g2_); \
        const double vg3 = (g3_); \
        const double vg01 = (g01_); \
        const double vp0 = (p0_); \
        const double vp01 = (p01_); \
        const double vp0123 = (p0123_); \
        const double vg01h = 0.5 * vg01; \
        const double vh0 = vg0 * vg01; \
        const double vh1 = 1.0 - vh0; \
        const double vg012 = -vg01 * (h0123_); \
        const double vh012 = vg01h * vg012; \
        const double vg0123h = 0.5 * (g0123_); \
        out[ 0][i] = vp01 - vp0; \
        out[ 1][i] = -vg012 * vp0123; \
        out[ 2][i] = -vh1; \
        out[ 3][i] = vh0 * vg012; \
        out[ 4][i] = vg01h; \
        out[ 5][i] = vh012; \
        out[ 6][i] = vh1; \
        out[ 7][i] = vh1 * vg012; \
        out[ 8][i] = vg01h; \
        out[ 9][i] = vh012; \
        out[10][i] = vg2 * (g0123_); \
        out[11][i] = vg0123h; \
        out[12][i] = vg3 * (g0123_); \
        out[13][i] = vg0123h; \
    }
#define INIT_VRR_COEFFS_D(out, g0_, g2_, g3_, g01_, f01_, p0_, p01_, p0123_, g0123_, h0123_) /* referrence variable: i */ \
    { \
        const double vg0 = (g0_); \
        const double vg2 = (g2_); \
        const double vg3 = (g3_); \
        const double vg01 = (g01_); \
        const double vf01 = (f01_); \
        const double vp0 = (p0_); \
        const double vp01 = (p01_); \
        const double vp0123 = (p0123_); \
        const double vg0d = 2.0 * vg0; \
        const double vf01d = -2.0 * vf01; \
        const double vh0 = vg0 * vg01; \
        const double vh012 = -vh0 * vg01 * (h0123_); \
        const double vh0123 = vg0 * (g0123_); \
        out[ 0][i] = vg0d * (vp01 - vp0); \
        out[ 1][i] = 2.0 * vh0 * (h0123_) * vp0123; \
        out[ 2][i] = vf01d; \
        out[ 3][i] = vg0d * vh012; \
        out[ 4][i] = vh0 - 1.0; \
        out[ 5][i] = vh012; \
        out[ 6][i] = -vf01d; \
        out[ 7][i] = vf01d * vg01 * (h0123_); \
        out[ 8][i] = vh0; \
        out[ 9][i] = vh012; \
        out[10][i] = 2.0 * vg2 * vh0123; \
        out[11][i] = vh0123; \
        out[12][i] = 2.0 * vg3 * vh0123; \
        out[13][i] = vh0123; \
    }
#define EXPAND_VRR(coe, xyz, i0, i1, i2, i3) /* referrence variables: itg, e, s, is, ng0123 */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng0123; i++) { \
                c[i] = coe[0][i]; \
                v[i] = 0.0; \
            } \
        } \
        { \
            stack_index_t t = s; \
            t.i.eri.m++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng0123; i++) { \
                c[i] = coe[1][i]; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.eri.d[i0].xyz > 0) { \
            { \
                stack_index_t t = s; \
                t.i.eri.d[i0].xyz--; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[2][i] * s.i.eri.d[i0].xyz; \
                    v[i] = 0.0; \
                } \
            } \
            { \
                stack_index_t t = s; \
                t.i.eri.d[i0].xyz--; \
                t.i.eri.m++; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[3][i] * s.i.eri.d[i0].xyz; \
                    v[i] = 0.0; \
                } \
            } \
        } \
        if (s.i.eri.a[i0].xyz > 0) { \
            { \
                stack_index_t t = s; \
                t.i.eri.a[i0].xyz--; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[4][i] * s.i.eri.a[i0].xyz; \
                    v[i] = 0.0; \
                } \
            } \
            { \
                stack_index_t t = s; \
                t.i.eri.a[i0].xyz--; \
                t.i.eri.m++; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[5][i] * s.i.eri.a[i0].xyz; \
                    v[i] = 0.0; \
                } \
            } \
        } \
        if (s.i.eri.d[i1].xyz > 0) { \
            { \
                stack_index_t t = s; \
                t.i.eri.d[i1].xyz--; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[6][i] * s.i.eri.d[i1].xyz; \
                    v[i] = 0.0; \
                } \
            } \
            { \
                stack_index_t t = s; \
                t.i.eri.d[i1].xyz--; \
                t.i.eri.m++; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[7][i] * s.i.eri.d[i1].xyz; \
                    v[i] = 0.0; \
                } \
            } \
        } \
        if (s.i.eri.a[i1].xyz > 0) { \
            { \
                stack_index_t t = s; \
                t.i.eri.a[i1].xyz--; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[8][i] * s.i.eri.a[i1].xyz; \
                    v[i] = 0.0; \
                } \
            } \
            { \
                stack_index_t t = s; \
                t.i.eri.a[i1].xyz--; \
                t.i.eri.m++; \
                if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
                size_t it; \
                gtoint__stack__write(&(itg->s), &t, &it); \
                double *const c = gtoint__stack__coefficients(&(itg->s), it); \
                double *const v = gtoint__stack__integrals(&(itg->s), it); \
                for (size_t i = 0; i < ng0123; i++) { \
                    c[i] = coe[9][i] * s.i.eri.a[i1].xyz; \
                    v[i] = 0.0; \
                } \
            } \
        } \
        if (s.i.eri.d[i2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.eri.d[i2].xyz--; \
            t.i.eri.m++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng0123; i++) { \
                c[i] = coe[10][i] * s.i.eri.d[i2].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.eri.a[i2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.eri.a[i2].xyz--; \
            t.i.eri.m++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng0123; i++) { \
                c[i] = coe[11][i] * s.i.eri.a[i2].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.eri.d[i3].xyz > 0) { \
            stack_index_t t = s; \
            t.i.eri.d[i3].xyz--; \
            t.i.eri.m++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng0123; i++) { \
                c[i] = coe[12][i] * s.i.eri.d[i3].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.eri.a[i3].xyz > 0) { \
            stack_index_t t = s; \
            t.i.eri.a[i3].xyz--; \
            t.i.eri.m++; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng0123; i++) { \
                c[i] = coe[13][i] * s.i.eri.a[i3].xyz; \
                v[i] = 0.0; \
            } \
        } \
    }
    if (!gtoint__double_array__resize(&(itg->v), na0 * na1 * na2 * na3 * nd)) return GTOINT_ERROR_MEMORY;
    const double r01x = p1->x - p0->x;
    const double r01y = p1->y - p0->y;
    const double r01z = p1->z - p0->z;
    const double r01w = r01x * r01x + r01y * r01y + r01z * r01z;
    const double r23x = p3->x - p2->x;
    const double r23y = p3->y - p2->y;
    const double r23z = p3->z - p2->z;
    const double r23w = r23x * r23x + r23y * r23y + r23z * r23z;
    int3_t max_a0 = { 0, 0, 0 };
    int3_t max_a1 = { 0, 0, 0 };
    int3_t max_a2 = { 0, 0, 0 };
    int3_t max_a3 = { 0, 0, 0 };
    int3_t max_d0 = { 0, 0, 0 };
    int3_t max_d1 = { 0, 0, 0 };
    int3_t max_d2 = { 0, 0, 0 };
    int3_t max_d3 = { 0, 0, 0 };
    for (size_t i = 0; i < na0; i++) {
        if (max_a0.x < a0[i].x) max_a0.x = a0[i].x;
        if (max_a0.y < a0[i].y) max_a0.y = a0[i].y;
        if (max_a0.z < a0[i].z) max_a0.z = a0[i].z;
    }
    for (size_t i = 0; i < na1; i++) {
        if (max_a1.x < a1[i].x) max_a1.x = a1[i].x;
        if (max_a1.y < a1[i].y) max_a1.y = a1[i].y;
        if (max_a1.z < a1[i].z) max_a1.z = a1[i].z;
    }
    for (size_t i = 0; i < na2; i++) {
        if (max_a2.x < a2[i].x) max_a2.x = a2[i].x;
        if (max_a2.y < a2[i].y) max_a2.y = a2[i].y;
        if (max_a2.z < a2[i].z) max_a2.z = a2[i].z;
    }
    for (size_t i = 0; i < na3; i++) {
        if (max_a3.x < a3[i].x) max_a3.x = a3[i].x;
        if (max_a3.y < a3[i].y) max_a3.y = a3[i].y;
        if (max_a3.z < a3[i].z) max_a3.z = a3[i].z;
    }
    for (size_t i = 0; i < nd; i++) {
        if (max_d0.x < d0[i].x) max_d0.x = d0[i].x;
        if (max_d0.y < d0[i].y) max_d0.y = d0[i].y;
        if (max_d0.z < d0[i].z) max_d0.z = d0[i].z;
        if (max_d1.x < d1[i].x) max_d1.x = d1[i].x;
        if (max_d1.y < d1[i].y) max_d1.y = d1[i].y;
        if (max_d1.z < d1[i].z) max_d1.z = d1[i].z;
        if (max_d2.x < d2[i].x) max_d2.x = d2[i].x;
        if (max_d2.y < d2[i].y) max_d2.y = d2[i].y;
        if (max_d2.z < d2[i].z) max_d2.z = d2[i].z;
        if (max_d3.x < d3[i].x) max_d3.x = d3[i].x;
        if (max_d3.y < d3[i].y) max_d3.y = d3[i].y;
        if (max_d3.z < d3[i].z) max_d3.z = d3[i].z;
    }
    const size_t ng01 = ng0 * ng1;
    const size_t ng23 = ng2 * ng3;
    const size_t ng0123 = ng01 * ng23;
    if (!gtoint__cache__reset(&(itg->c), ng0123)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_array__resize(&(itg->w), ng0123 * NVAR)) return GTOINT_ERROR_MEMORY;
    size_t ivar = 0;
    double *const g01   = itg->w.p + ng0123 * ivar++;
    double *const p01x  = itg->w.p + ng0123 * ivar++;
    double *const p01y  = itg->w.p + ng0123 * ivar++;
    double *const p01z  = itg->w.p + ng0123 * ivar++;
    double *const f01   = itg->w.p + ng0123 * ivar++;
    double *const o01   = itg->w.p + ng0123 * ivar++;
    double *const g23   = itg->w.p + ng0123 * ivar++;
    double *const p23x  = itg->w.p + ng0123 * ivar++;
    double *const p23y  = itg->w.p + ng0123 * ivar++;
    double *const p23z  = itg->w.p + ng0123 * ivar++;
    double *const f23   = itg->w.p + ng0123 * ivar++;
    double *const o23   = itg->w.p + ng0123 * ivar++;
    double *const o0123 = itg->w.p + ng0123 * ivar++;
    double *const q0123 = itg->w.p + ng0123 * ivar++;
    double *ca0x[14], *ca0y[14], *ca0z[14];
    double *ca1x[14], *ca1y[14], *ca1z[14];
    double *ca2x[14], *ca2y[14], *ca2z[14];
    double *ca3x[14], *ca3y[14], *ca3z[14];
    double *cd0x[14], *cd0y[14], *cd0z[14];
    double *cd1x[14], *cd1y[14], *cd1z[14];
    double *cd2x[14], *cd2y[14], *cd2z[14];
    double *cd3x[14], *cd3y[14], *cd3z[14];
    for (size_t i = 0; i < 14; i++) ca0x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca0y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca0z[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca1x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca1y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca1z[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca2x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca2y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca2z[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca3x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca3y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) ca3z[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd0x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd0y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd0z[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd1x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd1y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd1z[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd2x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd2y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd2z[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd3x[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd3y[i] = itg->w.p + ng0123 * ivar++;
    for (size_t i = 0; i < 14; i++) cd3z[i] = itg->w.p + ng0123 * ivar++;
    assert(ivar == NVAR);
    for (size_t i1 = 0; i1 < ng1; i1++) {
    for (size_t i0 = 0; i0 < ng0; i0++) {
        const double t_g01 = 1.0 / (g0[i0] + g1[i1]);
        const double t_p01x = (g0[i0] * p0->x + g1[i1] * p1->x) * t_g01;
        const double t_p01y = (g0[i0] * p0->y + g1[i1] * p1->y) * t_g01;
        const double t_p01z = (g0[i0] * p0->z + g1[i1] * p1->z) * t_g01;
        const double t_h01 = D_PI * t_g01;
        const double t_f01 = g0[i0] * g1[i1] * t_g01;
        const double t_o01 = t_h01 * sqrt(t_h01) * exp(-t_f01 * r01w);
        const size_t i01 = i0 + ng0 * i1;
        for (size_t i23 = 0; i23 < ng23; i23++) {
            const size_t i = i01 + ng01 * i23;
            g01[i]  = t_g01;
            p01x[i] = t_p01x;
            p01y[i] = t_p01y;
            p01z[i] = t_p01z;
            f01[i]  = t_f01;
            o01[i]  = t_o01;
        }
    }
    }
    for (size_t i3 = 0; i3 < ng3; i3++) {
    for (size_t i2 = 0; i2 < ng2; i2++) {
        const double t_g23 = 1.0 / (g2[i2] + g3[i3]);
        const double t_p23x = (g2[i2] * p2->x + g3[i3] * p3->x) * t_g23;
        const double t_p23y = (g2[i2] * p2->y + g3[i3] * p3->y) * t_g23;
        const double t_p23z = (g2[i2] * p2->z + g3[i3] * p3->z) * t_g23;
        const double t_h23 = D_PI * t_g23;
        const double t_f23 = g2[i2] * g3[i3] * t_g23;
        const double t_o23 = t_h23 * sqrt(t_h23) * exp(-t_f23 * r23w);
        const size_t i23 = ng01 * (i2 + ng2 * i3);
        for (size_t i01 = 0; i01 < ng01; i01++) {
            const size_t i = i23 + i01;
            g23[i]  = t_g23;
            p23x[i] = t_p23x;
            p23y[i] = t_p23y;
            p23z[i] = t_p23z;
            f23[i]  = t_f23;
            o23[i]  = t_o23;
        }
    }
    }
    for (size_t i3 = 0; i3 < ng3; i3++) {
    for (size_t i2 = 0; i2 < ng2; i2++) {
    for (size_t i1 = 0; i1 < ng1; i1++) {
    for (size_t i0 = 0; i0 < ng0; i0++) {
        const size_t i = i0 + ng0 * (i1 + ng1 * (i2 + ng2 * i3));
        const double g0123 = 1.0 / (g0[i0] + g1[i1] + g2[i2] + g3[i3]);
        const double h0123 = (g0[i0] + g1[i1]) * (g2[i2] + g3[i3]) * g0123;
        const double p0123x = p23x[i] - p01x[i];
        const double p0123y = p23y[i] - p01y[i];
        const double p0123z = p23z[i] - p01z[i];
        o0123[i] = D_2DRPI * sqrt(h0123) * o01[i] * o23[i];
        q0123[i] = h0123 * (p0123x * p0123x + p0123y * p0123y + p0123z * p0123z);
        if (max_a0.x > 0) {
            INIT_VRR_COEFFS_A(ca0x, g0[i0], g2[i2], g3[i3], g01[i], p0->x, p01x[i], p0123x, g0123, h0123);
        }
        if (max_a0.y > 0) {
            INIT_VRR_COEFFS_A(ca0y, g0[i0], g2[i2], g3[i3], g01[i], p0->y, p01y[i], p0123y, g0123, h0123);
        }
        if (max_a0.z > 0) {
            INIT_VRR_COEFFS_A(ca0z, g0[i0], g2[i2], g3[i3], g01[i], p0->z, p01z[i], p0123z, g0123, h0123);
        }
        if (max_a1.x > 0) {
            INIT_VRR_COEFFS_A(ca1x, g1[i1], g3[i3], g2[i2], g01[i], p1->x, p01x[i], p0123x, g0123, h0123);
        }
        if (max_a1.y > 0) {
            INIT_VRR_COEFFS_A(ca1y, g1[i1], g3[i3], g2[i2], g01[i], p1->y, p01y[i], p0123y, g0123, h0123);
        }
        if (max_a1.z > 0) {
            INIT_VRR_COEFFS_A(ca1z, g1[i1], g3[i3], g2[i2], g01[i], p1->z, p01z[i], p0123z, g0123, h0123);
        }
        if (max_a2.x > 0) {
            INIT_VRR_COEFFS_A(ca2x, g2[i2], g0[i0], g1[i1], g23[i], p2->x, p23x[i], -p0123x, g0123, h0123);
        }
        if (max_a2.y > 0) {
            INIT_VRR_COEFFS_A(ca2y, g2[i2], g0[i0], g1[i1], g23[i], p2->y, p23y[i], -p0123y, g0123, h0123);
        }
        if (max_a2.z > 0) {
            INIT_VRR_COEFFS_A(ca2z, g2[i2], g0[i0], g1[i1], g23[i], p2->z, p23z[i], -p0123z, g0123, h0123);
        }
        if (max_a3.x > 0) {
            INIT_VRR_COEFFS_A(ca3x, g3[i3], g1[i1], g0[i0], g23[i], p3->x, p23x[i], -p0123x, g0123, h0123);
        }
        if (max_a3.y > 0) {
            INIT_VRR_COEFFS_A(ca3y, g3[i3], g1[i1], g0[i0], g23[i], p3->y, p23y[i], -p0123y, g0123, h0123);
        }
        if (max_a3.z > 0) {
            INIT_VRR_COEFFS_A(ca3z, g3[i3], g1[i1], g0[i0], g23[i], p3->z, p23z[i], -p0123z, g0123, h0123);
        }
        if (max_d0.x > 0) {
            INIT_VRR_COEFFS_D(cd0x, g0[i0], g2[i2], g3[i3], g01[i], f01[i], p0->x, p01x[i], p0123x, g0123, h0123);
        }
        if (max_d0.y > 0) {
            INIT_VRR_COEFFS_D(cd0y, g0[i0], g2[i2], g3[i3], g01[i], f01[i], p0->y, p01y[i], p0123y, g0123, h0123);
        }
        if (max_d0.z > 0) {
            INIT_VRR_COEFFS_D(cd0z, g0[i0], g2[i2], g3[i3], g01[i], f01[i], p0->z, p01z[i], p0123z, g0123, h0123);
        }
        if (max_d1.x > 0) {
            INIT_VRR_COEFFS_D(cd1x, g1[i1], g3[i3], g2[i2], g01[i], f01[i], p1->x, p01x[i], p0123x, g0123, h0123);
        }
        if (max_d1.y > 0) {
            INIT_VRR_COEFFS_D(cd1y, g1[i1], g3[i3], g2[i2], g01[i], f01[i], p1->y, p01y[i], p0123y, g0123, h0123);
        }
        if (max_d1.z > 0) {
            INIT_VRR_COEFFS_D(cd1z, g1[i1], g3[i3], g2[i2], g01[i], f01[i], p1->z, p01z[i], p0123z, g0123, h0123);
        }
        if (max_d2.x > 0) {
            INIT_VRR_COEFFS_D(cd2x, g2[i2], g0[i0], g1[i1], g23[i], f23[i], p2->x, p23x[i], -p0123x, g0123, h0123);
        }
        if (max_d2.y > 0) {
            INIT_VRR_COEFFS_D(cd2y, g2[i2], g0[i0], g1[i1], g23[i], f23[i], p2->y, p23y[i], -p0123y, g0123, h0123);
        }
        if (max_d2.z > 0) {
            INIT_VRR_COEFFS_D(cd2z, g2[i2], g0[i0], g1[i1], g23[i], f23[i], p2->z, p23z[i], -p0123z, g0123, h0123);
        }
        if (max_d3.x > 0) {
            INIT_VRR_COEFFS_D(cd3x, g3[i3], g1[i1], g0[i0], g23[i], f23[i], p3->x, p23x[i], -p0123x, g0123, h0123);
        }
        if (max_d3.y > 0) {
            INIT_VRR_COEFFS_D(cd3y, g3[i3], g1[i1], g0[i0], g23[i], f23[i], p3->y, p23y[i], -p0123y, g0123, h0123);
        }
        if (max_d3.z > 0) {
            INIT_VRR_COEFFS_D(cd3z, g3[i3], g1[i1], g0[i0], g23[i], f23[i], p3->z, p23z[i], -p0123z, g0123, h0123);
        }
    }
    }
    }
    }
    size_t io = 0;
    for (size_t id = 0; id < nd; id++) {
    for (size_t ia3 = 0; ia3 < na3; ia3++) {
    for (size_t ia2 = 0; ia2 < na2; ia2++) {
    for (size_t ia1 = 0; ia1 < na1; ia1++) {
    for (size_t ia0 = 0; ia0 < na0; ia0++) {
        itg->v.p[io] = 0.0;
        gtoint__stack__reset(&(itg->s), ng0123);
        {
            stack_index_t t = { 0 };
            t.i.eri.a[0] = a0[ia0];
            t.i.eri.a[1] = a1[ia1];
            t.i.eri.a[2] = a2[ia2];
            t.i.eri.a[3] = a3[ia3];
            t.i.eri.d[0] = d0[id];
            t.i.eri.d[1] = d1[id];
            t.i.eri.d[2] = d2[id];
            t.i.eri.d[3] = d3[id];
            t.i.eri.m = 0;
            t.o = STACK_VOID_INDEX;
            t.b = false;
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY;
            size_t it;
            gtoint__stack__write(&(itg->s), &t, &it);
            double *const c = gtoint__stack__coefficients(&(itg->s), it);
            double *const v = gtoint__stack__integrals(&(itg->s), it);
            for (size_t i3 = 0; i3 < ng3; i3++) {
            for (size_t i2 = 0; i2 < ng2; i2++) {
            for (size_t i1 = 0; i1 < ng1; i1++) {
            for (size_t i0 = 0; i0 < ng0; i0++) {
                const size_t i = i0 + ng0 * (i1 + ng1 * (i2 + ng2 * i3));
                c[i] = c0[i0 + ng0 * ia0] * c1[i1 + ng1 * ia1] * c2[i2 + ng2 * ia2] * c3[i3 + ng3 * ia3];
                v[i] = 0.0;
            }
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
                gtoint__cache__reference_to_store_eri(&(itg->c), &s.i, &w);
                for (size_t i = 0; i < ng0123; i++) w[i] = v[i];
                if (s.o == STACK_VOID_INDEX) {
                    double o = 0.0;
                    for (size_t i = 0; i < ng0123; i++) o += c[i] * v[i];
                    itg->v.p[io] += o;
                }
                else {
                    double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                    for (size_t i = 0; i < ng0123; i++) o[i] += c[i] * v[i];
                }
                gtoint__stack__pop(&(itg->s));
            }
            else {
                const double *w;
                if (gtoint__cache__reference_to_fetch_eri(&(itg->c), &s.i, &w)) {
                    const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                    if (s.o == STACK_VOID_INDEX) {
                        double o = 0.0;
                        for (size_t i = 0; i < ng0123; i++) o += c[i] * w[i];
                        itg->v.p[io] += o;
                    }
                    else {
                        double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                        for (size_t i = 0; i < ng0123; i++) o[i] += c[i] * w[i];
                    }
                    gtoint__stack__pop(&(itg->s));
                }
                else {
                    {
                        s.b = true;
                        gtoint__stack__write(&(itg->s), &s, &is);
                        double *const v = gtoint__stack__integrals(&(itg->s), is);
                        for (size_t i = 0; i < ng0123; i++) v[i] = 0.0;
                    }
                    if (s.i.eri.a[0].x > 0) {
                        s.i.eri.a[0].x--;
                        EXPAND_VRR(ca0x, x, 0, 1, 2, 3);
                    }
                    else if (s.i.eri.a[0].y > 0) {
                        s.i.eri.a[0].y--;
                        EXPAND_VRR(ca0y, y, 0, 1, 2, 3);
                    }
                    else if (s.i.eri.a[0].z > 0) {
                        s.i.eri.a[0].z--;
                        EXPAND_VRR(ca0z, z, 0, 1, 2, 3);
                    }
                    else if (s.i.eri.a[1].x > 0) {
                        s.i.eri.a[1].x--;
                        EXPAND_VRR(ca1x, x, 1, 0, 3, 2);
                    }
                    else if (s.i.eri.a[1].y > 0) {
                        s.i.eri.a[1].y--;
                        EXPAND_VRR(ca1y, y, 1, 0, 3, 2);
                    }
                    else if (s.i.eri.a[1].z > 0) {
                        s.i.eri.a[1].z--;
                        EXPAND_VRR(ca1z, z, 1, 0, 3, 2);
                    }
                    else if (s.i.eri.a[2].x > 0) {
                        s.i.eri.a[2].x--;
                        EXPAND_VRR(ca2x, x, 2, 3, 0, 1);
                    }
                    else if (s.i.eri.a[2].y > 0) {
                        s.i.eri.a[2].y--;
                        EXPAND_VRR(ca2y, y, 2, 3, 0, 1);
                    }
                    else if (s.i.eri.a[2].z > 0) {
                        s.i.eri.a[2].z--;
                        EXPAND_VRR(ca2z, z, 2, 3, 0, 1);
                    }
                    else if (s.i.eri.a[3].x > 0) {
                        s.i.eri.a[3].x--;
                        EXPAND_VRR(ca3x, x, 3, 2, 1, 0);
                    }
                    else if (s.i.eri.a[3].y > 0) {
                        s.i.eri.a[3].y--;
                        EXPAND_VRR(ca3y, y, 3, 2, 1, 0);
                    }
                    else if (s.i.eri.a[3].z > 0) {
                        s.i.eri.a[3].z--;
                        EXPAND_VRR(ca3z, z, 3, 2, 1, 0);
                    }
                    else if (s.i.eri.d[0].x > 0) {
                        s.i.eri.d[0].x--;
                        EXPAND_VRR(cd0x, x, 0, 1, 2, 3);
                    }
                    else if (s.i.eri.d[0].y > 0) {
                        s.i.eri.d[0].y--;
                        EXPAND_VRR(cd0y, y, 0, 1, 2, 3);
                    }
                    else if (s.i.eri.d[0].z > 0) {
                        s.i.eri.d[0].z--;
                        EXPAND_VRR(cd0z, z, 0, 1, 2, 3);
                    }
                    else if (s.i.eri.d[1].x > 0) {
                        s.i.eri.d[1].x--;
                        EXPAND_VRR(cd1x, x, 1, 0, 3, 2);
                    }
                    else if (s.i.eri.d[1].y > 0) {
                        s.i.eri.d[1].y--;
                        EXPAND_VRR(cd1y, y, 1, 0, 3, 2);
                    }
                    else if (s.i.eri.d[1].z > 0) {
                        s.i.eri.d[1].z--;
                        EXPAND_VRR(cd1z, z, 1, 0, 3, 2);
                    }
                    else if (s.i.eri.d[2].x > 0) {
                        s.i.eri.d[2].x--;
                        EXPAND_VRR(cd2x, x, 2, 3, 0, 1);
                    }
                    else if (s.i.eri.d[2].y > 0) {
                        s.i.eri.d[2].y--;
                        EXPAND_VRR(cd2y, y, 2, 3, 0, 1);
                    }
                    else if (s.i.eri.d[2].z > 0) {
                        s.i.eri.d[2].z--;
                        EXPAND_VRR(cd2z, z, 2, 3, 0, 1);
                    }
                    else if (s.i.eri.d[3].x > 0) {
                        s.i.eri.d[3].x--;
                        EXPAND_VRR(cd3x, x, 3, 2, 1, 0);
                    }
                    else if (s.i.eri.d[3].y > 0) {
                        s.i.eri.d[3].y--;
                        EXPAND_VRR(cd3y, y, 3, 2, 1, 0);
                    }
                    else if (s.i.eri.d[3].z > 0) {
                        s.i.eri.d[3].z--;
                        EXPAND_VRR(cd3z, z, 3, 2, 1, 0);
                    }
                    else {
                        double *v;
                        gtoint__cache__reference_to_store_eri(&(itg->c), &s.i, &v);
                        for (size_t i = 0; i < ng0123; i++) {
                            v[i] = o0123[i] * gtoint__compute_boys_function(s.i.eri.m, q0123[i], itg->tol);
                        }
                        const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                        if (s.o == STACK_VOID_INDEX) {
                            double o = 0.0;
                            for (size_t i = 0; i < ng0123; i++) o += c[i] * v[i];
                            itg->v.p[io] += o;
                        }
                        else {
                            double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                            for (size_t i = 0; i < ng0123; i++) o[i] += c[i] * v[i];
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
    }
    }
    return GTOINT_ERROR_OK;
#undef NVAR
#undef INIT_VRR_COEFFS_A
#undef INIT_VRR_COEFFS_D
#undef INIT_VRR_COEFFS_C
#undef EXPAND_VRR
#undef EXPAND_VRR_C
}

gtoint_error_t gtoint_compute_electron_repulsion_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    const gtoint_double3_t *p2, gtoint_basis_shell_t bas2,
    const gtoint_double3_t *p3, gtoint_basis_shell_t bas3,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *d2, const gtoint_int3_t *d3,
    double *out
) {
    if (nd < 0) return GTOINT_ERROR_ARGUMENT;
    const gtoint_error_t e = compute_electron_repulsion_integrals_(
        itg,
        p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p,
        p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p,
        p2, bas2->spec.mom.n, bas2->spec.mom.p, bas2->prim.n, bas2->prim.expo.p, bas2->prim.coef.p,
        p3, bas3->spec.mom.n, bas3->spec.mom.p, bas3->prim.n, bas3->prim.expo.p, bas3->prim.coef.p,
        nd, d0, d1, d2, d3
    );
    if (e != GTOINT_ERROR_OK) return e;
    const size_t mm0 = bas0->spec.i.m.p[bas0->spec.n];
    const size_t mm1 = bas1->spec.i.m.p[bas1->spec.n];
    const size_t mm2 = bas2->spec.i.m.p[bas2->spec.n];
    const size_t mm3 = bas3->spec.i.m.p[bas3->spec.n];
    const size_t mb0 = bas0->spec.i.b.p[bas0->spec.n];
    const size_t mb1 = bas1->spec.i.b.p[bas1->spec.n];
    const size_t mb2 = bas2->spec.i.b.p[bas2->spec.n];
    const size_t mb3 = bas3->spec.i.b.p[bas3->spec.n];
    for (int id = 0; id < nd; id++) {
        for (size_t ia3 = 0; ia3 < bas3->spec.n; ia3++) {
            const size_t im3 = bas3->spec.i.m.p[ia3];
            const size_t ib3 = bas3->spec.i.b.p[ia3];
            const spherical_harmonics_t *const s3 = &(bas3->spec.sph.p[ia3]);
        for (size_t ia2 = 0; ia2 < bas2->spec.n; ia2++) {
            const size_t im2 = bas2->spec.i.m.p[ia2];
            const size_t ib2 = bas2->spec.i.b.p[ia2];
            const spherical_harmonics_t *const s2 = &(bas2->spec.sph.p[ia2]);
        for (size_t ia1 = 0; ia1 < bas1->spec.n; ia1++) {
            const size_t im1 = bas1->spec.i.m.p[ia1];
            const size_t ib1 = bas1->spec.i.b.p[ia1];
            const spherical_harmonics_t *const s1 = &(bas1->spec.sph.p[ia1]);
        for (size_t ia0 = 0; ia0 < bas0->spec.n; ia0++) {
            const size_t im0 = bas0->spec.i.m.p[ia0];
            const size_t ib0 = bas0->spec.i.b.p[ia0];
            const spherical_harmonics_t *const s0 = &(bas0->spec.sph.p[ia0]);
            for (size_t jb3 = 0; jb3 < s3->n; jb3++) {
            for (size_t jb2 = 0; jb2 < s2->n; jb2++) {
            for (size_t jb1 = 0; jb1 < s1->n; jb1++) {
            for (size_t jb0 = 0; jb0 < s0->n; jb0++) {
                double v3 = 0.0;
                for (size_t jm3 = 0; jm3 < s3->l.p[jb3]; jm3++) {
                    const size_t k3 = jm3 + s3->m * jb3;
                    double v2 = 0.0;
                    for (size_t jm2 = 0; jm2 < s2->l.p[jb2]; jm2++) {
                        const size_t k2 = jm2 + s2->m * jb2;
                        double v1 = 0.0;
                        for (size_t jm1 = 0; jm1 < s1->l.p[jb1]; jm1++) {
                            const size_t k1 = jm1 + s1->m * jb1;
                            double v0 = 0.0;
                            for (size_t jm0 = 0; jm0 < s0->l.p[jb0]; jm0++) {
                                const size_t k0 = jm0 + s0->m * jb0;
                                v0 += s0->c.p[k0] * itg->v.p[
                                    (im0 + s0->i.p[k0]) + mm0 * ((im1 + s1->i.p[k1]) + mm1 * ((im2 + s2->i.p[k2]) + mm2 * ((im3 + s3->i.p[k3]) + mm3 * id)))
                                ];
                            }
                            v1 += s1->c.p[k1] * v0;
                        }
                        v2 += s2->c.p[k2] * v1;
                    }
                    v3 += s3->c.p[k3] * v2;
                }
                out[(ib0 + jb0) + mb0 * ((ib1 + jb1) + mb1 * ((ib2 + jb2) + mb2 * ((ib3 + jb3) + mb3 * id)))] = v3;
            }
            }
            }
            }
        }
        }
        }
        }
    }
    return GTOINT_ERROR_OK;
}
