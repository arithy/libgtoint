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

#include "gtoint-private.h"

#include <math.h>
#include <assert.h>

#define D_RPI3 5.5683279968317078452848179821188 /* PI^(3/2) */

gtoint_error_t gtoint__compute_scalar_ecp_type1_integrals_2(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    const double3_t *pc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
) {
#define NTERM 6
#define NPTR (3 * 5)
#define NVAR (1 + NTERM * NPTR)
#define INIT_VRR_COEFFS_A(out, xyz, g0_, g1_, g2_, g012_, p0_, p1_, p2_) /* reference variable: i */ \
    { \
        out[0][i] = ((g1_) * ((p1_)->xyz - (p0_)->xyz) + (g2_) * ((p2_)->xyz - (p0_)->xyz)) * (g012_); \
        out[1][i] = -((g1_) + (g2_)) * (g012_); \
        out[2][i] = 0.5 * (g012_); \
        out[3][i] = (g1_)* (g012_); \
        out[4][i] = 0.5 * (g012_); \
        out[5][i] = (g2_)* (g012_); \
    }
#define INIT_VRR_COEFFS_D(out, xyz, g0_, g1_, g2_, g012_, p0_, p1_, p2_) /* reference variable: i */ \
    { \
        out[0][i] = 2.0 * (g0_) * ((g1_) * ((p1_)->xyz - (p0_)->xyz) + (g2_) * ((p2_)->xyz - (p0_)->xyz)) * (g012_); \
        out[1][i] = -2.0 * (g0_) * ((g1_) + (g2_)) * (g012_); \
        out[2][i] = -((g1_) + (g2_)) * (g012_); \
        out[3][i] = 2.0 * (g0_) * (g1_) * (g012_); \
        out[4][i] = (g0_) * (g012_); \
        out[5][i] = 2.0 * (g0_) * (g2_) * (g012_); \
    }
#define INIT_VRR_COEFFS_C(out, xyz, g0_, g1_, g2_, g012_, p0_, p1_, p2_) /* reference variable: i */ \
    { \
        out[0][i] = 2.0 * (g2_) * ((g1_) * ((p1_)->xyz - (p2_)->xyz) + (g0_) * ((p0_)->xyz - (p2_)->xyz)) * (g012_); \
        out[1][i] = 2.0 * (g0_) * (g2_) * (g012_); \
        out[2][i] = (g2_) * (g012_); \
        out[3][i] = 2.0 * (g1_) * (g2_) * (g012_); \
        out[4][i] = (g2_) * (g012_); \
        out[5][i] = -2.0 * (g2_) * ((g0_) + (g1_)) * (g012_); \
    }
#define EXPAND_VRR(coe, xyz, i0, i1) /* reference variables: itg, e, s, is, ng01c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[0][i]; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.d[i0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.d[i0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[1][i] * s.i.ecp0.d[i0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.a[i0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.a[i0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[2][i] * s.i.ecp0.a[i0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.d[i1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.d[i1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[3][i] * s.i.ecp0.d[i1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.a[i1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.a[i1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[4][i] * s.i.ecp0.a[i1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[5][i] * s.i.ecp0.d[2].xyz; \
                v[i] = 0.0; \
            } \
        } \
    }
    if (!gtoint__double_array__resize(&(itg->v), na0 * na1 * nd)) return GTOINT_ERROR_MEMORY;
    const double r01x = p1->x - p0->x;
    const double r01y = p1->y - p0->y;
    const double r01z = p1->z - p0->z;
    const double r01w = r01x * r01x + r01y * r01y + r01z * r01z;
    const double rc0x = p0->x - pc->x;
    const double rc0y = p0->y - pc->y;
    const double rc0z = p0->z - pc->z;
    const double rc0w = rc0x * rc0x + rc0y * rc0y + rc0z * rc0z;
    const double rc1x = p1->x - pc->x;
    const double rc1y = p1->y - pc->y;
    const double rc1z = p1->z - pc->z;
    const double rc1w = rc1x * rc1x + rc1y * rc1y + rc1z * rc1z;
    int3_t max_a0 = { 0, 0, 0 };
    int3_t max_a1 = { 0, 0, 0 };
    int3_t max_d0 = { 0, 0, 0 };
    int3_t max_d1 = { 0, 0, 0 };
    int3_t max_dc = { 0, 0, 0 };
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
    for (size_t i = 0; i < nd; i++) {
        if (max_d0.x < d0[i].x) max_d0.x = d0[i].x;
        if (max_d0.y < d0[i].y) max_d0.y = d0[i].y;
        if (max_d0.z < d0[i].z) max_d0.z = d0[i].z;
        if (max_d1.x < d1[i].x) max_d1.x = d1[i].x;
        if (max_d1.y < d1[i].y) max_d1.y = d1[i].y;
        if (max_d1.z < d1[i].z) max_d1.z = d1[i].z;
        if (max_dc.x < dc[i].x) max_dc.x = dc[i].x;
        if (max_dc.y < dc[i].y) max_dc.y = dc[i].y;
        if (max_dc.z < dc[i].z) max_dc.z = dc[i].z;
    }
    const size_t ng01c = ng0 * ng1 * ngc;
    if (!gtoint__cache__reset(&(itg->c), ng01c)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_array__resize(&(itg->w), ng01c * NVAR)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_pointer_array__resize(&(itg->p), NTERM * NPTR)) return GTOINT_ERROR_MEMORY;
    size_t ivar = 0;
    double *const o01c  = itg->w.p + ng01c * ivar++;
    size_t iptr = 0;
    double **const ca0x = itg->p.p + NTERM * iptr++;
    double **const ca0y = itg->p.p + NTERM * iptr++;
    double **const ca0z = itg->p.p + NTERM * iptr++;
    double **const ca1x = itg->p.p + NTERM * iptr++;
    double **const ca1y = itg->p.p + NTERM * iptr++;
    double **const ca1z = itg->p.p + NTERM * iptr++;
    double **const cd0x = itg->p.p + NTERM * iptr++;
    double **const cd0y = itg->p.p + NTERM * iptr++;
    double **const cd0z = itg->p.p + NTERM * iptr++;
    double **const cd1x = itg->p.p + NTERM * iptr++;
    double **const cd1y = itg->p.p + NTERM * iptr++;
    double **const cd1z = itg->p.p + NTERM * iptr++;
    double **const cdcx = itg->p.p + NTERM * iptr++;
    double **const cdcy = itg->p.p + NTERM * iptr++;
    double **const cdcz = itg->p.p + NTERM * iptr++;
    for (size_t i = 0; i < NTERM; i++) ca0x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca0y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca0z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca1x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca1y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca1z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd0x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd0y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd0z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd1x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd1y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd1z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cdcx[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cdcy[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cdcz[i] = itg->w.p + ng01c * ivar++;
    assert(ivar == NVAR);
    assert(iptr == NPTR);
    for (size_t ic = 0; ic < ngc; ic++) {
    for (size_t i1 = 0; i1 < ng1; i1++) {
    for (size_t i0 = 0; i0 < ng0; i0++) {
        const size_t i = i0 + ng0 * (i1 + ng1 * ic);
        const double g01c = 1.0 / (g0[i0] + g1[i1] + gc[ic]);
        o01c[i] = D_RPI3 * g01c * sqrt(g01c) * exp(-(
            g0[i0] * g1[i1] * r01w + g1[i1] * gc[ic] * rc1w + gc[ic] * g0[i0] * rc0w
        ) * g01c);
        if (max_a0.x > 0) {
            INIT_VRR_COEFFS_A(ca0x, x, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_a0.y > 0) {
            INIT_VRR_COEFFS_A(ca0y, y, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_a0.z > 0) {
            INIT_VRR_COEFFS_A(ca0z, z, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_a1.x > 0) {
            INIT_VRR_COEFFS_A(ca1x, x, g1[i1], g0[i0], gc[ic], g01c, p1, p0, pc);
        }
        if (max_a1.y > 0) {
            INIT_VRR_COEFFS_A(ca1y, y, g1[i1], g0[i0], gc[ic], g01c, p1, p0, pc);
        }
        if (max_a1.z > 0) {
            INIT_VRR_COEFFS_A(ca1z, z, g1[i1], g0[i0], gc[ic], g01c, p1, p0, pc);
        }
        if (max_d0.x > 0) {
            INIT_VRR_COEFFS_D(cd0x, x, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_d0.y > 0) {
            INIT_VRR_COEFFS_D(cd0y, y, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_d0.z > 0) {
            INIT_VRR_COEFFS_D(cd0z, z, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_d1.x > 0) {
            INIT_VRR_COEFFS_D(cd1x, x, g1[i1], g0[i0], gc[ic], g01c, p1, p0, pc);
        }
        if (max_d1.y > 0) {
            INIT_VRR_COEFFS_D(cd1y, y, g1[i1], g0[i0], gc[ic], g01c, p1, p0, pc);
        }
        if (max_d1.z > 0) {
            INIT_VRR_COEFFS_D(cd1z, z, g1[i1], g0[i0], gc[ic], g01c, p1, p0, pc);
        }
        if (max_dc.x > 0) {
            INIT_VRR_COEFFS_C(cdcx, x, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_dc.y > 0) {
            INIT_VRR_COEFFS_C(cdcy, y, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
        if (max_dc.z > 0) {
            INIT_VRR_COEFFS_C(cdcz, z, g0[i0], g1[i1], gc[ic], g01c, p0, p1, pc);
        }
    }
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
            t.i.ecp0.a[0] = a0[ia0];
            t.i.ecp0.a[1] = a1[ia1];
            t.i.ecp0.d[0] = d0[id];
            t.i.ecp0.d[1] = d1[id];
            t.i.ecp0.d[2] = dc[id];
            t.o = STACK_VOID_INDEX;
            t.b = false;
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY;
            size_t it;
            gtoint__stack__write(&(itg->s), &t, &it);
            double *const c = gtoint__stack__coefficients(&(itg->s), it);
            double *const v = gtoint__stack__integrals(&(itg->s), it);
            for (size_t ic = 0; ic < ngc; ic++) {
            for (size_t i1 = 0; i1 < ng1; i1++) {
            for (size_t i0 = 0; i0 < ng0; i0++) {
                const size_t i = i0 + ng0 * (i1 + ng1 * ic);
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
                if (!gtoint__cache__reference_to_store_ecp0(&(itg->c), &s.i, &w)) return GTOINT_ERROR_MEMORY;
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
                if (gtoint__cache__reference_to_fetch_ecp0(&(itg->c), &s.i, &w)) {
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
                    if (s.i.ecp0.a[0].x > 0) {
                        s.i.ecp0.a[0].x--;
                        EXPAND_VRR(ca0x, x, 0, 1);
                    }
                    else if (s.i.ecp0.a[0].y > 0) {
                        s.i.ecp0.a[0].y--;
                        EXPAND_VRR(ca0y, y, 0, 1);
                    }
                    else if (s.i.ecp0.a[0].z > 0) {
                        s.i.ecp0.a[0].z--;
                        EXPAND_VRR(ca0z, z, 0, 1);
                    }
                    else if (s.i.ecp0.a[1].x > 0) {
                        s.i.ecp0.a[1].x--;
                        EXPAND_VRR(ca1x, x, 1, 0);
                    }
                    else if (s.i.ecp0.a[1].y > 0) {
                        s.i.ecp0.a[1].y--;
                        EXPAND_VRR(ca1y, y, 1, 0);
                    }
                    else if (s.i.ecp0.a[1].z > 0) {
                        s.i.ecp0.a[1].z--;
                        EXPAND_VRR(ca1z, z, 1, 0);
                    }
                    else if (s.i.ecp0.d[0].x > 0) {
                        s.i.ecp0.d[0].x--;
                        EXPAND_VRR(cd0x, x, 0, 1);
                    }
                    else if (s.i.ecp0.d[0].y > 0) {
                        s.i.ecp0.d[0].y--;
                        EXPAND_VRR(cd0y, y, 0, 1);
                    }
                    else if (s.i.ecp0.d[0].z > 0) {
                        s.i.ecp0.d[0].z--;
                        EXPAND_VRR(cd0z, z, 0, 1);
                    }
                    else if (s.i.ecp0.d[1].x > 0) {
                        s.i.ecp0.d[1].x--;
                        EXPAND_VRR(cd1x, x, 1, 0);
                    }
                    else if (s.i.ecp0.d[1].y > 0) {
                        s.i.ecp0.d[1].y--;
                        EXPAND_VRR(cd1y, y, 1, 0);
                    }
                    else if (s.i.ecp0.d[1].z > 0) {
                        s.i.ecp0.d[1].z--;
                        EXPAND_VRR(cd1z, z, 1, 0);
                    }
                    else if (s.i.ecp0.d[2].x > 0) {
                        s.i.ecp0.d[2].x--;
                        EXPAND_VRR(cdcx, x, 0, 1);
                    }
                    else if (s.i.ecp0.d[2].y > 0) {
                        s.i.ecp0.d[2].y--;
                        EXPAND_VRR(cdcy, y, 0, 1);
                    }
                    else if (s.i.ecp0.d[2].z > 0) {
                        s.i.ecp0.d[2].z--;
                        EXPAND_VRR(cdcz, z, 0, 1);
                    }
                    else {
                        double *v;
                        if (!gtoint__cache__reference_to_store_ecp0(&(itg->c), &s.i, &v)) return GTOINT_ERROR_MEMORY;
                        for (size_t i = 0; i < ng01c; i++) {
                            v[i] = o01c[i];
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
#undef NTERM
#undef NPTR
#undef NVAR
#undef INIT_VRR_COEFFS_A
#undef INIT_VRR_COEFFS_D
#undef INIT_VRR_COEFFS_C
#undef EXPAND_VRR
}

gtoint_error_t gtoint__compute_weighted_scalar_ecp_type1_integrals_2(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0, const double3_t *pw0, double gw0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1, const double3_t *pw1, double gw1,
    const double3_t *pc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
) {
#define NTERM 6
#define NPTR (3 * 5)
#define NVAR (1 + NTERM * NPTR)
#define INIT_VRR_COEFFS_A(out, g0_, g1_, gc_, r01_, r0c_, r0w_, r1c_, rwc_, g01wc_) /* reference variable: i, gw */ \
    { \
        out[0][i] = ((g1_) * (r01_) + (gc_) * (r0c_) + gw * (r0w_)) * (g01wc_); \
        out[1][i] = -((g1_) + (gc_) + gw) * (g01wc_); \
        out[2][i] = 0.5 * (g01wc_); \
        out[3][i] = (g1_) * (g01wc_); \
        out[4][i] = 0.5 * (g01wc_); \
        out[5][i] = (gc_) * (g01wc_); \
    }
#define INIT_VRR_COEFFS_D(out, g0_, g1_, gc_, r01_, r0c_, r0w_, r1c_, rwc_, g01wc_) /* reference variable: i, gw */ \
    { \
        out[0][i] = 2.0 * (g0_) * ((g1_) * (r01_) + (gc_) * (r0c_) + gw * (r0w_)) * (g01wc_); \
        out[1][i] = -2.0 * (g0_) * ((g1_) + (gc_) + gw) * (g01wc_); \
        out[2][i] = -((g1_) + (gc_) + gw) * (g01wc_); \
        out[3][i] = 2.0 * (g0_) * (g1_) * (g01wc_); \
        out[4][i] = (g0_) * (g01wc_); \
        out[5][i] = 2.0 * (g0_) * (gc_) * (g01wc_); \
    }
#define INIT_VRR_COEFFS_C(out, g0_, g1_, gc_, r01_, r0c_, r0w_, r1c_, rwc_, g01wc_) /* reference variable: i, gw */ \
    { \
        out[0][i] = -2.0 * (gc_) * ((g0_) * (r0c_) + (g1_) * (r1c_) + gw * (rwc_)) * (g01wc_); \
        out[1][i] = 2.0 * (g0_) * (gc_) * (g01wc_); \
        out[2][i] = (gc_) * (g01wc_); \
        out[3][i] = 2.0 * (g1_) * (gc_) * (g01wc_); \
        out[4][i] = (gc_) * (g01wc_); \
        out[5][i] = -2.0 * ((g0_) + (g1_) + gw) * (gc_) * (g01wc_); \
    }
#define EXPAND_VRR(coe, xyz, i0, i1) /* reference variables: itg, e, s, is, ng01c */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[0][i]; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.d[i0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.d[i0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[1][i] * s.i.ecp0.d[i0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.a[i0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.a[i0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[2][i] * s.i.ecp0.a[i0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.d[i1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.d[i1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[3][i] * s.i.ecp0.d[i1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.a[i1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.a[i1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[4][i] * s.i.ecp0.a[i1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.ecp0.d[2].xyz > 0) { \
            stack_index_t t = s; \
            t.i.ecp0.d[2].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01c; i++) { \
                c[i] = coe[5][i] * s.i.ecp0.d[2].xyz; \
                v[i] = 0.0; \
            } \
        } \
    }
    if (!gtoint__double_array__resize(&(itg->v), na0 * na1 * nd)) return GTOINT_ERROR_MEMORY;
    const double gw = gw0 + gw1;
    const double pwx = (gw0 * pw0->x + gw1 * pw1->x) / gw;
    const double pwy = (gw0 * pw0->y + gw1 * pw1->y) / gw;
    const double pwz = (gw0 * pw0->z + gw1 * pw1->z) / gw;
    const double rwwx = pw1->x - pw0->x;
    const double rwwy = pw1->y - pw0->y;
    const double rwwz = pw1->z - pw0->z;
    const double rwww = rwwx * rwwx + rwwy * rwwy + rwwz * rwwz;
    const double cw = exp(-gw0 * gw1 * rwww / gw);
    const double r01x = p1->x - p0->x;
    const double r01y = p1->y - p0->y;
    const double r01z = p1->z - p0->z;
    const double r01w = r01x * r01x + r01y * r01y + r01z * r01z;
    const double r0cx = pc->x - p0->x;
    const double r0cy = pc->y - p0->y;
    const double r0cz = pc->z - p0->z;
    const double r0cw = r0cx * r0cx + r0cy * r0cy + r0cz * r0cz;
    const double r1cx = pc->x - p1->x;
    const double r1cy = pc->y - p1->y;
    const double r1cz = pc->z - p1->z;
    const double r1cw = r1cx * r1cx + r1cy * r1cy + r1cz * r1cz;
    const double r0wx = pwx - p0->x;
    const double r0wy = pwy - p0->y;
    const double r0wz = pwz - p0->z;
    const double r0ww = r0wx * r0wx + r0wy * r0wy + r0wz * r0wz;
    const double r1wx = pwx - p1->x;
    const double r1wy = pwy - p1->y;
    const double r1wz = pwz - p1->z;
    const double r1ww = r1wx * r1wx + r1wy * r1wy + r1wz * r1wz;
    const double rwcx = pc->x - pwx;
    const double rwcy = pc->y - pwy;
    const double rwcz = pc->z - pwz;
    const double rwcw = rwcx * rwcx + rwcy * rwcy + rwcz * rwcz;
    int3_t max_a0 = { 0, 0, 0 };
    int3_t max_a1 = { 0, 0, 0 };
    int3_t max_d0 = { 0, 0, 0 };
    int3_t max_d1 = { 0, 0, 0 };
    int3_t max_dc = { 0, 0, 0 };
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
    for (size_t i = 0; i < nd; i++) {
        if (max_d0.x < d0[i].x) max_d0.x = d0[i].x;
        if (max_d0.y < d0[i].y) max_d0.y = d0[i].y;
        if (max_d0.z < d0[i].z) max_d0.z = d0[i].z;
        if (max_d1.x < d1[i].x) max_d1.x = d1[i].x;
        if (max_d1.y < d1[i].y) max_d1.y = d1[i].y;
        if (max_d1.z < d1[i].z) max_d1.z = d1[i].z;
        if (max_dc.x < dc[i].x) max_dc.x = dc[i].x;
        if (max_dc.y < dc[i].y) max_dc.y = dc[i].y;
        if (max_dc.z < dc[i].z) max_dc.z = dc[i].z;
    }
    const size_t ng01c = ng0 * ng1 * ngc;
    if (!gtoint__cache__reset(&(itg->c), ng01c)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_array__resize(&(itg->w), ng01c * NVAR)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_pointer_array__resize(&(itg->p), NTERM * NPTR)) return GTOINT_ERROR_MEMORY;
    size_t ivar = 0;
    double *const o01c  = itg->w.p + ng01c * ivar++;
    size_t iptr = 0;
    double **const ca0x = itg->p.p + NTERM * iptr++;
    double **const ca0y = itg->p.p + NTERM * iptr++;
    double **const ca0z = itg->p.p + NTERM * iptr++;
    double **const ca1x = itg->p.p + NTERM * iptr++;
    double **const ca1y = itg->p.p + NTERM * iptr++;
    double **const ca1z = itg->p.p + NTERM * iptr++;
    double **const cd0x = itg->p.p + NTERM * iptr++;
    double **const cd0y = itg->p.p + NTERM * iptr++;
    double **const cd0z = itg->p.p + NTERM * iptr++;
    double **const cd1x = itg->p.p + NTERM * iptr++;
    double **const cd1y = itg->p.p + NTERM * iptr++;
    double **const cd1z = itg->p.p + NTERM * iptr++;
    double **const cdcx = itg->p.p + NTERM * iptr++;
    double **const cdcy = itg->p.p + NTERM * iptr++;
    double **const cdcz = itg->p.p + NTERM * iptr++;
    for (size_t i = 0; i < NTERM; i++) ca0x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca0y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca0z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca1x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca1y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) ca1z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd0x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd0y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd0z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd1x[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd1y[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cd1z[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cdcx[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cdcy[i] = itg->w.p + ng01c * ivar++;
    for (size_t i = 0; i < NTERM; i++) cdcz[i] = itg->w.p + ng01c * ivar++;
    assert(ivar == NVAR);
    assert(iptr == NPTR);
    for (size_t ic = 0; ic < ngc; ic++) {
    for (size_t i1 = 0; i1 < ng1; i1++) {
    for (size_t i0 = 0; i0 < ng0; i0++) {
        const size_t i = i0 + ng0 * (i1 + ng1 * ic);
        const double g01wc = 1.0 / (g0[i0] + g1[i1] + gc[ic] + gw);
        o01c[i] = D_RPI3 * g01wc * sqrt(g01wc) * exp(-(
            g0[i0] * g1[i1] * r01w + g1[i1] * gc[ic] * r1cw + gc[ic] * g0[i0] * r0cw +
            gw * (g0[i0] * r0ww + g1[i1] * r1ww + gc[ic] * rwcw)
        ) * g01wc) * cw;
        if (max_a0.x > 0) {
            INIT_VRR_COEFFS_A(ca0x, g0[i0], g1[i1], gc[ic], r01x, r0cx, r0wx, r1cx, rwcx, g01wc);
        }
        if (max_a0.y > 0) {
            INIT_VRR_COEFFS_A(ca0y, g0[i0], g1[i1], gc[ic], r01y, r0cy, r0wy, r1cy, rwcy, g01wc);
        }
        if (max_a0.z > 0) {
            INIT_VRR_COEFFS_A(ca0z, g0[i0], g1[i1], gc[ic], r01z, r0cz, r0wz, r1cz, rwcz, g01wc);
        }
        if (max_a1.x > 0) {
            INIT_VRR_COEFFS_A(ca1x, g1[i1], g0[i0], gc[ic], -r01x, r1cx, r1wx, r0cx, rwcx, g01wc);
        }
        if (max_a1.y > 0) {
            INIT_VRR_COEFFS_A(ca1y, g1[i1], g0[i0], gc[ic], -r01y, r1cy, r1wy, r0cy, rwcy, g01wc);
        }
        if (max_a1.z > 0) {
            INIT_VRR_COEFFS_A(ca1z, g1[i1], g0[i0], gc[ic], -r01z, r1cz, r1wz, r0cz, rwcz, g01wc);
        }
        if (max_d0.x > 0) {
            INIT_VRR_COEFFS_D(cd0x, g0[i0], g1[i1], gc[ic], r01x, r0cx, r0wx, r1cx, rwcx, g01wc);
        }
        if (max_d0.y > 0) {
            INIT_VRR_COEFFS_D(cd0y, g0[i0], g1[i1], gc[ic], r01y, r0cy, r0wy, r1cy, rwcy, g01wc);
        }
        if (max_d0.z > 0) {
            INIT_VRR_COEFFS_D(cd0z, g0[i0], g1[i1], gc[ic], r01z, r0cz, r0wz, r1cz, rwcz, g01wc);
        }
        if (max_d1.x > 0) {
            INIT_VRR_COEFFS_D(cd1x, g1[i1], g0[i0], gc[ic], -r01x, r1cx, r1wx, r0cx, rwcx, g01wc);
        }
        if (max_d1.y > 0) {
            INIT_VRR_COEFFS_D(cd1y, g1[i1], g0[i0], gc[ic], -r01y, r1cy, r1wy, r0cy, rwcy, g01wc);
        }
        if (max_d1.z > 0) {
            INIT_VRR_COEFFS_D(cd1z, g1[i1], g0[i0], gc[ic], -r01z, r1cz, r1wz, r0cz, rwcz, g01wc);
        }
        if (max_dc.x > 0) {
            INIT_VRR_COEFFS_C(cdcx, g0[i0], g1[i1], gc[ic], r01x, r0cx, r0wx, r1cx, rwcx, g01wc);
        }
        if (max_dc.y > 0) {
            INIT_VRR_COEFFS_C(cdcy, g0[i0], g1[i1], gc[ic], r01y, r0cy, r0wy, r1cy, rwcy, g01wc);
        }
        if (max_dc.z > 0) {
            INIT_VRR_COEFFS_C(cdcz, g0[i0], g1[i1], gc[ic], r01z, r0cz, r0wz, r1cz, rwcz, g01wc);
        }
    }
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
            t.i.ecp0.a[0] = a0[ia0];
            t.i.ecp0.a[1] = a1[ia1];
            t.i.ecp0.d[0] = d0[id];
            t.i.ecp0.d[1] = d1[id];
            t.i.ecp0.d[2] = dc[id];
            t.o = STACK_VOID_INDEX;
            t.b = false;
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY;
            size_t it;
            gtoint__stack__write(&(itg->s), &t, &it);
            double *const c = gtoint__stack__coefficients(&(itg->s), it);
            double *const v = gtoint__stack__integrals(&(itg->s), it);
            for (size_t ic = 0; ic < ngc; ic++) {
            for (size_t i1 = 0; i1 < ng1; i1++) {
            for (size_t i0 = 0; i0 < ng0; i0++) {
                const size_t i = i0 + ng0 * (i1 + ng1 * ic);
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
                if (!gtoint__cache__reference_to_store_ecp0(&(itg->c), &s.i, &w)) return GTOINT_ERROR_MEMORY;
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
                if (gtoint__cache__reference_to_fetch_ecp0(&(itg->c), &s.i, &w)) {
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
                    if (s.i.ecp0.a[0].x > 0) {
                        s.i.ecp0.a[0].x--;
                        EXPAND_VRR(ca0x, x, 0, 1);
                    }
                    else if (s.i.ecp0.a[0].y > 0) {
                        s.i.ecp0.a[0].y--;
                        EXPAND_VRR(ca0y, y, 0, 1);
                    }
                    else if (s.i.ecp0.a[0].z > 0) {
                        s.i.ecp0.a[0].z--;
                        EXPAND_VRR(ca0z, z, 0, 1);
                    }
                    else if (s.i.ecp0.a[1].x > 0) {
                        s.i.ecp0.a[1].x--;
                        EXPAND_VRR(ca1x, x, 1, 0);
                    }
                    else if (s.i.ecp0.a[1].y > 0) {
                        s.i.ecp0.a[1].y--;
                        EXPAND_VRR(ca1y, y, 1, 0);
                    }
                    else if (s.i.ecp0.a[1].z > 0) {
                        s.i.ecp0.a[1].z--;
                        EXPAND_VRR(ca1z, z, 1, 0);
                    }
                    else if (s.i.ecp0.d[0].x > 0) {
                        s.i.ecp0.d[0].x--;
                        EXPAND_VRR(cd0x, x, 0, 1);
                    }
                    else if (s.i.ecp0.d[0].y > 0) {
                        s.i.ecp0.d[0].y--;
                        EXPAND_VRR(cd0y, y, 0, 1);
                    }
                    else if (s.i.ecp0.d[0].z > 0) {
                        s.i.ecp0.d[0].z--;
                        EXPAND_VRR(cd0z, z, 0, 1);
                    }
                    else if (s.i.ecp0.d[1].x > 0) {
                        s.i.ecp0.d[1].x--;
                        EXPAND_VRR(cd1x, x, 1, 0);
                    }
                    else if (s.i.ecp0.d[1].y > 0) {
                        s.i.ecp0.d[1].y--;
                        EXPAND_VRR(cd1y, y, 1, 0);
                    }
                    else if (s.i.ecp0.d[1].z > 0) {
                        s.i.ecp0.d[1].z--;
                        EXPAND_VRR(cd1z, z, 1, 0);
                    }
                    else if (s.i.ecp0.d[2].x > 0) {
                        s.i.ecp0.d[2].x--;
                        EXPAND_VRR(cdcx, x, 0, 1);
                    }
                    else if (s.i.ecp0.d[2].y > 0) {
                        s.i.ecp0.d[2].y--;
                        EXPAND_VRR(cdcy, y, 0, 1);
                    }
                    else if (s.i.ecp0.d[2].z > 0) {
                        s.i.ecp0.d[2].z--;
                        EXPAND_VRR(cdcz, z, 0, 1);
                    }
                    else {
                        double *v;
                        if (!gtoint__cache__reference_to_store_ecp0(&(itg->c), &s.i, &v)) return GTOINT_ERROR_MEMORY;
                        for (size_t i = 0; i < ng01c; i++) {
                            v[i] = o01c[i];
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
#undef NTERM
#undef NPTR
#undef NVAR
#undef INIT_VRR_COEFFS_A
#undef INIT_VRR_COEFFS_D
#undef INIT_VRR_COEFFS_C
#undef EXPAND_VRR
}
