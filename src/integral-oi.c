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

#include <math.h>
#include <assert.h>

#define D_PI 3.1415926535897932384626433832795

static gtoint_error_t compute_overlap_integrals_(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    size_t nd, const int3_t *d0, const int3_t *d1
) {
#define NVAR (1 + 5 * 3 * 4)
#define INIT_VRR_COEFFS_A(out, g1_, r01_, g01_) /* referrence variable: i */ \
    { \
        const double vg1 = (g1_); \
        const double vr01 = (r01_); \
        const double vg01h = 0.5 * (g01_); \
        const double vh1 = vg1 * (g01_); \
        out[0][i] = vh1 * vr01; \
        out[1][i] = -vh1; \
        out[2][i] = vg01h; \
        out[3][i] = vh1; \
        out[4][i] = vg01h; \
    }
#define INIT_VRR_COEFFS_D(out, g0_, r01_, g01_, f01_) /* referrence variable: i */ \
    { \
        const double vg0 = (g0_); \
        const double vr01 = (r01_); \
        const double vh0 = vg0 * (g01_); \
        const double vh1 = vh0 - 1.0; \
        const double vf01d = 2.0 * (f01_); \
        out[0][i] = vf01d * vr01; \
        out[1][i] = -vf01d; \
        out[2][i] = vh1; \
        out[3][i] = vf01d; \
        out[4][i] = vh0; \
    }
#define EXPAND_VRR(coe, xyz, i0, i1) /* referrence variables: itg, e, s, is, ng01 */ \
    { \
        s.o = is; \
        s.b = false; \
        { \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &s, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01; i++) { \
                c[i] = coe[0][i]; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.oi.d[i0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.oi.d[i0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01; i++) { \
                c[i] = coe[1][i] * s.i.oi.d[i0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.oi.a[i0].xyz > 0) { \
            stack_index_t t = s; \
            t.i.oi.a[i0].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01; i++) { \
                c[i] = coe[2][i] * s.i.oi.a[i0].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.oi.d[i1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.oi.d[i1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01; i++) { \
                c[i] = coe[3][i] * s.i.oi.d[i1].xyz; \
                v[i] = 0.0; \
            } \
        } \
        if (s.i.oi.a[i1].xyz > 0) { \
            stack_index_t t = s; \
            t.i.oi.a[i1].xyz--; \
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY; \
            size_t it; \
            gtoint__stack__write(&(itg->s), &t, &it); \
            double *const c = gtoint__stack__coefficients(&(itg->s), it); \
            double *const v = gtoint__stack__integrals(&(itg->s), it); \
            for (size_t i = 0; i < ng01; i++) { \
                c[i] = coe[4][i] * s.i.oi.a[i1].xyz; \
                v[i] = 0.0; \
            } \
        } \
    }
    if (!gtoint__double_array__resize(&(itg->v), na0 * na1 * nd)) return GTOINT_ERROR_MEMORY;
    const double r01x = p1->x - p0->x;
    const double r01y = p1->y - p0->y;
    const double r01z = p1->z - p0->z;
    const double r01w = r01x * r01x + r01y * r01y + r01z * r01z;
    int3_t max_a0 = { 0, 0, 0 };
    int3_t max_a1 = { 0, 0, 0 };
    int3_t max_d0 = { 0, 0, 0 };
    int3_t max_d1 = { 0, 0, 0 };
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
    }
    const size_t ng01 = ng0 * ng1;
    if (!gtoint__cache__reset(&(itg->c), ng01)) return GTOINT_ERROR_MEMORY;
    if (!gtoint__double_array__resize(&(itg->w), ng01 * NVAR)) return GTOINT_ERROR_MEMORY;
    size_t ivar = 0;
    double *const o01 = itg->w.p + ng01 * ivar++;
    double *ca0x[5], *ca0y[5], *ca0z[5];
    double *ca1x[5], *ca1y[5], *ca1z[5];
    double *cd0x[5], *cd0y[5], *cd0z[5];
    double *cd1x[5], *cd1y[5], *cd1z[5];
    for (size_t i = 0; i < 5; i++) ca0x[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) ca0y[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) ca0z[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) ca1x[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) ca1y[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) ca1z[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) cd0x[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) cd0y[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) cd0z[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) cd1x[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) cd1y[i] = itg->w.p + ng01 * ivar++;
    for (size_t i = 0; i < 5; i++) cd1z[i] = itg->w.p + ng01 * ivar++;
    assert(ivar == NVAR);
    for (size_t i1 = 0; i1 < ng1; i1++) {
    for (size_t i0 = 0; i0 < ng0; i0++) {
        const size_t i = i0 + ng0 * i1;
        const double g01 = 1.0 / (g0[i0] + g1[i1]);
        const double h01 = D_PI * g01;
        const double f01 = g0[i0] * g1[i1] * g01;
        o01[i] = h01 * sqrt(h01) * exp(-f01 * r01w);
        if (max_a0.x > 0) {
            INIT_VRR_COEFFS_A(ca0x, g1[i1], r01x, g01);
        }
        if (max_a0.y > 0) {
            INIT_VRR_COEFFS_A(ca0y, g1[i1], r01y, g01);
        }
        if (max_a0.z > 0) {
            INIT_VRR_COEFFS_A(ca0z, g1[i1], r01z, g01);
        }
        if (max_a1.x > 0) {
            INIT_VRR_COEFFS_A(ca1x, g0[i0], -r01x, g01);
        }
        if (max_a1.y > 0) {
            INIT_VRR_COEFFS_A(ca1y, g0[i0], -r01y, g01);
        }
        if (max_a1.z > 0) {
            INIT_VRR_COEFFS_A(ca1z, g0[i0], -r01z, g01);
        }
        if (max_d0.x > 0) {
            INIT_VRR_COEFFS_D(cd0x, g0[i0], r01x, g01, f01);
        }
        if (max_d0.y > 0) {
            INIT_VRR_COEFFS_D(cd0y, g0[i0], r01y, g01, f01);
        }
        if (max_d0.z > 0) {
            INIT_VRR_COEFFS_D(cd0z, g0[i0], r01z, g01, f01);
        }
        if (max_d1.x > 0) {
            INIT_VRR_COEFFS_D(cd1x, g1[i1], -r01x, g01, f01);
        }
        if (max_d1.y > 0) {
            INIT_VRR_COEFFS_D(cd1y, g1[i1], -r01y, g01, f01);
        }
        if (max_d1.z > 0) {
            INIT_VRR_COEFFS_D(cd1z, g1[i1], -r01z, g01, f01);
        }
    }
    }
    size_t io = 0;
    for (size_t id = 0; id < nd; id++) {
    for (size_t ia1 = 0; ia1 < na1; ia1++) {
    for (size_t ia0 = 0; ia0 < na0; ia0++) {
        itg->v.p[io] = 0.0;
        gtoint__stack__reset(&(itg->s), ng01);
        {
            stack_index_t t = { 0 };
            t.i.oi.a[0] = a0[ia0];
            t.i.oi.a[1] = a1[ia1];
            t.i.oi.d[0] = d0[id];
            t.i.oi.d[1] = d1[id];
            t.o = STACK_VOID_INDEX;
            t.b = false;
            if (!gtoint__stack__push(&(itg->s))) return GTOINT_ERROR_MEMORY;
            size_t it;
            gtoint__stack__write(&(itg->s), &t, &it);
            double *const c = gtoint__stack__coefficients(&(itg->s), it);
            double *const v = gtoint__stack__integrals(&(itg->s), it);
            for (size_t i1 = 0; i1 < ng1; i1++) {
            for (size_t i0 = 0; i0 < ng0; i0++) {
                const size_t i = i0 + ng0 * i1;
                c[i] = c0[i0 + ng0 * ia0] * c1[i1 + ng1 * ia1];
                v[i] = 0.0;
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
                if (!gtoint__cache__reference_to_store_oi(&(itg->c), &s.i, &w)) return GTOINT_ERROR_MEMORY;
                for (size_t i = 0; i < ng01; i++) w[i] = v[i];
                if (s.o == STACK_VOID_INDEX) {
                    double o = 0.0;
                    for (size_t i = 0; i < ng01; i++) o += c[i] * v[i];
                    itg->v.p[io] += o;
                }
                else {
                    double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                    for (size_t i = 0; i < ng01; i++) o[i] += c[i] * v[i];
                }
                gtoint__stack__pop(&(itg->s));
            }
            else {
                const double *w;
                if (gtoint__cache__reference_to_fetch_oi(&(itg->c), &s.i, &w)) {
                    const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                    if (s.o == STACK_VOID_INDEX) {
                        double o = 0.0;
                        for (size_t i = 0; i < ng01; i++) o += c[i] * w[i];
                        itg->v.p[io] += o;
                    }
                    else {
                        double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                        for (size_t i = 0; i < ng01; i++) o[i] += c[i] * w[i];
                    }
                    gtoint__stack__pop(&(itg->s));
                }
                else {
                    {
                        s.b = true;
                        gtoint__stack__write(&(itg->s), &s, &is);
                        double *const v = gtoint__stack__integrals(&(itg->s), is);
                        for (size_t i = 0; i < ng01; i++) v[i] = 0.0;
                    }
                    if (s.i.oi.a[0].x > 0) {
                        s.i.oi.a[0].x--;
                        EXPAND_VRR(ca0x, x, 0, 1);
                    }
                    else if (s.i.oi.a[0].y > 0) {
                        s.i.oi.a[0].y--;
                        EXPAND_VRR(ca0y, y, 0, 1);
                    }
                    else if (s.i.oi.a[0].z > 0) {
                        s.i.oi.a[0].z--;
                        EXPAND_VRR(ca0z, z, 0, 1);
                    }
                    else if (s.i.oi.a[1].x > 0) {
                        s.i.oi.a[1].x--;
                        EXPAND_VRR(ca1x, x, 1, 0);
                    }
                    else if (s.i.oi.a[1].y > 0) {
                        s.i.oi.a[1].y--;
                        EXPAND_VRR(ca1y, y, 1, 0);
                    }
                    else if (s.i.oi.a[1].z > 0) {
                        s.i.oi.a[1].z--;
                        EXPAND_VRR(ca1z, z, 1, 0);
                    }
                    else if (s.i.oi.d[0].x > 0) {
                        s.i.oi.d[0].x--;
                        EXPAND_VRR(cd0x, x, 0, 1);
                    }
                    else if (s.i.oi.d[0].y > 0) {
                        s.i.oi.d[0].y--;
                        EXPAND_VRR(cd0y, y, 0, 1);
                    }
                    else if (s.i.oi.d[0].z > 0) {
                        s.i.oi.d[0].z--;
                        EXPAND_VRR(cd0z, z, 0, 1);
                    }
                    else if (s.i.oi.d[1].x > 0) {
                        s.i.oi.d[1].x--;
                        EXPAND_VRR(cd1x, x, 1, 0);
                    }
                    else if (s.i.oi.d[1].y > 0) {
                        s.i.oi.d[1].y--;
                        EXPAND_VRR(cd1y, y, 1, 0);
                    }
                    else if (s.i.oi.d[1].z > 0) {
                        s.i.oi.d[1].z--;
                        EXPAND_VRR(cd1z, z, 1, 0);
                    }
                    else {
                        double *v;
                        if (!gtoint__cache__reference_to_store_oi(&(itg->c), &s.i, &v)) return GTOINT_ERROR_MEMORY;
                        for (size_t i = 0; i < ng01; i++) {
                            v[i] = o01[i];
                        }
                        const double *const c = gtoint__stack__coefficients(&(itg->s), is);
                        if (s.o == STACK_VOID_INDEX) {
                            double o = 0.0;
                            for (size_t i = 0; i < ng01; i++) o += c[i] * v[i];
                            itg->v.p[io] += o;
                        }
                        else {
                            double *const o = gtoint__stack__integrals(&(itg->s), s.o);
                            for (size_t i = 0; i < ng01; i++) o[i] += c[i] * v[i];
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
#undef INIT_VRR_COEFFS_A
#undef INIT_VRR_COEFFS_D
#undef EXPAND_VRR
}

gtoint_error_t gtoint__normalize_basis_shell(gtoint_integrator_t itg, gtoint_basis_shell_t bas) {
    static const double3_t p = { 0.0, 0.0, 0.0 };
    static const int3_t d = { 0, 0, 0 };
    for (size_t ia = 0; ia < bas->spec.n; ia++) {
        const size_t im = bas->spec.i.m.p[ia];
        const size_t nm = bas->spec.i.m.p[ia + 1] - im;
        const size_t ng = bas->prim.n;
        const gtoint_error_t e = compute_overlap_integrals_(
            itg,
            &p, nm, bas->spec.mom.p + im, ng, bas->prim.expo.p, bas->prim.coef.p + ng * im,
            &p, nm, bas->spec.mom.p + im, ng, bas->prim.expo.p, bas->prim.coef.p + ng * im,
            1, &d, &d
        );
        if (e != GTOINT_ERROR_OK) return e;
        spherical_harmonics_t *const s = &(bas->spec.sph.p[ia]);
        for (size_t ib = 0; ib < s->n; ib++) {
            double v = 0.0;
            for (size_t jm1 = 0; jm1 < s->l.p[ib]; jm1++) {
                const size_t k1 = jm1 + s->m * ib;
                double u = 0.0;
                for (size_t jm0 = 0; jm0 < s->l.p[ib]; jm0++) {
                    const size_t k0 = jm0 + s->m * ib;
                    u += s->c.p[k0] * itg->v.p[s->i.p[k0] + nm * s->i.p[k1]];
                }
                v += s->c.p[k1] * u;
            }
            v = 1.0 / sqrt(v);
            for (size_t jm = 0; jm < s->l.p[ib]; jm++) {
                s->c.p[jm + s->m * ib] *= v;
            }
        }
    }
    return GTOINT_ERROR_OK;
}

gtoint_error_t gtoint_compute_overlap_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
    double *out
) {
    if (nd < 0) return GTOINT_ERROR_ARGUMENT;
    const gtoint_error_t e = compute_overlap_integrals_(
        itg,
        p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p,
        p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p,
        nd, d0, d1
    );
    if (e != GTOINT_ERROR_OK) return e;
    const size_t mm0 = bas0->spec.i.m.p[bas0->spec.n];
    const size_t mm1 = bas1->spec.i.m.p[bas1->spec.n];
    const size_t mb0 = bas0->spec.i.b.p[bas0->spec.n];
    const size_t mb1 = bas1->spec.i.b.p[bas1->spec.n];
    for (int id = 0; id < nd; id++) {
        for (size_t ia1 = 0; ia1 < bas1->spec.n; ia1++) {
            const size_t im1 = bas1->spec.i.m.p[ia1];
            const size_t ib1 = bas1->spec.i.b.p[ia1];
            const spherical_harmonics_t *const s1 = &(bas1->spec.sph.p[ia1]);
        for (size_t ia0 = 0; ia0 < bas0->spec.n; ia0++) {
            const size_t im0 = bas0->spec.i.m.p[ia0];
            const size_t ib0 = bas0->spec.i.b.p[ia0];
            const spherical_harmonics_t *const s0 = &(bas0->spec.sph.p[ia0]);
            for (size_t jb1 = 0; jb1 < s1->n; jb1++) {
            for (size_t jb0 = 0; jb0 < s0->n; jb0++) {
                double v1 = 0.0;
                for (size_t jm1 = 0; jm1 < s1->l.p[jb1]; jm1++) {
                    const size_t k1 = jm1 + s1->m * jb1;
                    double v0 = 0.0;
                    for (size_t jm0 = 0; jm0 < s0->l.p[jb0]; jm0++) {
                        const size_t k0 = jm0 + s0->m * jb0;
                        v0 += s0->c.p[k0] * itg->v.p[(im0 + s0->i.p[k0]) + mm0 * ((im1 + s1->i.p[k1]) + mm1 * id)];
                    }
                    v1 += s1->c.p[k1] * v0;
                }
                out[(ib0 + jb0) + mb0 * ((ib1 + jb1) + mb1 * id)] = v1;
            }
            }
        }
        }
    }
    return GTOINT_ERROR_OK;
}
