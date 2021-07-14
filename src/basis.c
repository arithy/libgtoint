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

#include <stdlib.h>
#include <math.h>

gtoint_error_t gtoint_basis_shell_create(gtoint_basis_shell_t *bas, gtoint_integrator_t itg, double sf, int na, const int *a, int ng, const double *g, const double *c) {
    if (na < 0 || ng < 0) return GTOINT_ERROR_ARGUMENT;
    struct gtoint_basis_shell_tag *const p = (struct gtoint_basis_shell_tag *)malloc(sizeof(struct gtoint_basis_shell_tag));
    if (p == NULL) return GTOINT_ERROR_MEMORY;
    gtoint__int_array__initialize(&(p->spec.azim));
    gtoint__int3_array__initialize(&(p->spec.mom));
    gtoint__spherical_harmonics_array__initialize(&(p->spec.sph));
    gtoint__size_t_array__initialize(&(p->spec.i.m));
    gtoint__size_t_array__initialize(&(p->spec.i.b));
    gtoint__double_array__initialize(&(p->prim.expo));
    gtoint__double_array__initialize(&(p->prim.coef));
    {
        size_t nm = 0;
        for (int ia = 0; ia < na; ia++) {
            int l = a[ia];
            if (l < 0) l = ~l;
            nm += (size_t)(((l + 1) * (l + 2)) >> 1);
        }
        if (!gtoint__int_array__resize(&(p->spec.azim), na)) goto ERROR;
        if (!gtoint__int3_array__resize(&(p->spec.mom), nm)) goto ERROR;
        if (!gtoint__spherical_harmonics_array__resize(&(p->spec.sph), na)) goto ERROR;
        if (!gtoint__size_t_array__resize(&(p->spec.i.m), na + 1)) goto ERROR;
        if (!gtoint__size_t_array__resize(&(p->spec.i.b), na + 1)) goto ERROR;
        if (!gtoint__double_array__resize(&(p->prim.expo), ng)) goto ERROR;
        if (!gtoint__double_array__resize(&(p->prim.coef), ng * nm)) goto ERROR;
        p->spec.n = (size_t)na;
        p->prim.n = (size_t)ng;
    }
    {
        const double f = sf * sf;
        for (size_t ig = 0; ig < p->prim.n; ig++) {
            p->prim.expo.p[ig] = g[ig] * f;
        }
    }
    {
        size_t im = 0, ib = 0;
        for (size_t ia = 0; ia < p->spec.n; ia++) {
            p->spec.azim.p[ia] = a[ia];
            p->spec.i.m.p[ia] = im;
            p->spec.i.b.p[ia] = ib;
            if (a[ia] >= 0) {
                if (!gtoint__spherical_harmonics__makeup_cartesian(&(p->spec.sph.p[ia]), a[ia])) goto ERROR;
                const int3_array_t *const m = &(p->spec.sph.p[ia].a);
                for (size_t jm = 0; jm < m->n; jm++) {
                    p->spec.mom.p[im + jm] = m->p[jm];
                }
                for (size_t ig = 0; ig < p->prim.n; ig++) {
                    for (size_t jm = 0; jm < m->n; jm++) {
                        p->prim.coef.p[ig + p->prim.n * (im + jm)] =
                            c[ig + p->prim.n * ia] * gtoint__compute_cartesian_normalization_constant(m->p[jm], p->prim.expo.p[ig]);
                    }
                }
            }
            else {
                if (!gtoint__spherical_harmonics_database__fetch(&(itg->h), ~a[ia], &(p->spec.sph.p[ia]))) goto ERROR;
                const int3_array_t *const m = &(p->spec.sph.p[ia].a);
                for (size_t jm = 0; jm < m->n; jm++) {
                    p->spec.mom.p[im + jm] = m->p[jm];
                }
                for (size_t ig = 0; ig < p->prim.n; ig++) {
                    const double cn = gtoint__compute_spherical_harmonics_normalization_constant(~a[ia], p->prim.expo.p[ig]);
                    for (size_t jm = 0; jm < m->n; jm++) {
                        p->prim.coef.p[ig + p->prim.n * (im + jm)] = c[ig + p->prim.n * ia] * cn;
                    }
                }
            }
            im += p->spec.sph.p[ia].a.n;
            ib += p->spec.sph.p[ia].n;
        }
        p->spec.i.m.p[p->spec.n] = im;
        p->spec.i.b.p[p->spec.n] = ib;
    }
    {
        const gtoint_error_t e = gtoint__normalize_basis_shell(itg, p);
        if (e != GTOINT_ERROR_OK) {
            gtoint_basis_shell_destroy(p);
            return e;
        }
    }
    *bas = p;
    return GTOINT_ERROR_OK;

ERROR:;
    gtoint_basis_shell_destroy(p);
    return GTOINT_ERROR_MEMORY;
}

void gtoint_basis_shell_destroy(gtoint_basis_shell_t bas) {
    if (!bas) return;
    gtoint__int_array__finalize(&(bas->spec.azim));
    gtoint__int3_array__finalize(&(bas->spec.mom));
    gtoint__spherical_harmonics_array__finalize(&(bas->spec.sph));
    gtoint__size_t_array__finalize(&(bas->spec.i.m));
    gtoint__size_t_array__finalize(&(bas->spec.i.b));
    gtoint__double_array__finalize(&(bas->prim.expo));
    gtoint__double_array__finalize(&(bas->prim.coef));
    free(bas);
}

int gtoint_basis_shell_get_count(gtoint_basis_shell_t bas) {
    return (int)bas->spec.i.b.p[bas->spec.n];
}
