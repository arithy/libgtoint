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

gtoint_error_t gtoint_compute_ecp_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    const gtoint_double3_t *pc, gtoint_ecp_shell_t ecp,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
    double *out
) {
    if (nd < 0) return GTOINT_ERROR_ARGUMENT;
    const gtoint_error_t e =
        (ecp->azim >= 0) ? gtoint__compute_scalar_ecp_type2_integrals(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p,
            pc, ecp->azim, ecp->rpow, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) :
        (ecp->rpow == 0) ? gtoint__compute_scalar_ecp_type1_integrals_0(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p,
            pc, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) :
        (ecp->rpow == 1) ? gtoint__compute_scalar_ecp_type1_integrals_1(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p,
            pc, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) :
        (ecp->rpow == 2) ? gtoint__compute_scalar_ecp_type1_integrals_2(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p,
            pc, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) : GTOINT_ERROR_ARGUMENT;
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

gtoint_error_t gtoint_compute_weighted_ecp_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0, const gtoint_double3_t *pw0, double gw0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1, const gtoint_double3_t *pw1, double gw1,
    const gtoint_double3_t *pc, gtoint_ecp_shell_t ecp,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
    double *out
) {
    if (nd < 0) return GTOINT_ERROR_ARGUMENT;
    const gtoint_error_t e =
        (ecp->azim >= 0) ? gtoint__compute_weighted_scalar_ecp_type2_integrals(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p, pw0, gw0,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p, pw1, gw1,
            pc, ecp->azim, ecp->rpow, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) :
        (ecp->rpow == 0) ? gtoint__compute_weighted_scalar_ecp_type1_integrals_0(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p, pw0, gw0,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p, pw1, gw1,
            pc, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) :
        (ecp->rpow == 1) ? gtoint__compute_weighted_scalar_ecp_type1_integrals_1(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p, pw0, gw0,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p, pw1, gw1,
            pc, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) :
        (ecp->rpow == 2) ? gtoint__compute_weighted_scalar_ecp_type1_integrals_2(
            itg,
            p0, bas0->spec.mom.n, bas0->spec.mom.p, bas0->prim.n, bas0->prim.expo.p, bas0->prim.coef.p, pw0, gw0,
            p1, bas1->spec.mom.n, bas1->spec.mom.p, bas1->prim.n, bas1->prim.expo.p, bas1->prim.coef.p, pw1, gw1,
            pc, ecp->prim.n, ecp->prim.expo.p, ecp->prim.coef.p,
            nd, d0, d1, dc
        ) : GTOINT_ERROR_ARGUMENT;
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
