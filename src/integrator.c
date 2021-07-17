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

#include <stdlib.h>

gtoint_error_t gtoint_integrator_create(gtoint_integrator_t *itg) {
    struct gtoint_integrator_tag *const p = (struct gtoint_integrator_tag *)malloc(sizeof(struct gtoint_integrator_tag));
    if (p == NULL) return GTOINT_ERROR_MEMORY;
    p->tol = 1e-10;
    p->cut = 1e-15;
    gtoint__spherical_harmonics_database__initialize(&(p->h));
    gtoint__stack__initialize(&(p->s));
    gtoint__cache__initialize(&(p->c));
    gtoint__double_array__initialize(&(p->v));
    gtoint__double_array__initialize(&(p->w));
    gtoint__spherical_harmonics_database__initialize(&(p->ecp.h));
    gtoint__spherical_harmonics_database__reset(&(p->ecp.h), false);
    gtoint__ecp_type2_angular_integral_database__initialize(&(p->ecp.a));
    gtoint__ecp_type2_spherical_factor_database__initialize(&(p->ecp.s[0]));
    gtoint__ecp_type2_spherical_factor_database__initialize(&(p->ecp.s[1]));
    gtoint__ecp_type2_radial_integral_database__initialize(&(p->ecp.r));
    p->ecp.r.tol = p->tol;
    *itg = p;
    return GTOINT_ERROR_OK;
}

void gtoint_integrator_destroy(gtoint_integrator_t itg) {
    if (!itg) return;
    gtoint__spherical_harmonics_database__finalize(&(itg->h));
    gtoint__stack__finalize(&(itg->s));
    gtoint__cache__finalize(&(itg->c));
    gtoint__double_array__finalize(&(itg->v));
    gtoint__double_array__finalize(&(itg->w));
    gtoint__spherical_harmonics_database__finalize(&(itg->ecp.h));
    gtoint__ecp_type2_angular_integral_database__finalize(&(itg->ecp.a));
    gtoint__ecp_type2_spherical_factor_database__finalize(&(itg->ecp.s[0]));
    gtoint__ecp_type2_spherical_factor_database__finalize(&(itg->ecp.s[1]));
    gtoint__ecp_type2_radial_integral_database__finalize(&(itg->ecp.r));
    free(itg);
}

void gtoint_integrator_cleanup_memory(gtoint_integrator_t itg) {
    gtoint__stack__reset(&(itg->s), 0);
    gtoint__stack__compact(&(itg->s));
    gtoint__cache__clear(&(itg->c));
    gtoint__cache__compact(&(itg->c));
    gtoint__double_array__resize(&(itg->v), 0);
    gtoint__double_array__compact(&(itg->v));
    gtoint__double_array__resize(&(itg->w), 0);
    gtoint__double_array__compact(&(itg->w));
    gtoint__ecp_type2_spherical_factor_database__clear(&(itg->ecp.s[0]));
    gtoint__ecp_type2_spherical_factor_database__compact(&(itg->ecp.s[0]));
    gtoint__ecp_type2_spherical_factor_database__clear(&(itg->ecp.s[1]));
    gtoint__ecp_type2_spherical_factor_database__compact(&(itg->ecp.s[1]));
    gtoint__ecp_type2_radial_integral_database__clear(&(itg->ecp.r));
    gtoint__ecp_type2_radial_integral_database__compact(&(itg->ecp.r));
}

void gtoint_integrator_set_error_tolerance(gtoint_integrator_t itg, double tol) {
    itg->tol = tol;
    itg->ecp.r.tol = tol;
}

double gtoint_integrator_get_error_tolerance(gtoint_integrator_t itg) {
    return itg->tol;
}

void gtoint_integrator_set_cutoff(gtoint_integrator_t itg, double cut) {
    itg->cut = cut;
}

double gtoint_integrator_get_cutoff(gtoint_integrator_t itg) {
    return itg->cut;
}
