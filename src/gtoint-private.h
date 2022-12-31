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

#ifndef GTOINT_INCLUDED_PRIVATE_H
#define GTOINT_INCLUDED_PRIVATE_H

#include "gtoint.h"

#include "type.h"
#include "array.h"
#include "stack.h"
#include "cache.h"
#include "spherical.h"

#include "integral-ecp-2.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gtoint_integrator_tag {
    double tol; /* The error tolerance of integrals (NAIs and ERIs). */
    double cut; /* The cutoff value (ECP integrals). */
    spherical_harmonics_database_t h; /* The spherical harmonics transform coefficients. */
    stack_t s; /* The stack for expanding integrals. */
    cache_t c; /* The cache to memorize the integral values. */
    double_array_t v; /* The Cartesian integral values */
    double_array_t w; /* The work array. */
    double_pointer_array_t p; /* The array for temporary pointers. */
    struct gtoint_integrator_ecp_tag {
        spherical_harmonics_database_t h;
        ecp_type2_angular_integral_database_t a;
        ecp_type2_spherical_factor_database_t s[2]; /* Must be cleared whenever center positions are changed. */
        ecp_type2_radial_integral_database_t r;     /* Must be cleared whenever center positions and exponents are changed. */
    } ecp;
};

struct gtoint_basis_shell_tag {
    struct basis_shell_spec_tag {
        size_t n; /* The number of azimuths. */
        int_array_t azim; /* The azimuthal numbers (0: s, 1: p, 2: d, ...; bit-inverted if a spherical type). */
        int3_array_t mom; /* The moments. */
        spherical_harmonics_array_t sph; /* The combinations of moments for the respective spherical bases. */
        struct basis_shell_spec_i_tag {
            size_t_array_t m; /* The indices of the first moments in the respective azimuths. */
            size_t_array_t b; /* The indices of the first bases in the respective azimuths. */
        } i;
    } spec;
    struct basis_shell_prim_tag {
        size_t n; /* The number of primitives. */
        double_array_t expo; /* The exponents of primitives applied with the scale factor. [prim.n] */
        double_array_t coef; /* The normalized coefficients of primitives. [spec.mom.n][prim.n] */
    } prim;
};

struct gtoint_ecp_shell_tag {
    int azim; /* The azimuthal number (-1: ul, 0: s, 1: p, 2: d, ...). */
    int rpow; /* The exponent of radii. */
    struct ecp_shell_prim_tag {
        size_t n; /* The number of primitives. */
        double_array_t expo; /* The exponents of primitives. [prim.n] */
        double_array_t coef; /* The coefficients of primitives. [prim.n] */
    } prim;
};

gtoint_error_t gtoint__normalize_basis_shell(gtoint_integrator_t itg, gtoint_basis_shell_t bas);

gtoint_error_t gtoint__compute_scalar_ecp_type1_integrals_0(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    const double3_t *pc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
);
gtoint_error_t gtoint__compute_scalar_ecp_type1_integrals_1(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    const double3_t *pc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
);
gtoint_error_t gtoint__compute_scalar_ecp_type1_integrals_2(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    const double3_t *pc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
);
gtoint_error_t gtoint__compute_scalar_ecp_type2_integrals(
    gtoint_integrator_t itg,
    const double3_t *p0, size_t na0, const int3_t *a0, size_t ng0, const double *g0, const double *c0,
    const double3_t *p1, size_t na1, const int3_t *a1, size_t ng1, const double *g1, const double *c1,
    const double3_t *pc, int ac, int rc, size_t ngc, const double *gc, const double *cc,
    size_t nd, const int3_t *d0, const int3_t *d1, const int3_t *dc
);

#ifdef __cplusplus
}
#endif

#endif /* !GTOINT_INCLUDED_PRIVATE_H */
