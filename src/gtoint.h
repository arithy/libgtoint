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

#ifndef GTOINT_INCLUDED
#define GTOINT_INCLUDED

#define GTOINT_NULL ((void *)0)

#ifdef __cplusplus
extern "C" {
#endif

typedef enum gtoint_error_tag {
    GTOINT_ERROR_OK          = 0,
    GTOINT_ERROR_ARGUMENT    = 1,
    GTOINT_ERROR_MEMORY      = 2,
    GTOINT_ERROR_UNSUPPORTED = -2,
    GTOINT_ERROR_INTERNAL    = -1
} gtoint_error_t;

typedef struct gtoint_int3_tag {
    int x, y, z;
} gtoint_int3_t;

typedef struct gtoint_double3_tag {
    double x, y, z;
} gtoint_double3_t;

typedef struct gtoint_integrator_tag *gtoint_integrator_t;
typedef struct gtoint_basis_shell_tag *gtoint_basis_shell_t;
typedef struct gtoint_ecp_shell_tag *gtoint_ecp_shell_t;

gtoint_error_t gtoint_integrator_create(gtoint_integrator_t *itg);
void gtoint_integrator_destroy(gtoint_integrator_t itg);
gtoint_error_t gtoint_integrator_copy(gtoint_integrator_t *itg, gtoint_integrator_t src);
void gtoint_integrator_cleanup_memory(gtoint_integrator_t itg);
void gtoint_integrator_set_error_tolerance(gtoint_integrator_t itg, double tol);
double gtoint_integrator_get_error_tolerance(gtoint_integrator_t itg);
void gtoint_integrator_set_cutoff(gtoint_integrator_t itg, double cut);
double gtoint_integrator_get_cutoff(gtoint_integrator_t itg);

gtoint_error_t gtoint_basis_shell_create(
    gtoint_basis_shell_t *bas, gtoint_integrator_t itg,
    double sf, int na, const int *a, int ng, const double *g, const double *c
);
void gtoint_basis_shell_destroy(gtoint_basis_shell_t bas);
gtoint_error_t gtoint_basis_shell_copy(gtoint_basis_shell_t *bas, gtoint_basis_shell_t src);
int gtoint_basis_shell_get_count(gtoint_basis_shell_t bas);
const int *gtoint_basis_shell_get_angular_numbers(gtoint_basis_shell_t bas, int *na);

gtoint_error_t gtoint_ecp_shell_create(
    gtoint_ecp_shell_t *ecp, gtoint_integrator_t itg,
    int a, int r, int ng, const double *g, const double *c
);
void gtoint_ecp_shell_destroy(gtoint_ecp_shell_t ecp);
gtoint_error_t gtoint_ecp_shell_copy(gtoint_ecp_shell_t *ecp, gtoint_ecp_shell_t src);
int gtoint_ecp_shell_get_angular_number(gtoint_ecp_shell_t ecp);

gtoint_error_t gtoint_compute_overlap_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
    double *out
);
gtoint_error_t gtoint_compute_kinetic_energy_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
    double *out
);
gtoint_error_t gtoint_compute_nuclear_attraction_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    const gtoint_double3_t *pc,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
    double *out
);
gtoint_error_t gtoint_compute_multipole_moment_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    const gtoint_double3_t *pm, int nam, const gtoint_int3_t *am,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
    double *out
);
gtoint_error_t gtoint_compute_electron_repulsion_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    const gtoint_double3_t *p2, gtoint_basis_shell_t bas2,
    const gtoint_double3_t *p3, gtoint_basis_shell_t bas3,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *d2, const gtoint_int3_t *d3,
    double *out
);
gtoint_error_t gtoint_compute_ecp_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
    const gtoint_double3_t *pc, gtoint_ecp_shell_t ecp,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
    double *out
);

gtoint_error_t gtoint_compute_weighted_overlap_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0, const gtoint_double3_t *pw0, double gw0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1, const gtoint_double3_t *pw1, double gw1,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
    double *out
);
gtoint_error_t gtoint_compute_weighted_kinetic_energy_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0, const gtoint_double3_t *pw0, double gw0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1, const gtoint_double3_t *pw1, double gw1,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
    double *out
);
gtoint_error_t gtoint_compute_weighted_nuclear_attraction_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0, const gtoint_double3_t *pw0, double gw0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1, const gtoint_double3_t *pw1, double gw1,
    const gtoint_double3_t *pc,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
    double *out
);
gtoint_error_t gtoint_compute_weighted_multipole_moment_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0, const gtoint_double3_t *pw0, double gw0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1, const gtoint_double3_t *pw1, double gw1,
    const gtoint_double3_t *pm, int nam, const gtoint_int3_t *am,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
    double *out
);
gtoint_error_t gtoint_compute_weighted_electron_repulsion_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0, const gtoint_double3_t *pw0, double gw0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1, const gtoint_double3_t *pw1, double gw1,
    const gtoint_double3_t *p2, gtoint_basis_shell_t bas2, const gtoint_double3_t *pw2, double gw2,
    const gtoint_double3_t *p3, gtoint_basis_shell_t bas3, const gtoint_double3_t *pw3, double gw3,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *d2, const gtoint_int3_t *d3,
    double *out
);
gtoint_error_t gtoint_compute_weighted_ecp_integrals(
    gtoint_integrator_t itg,
    const gtoint_double3_t *p0, gtoint_basis_shell_t bas0, const gtoint_double3_t *pw0, double gw0,
    const gtoint_double3_t *p1, gtoint_basis_shell_t bas1, const gtoint_double3_t *pw1, double gw1,
    const gtoint_double3_t *pc, gtoint_ecp_shell_t ecp,
    int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
    double *out
);

#ifdef __cplusplus
}
#endif

#endif /* !GTOINT_INCLUDED */
