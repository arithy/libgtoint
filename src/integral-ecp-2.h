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

#ifndef GTOINT_INCLUDED_INTEGRAL_ECP_2_H
#define GTOINT_INCLUDED_INTEGRAL_ECP_2_H

#include "type.h"
#include "array.h"
#include "spherical.h"

#define DATABASE_VOID_INDEX (~(~(size_t)0 >> 1))

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ecp_type2_angular_integral_index_tag {
    int l[2];
    int3_t k;
} ecp_type2_angular_integral_index_t;

typedef struct ecp_type2_angular_integral_entry_tag {
    ecp_type2_angular_integral_index_t i;
    double_array_t v; /* [2 * i.l[1] + 1][2 * i.l[0] + 1] */
} ecp_type2_angular_integral_entry_t;

typedef struct ecp_type2_angular_integral_array_tag {
    size_t m, n;
    ecp_type2_angular_integral_entry_t *p;
} ecp_type2_angular_integral_array_t;

typedef struct ecp_type2_angular_integral_database_tag {
    size_t n; /* The number of the entries. */
    size_t_array_t h; /* The hash table. [HASH_TABLE_SIZE] */
    size_t_array_t c; /* The next conflicted entry. [n] */
    ecp_type2_angular_integral_array_t a; /* The entries. [n] */
} ecp_type2_angular_integral_database_t;

typedef struct ecp_type2_spherical_factor_entry_tag {
    int l;
    double_array_t v; /* [2 * l + 1] */
} ecp_type2_spherical_factor_entry_t;

typedef struct ecp_type2_spherical_factor_array_tag {
    size_t m, n;
    ecp_type2_spherical_factor_entry_t *p;
} ecp_type2_spherical_factor_array_t;

typedef struct ecp_type2_spherical_factor_database_tag {
    double3_t r; /* The angular unit vector. */
    size_t n; /* The number of the entries. */
    ecp_type2_spherical_factor_array_t a; /* The entries. */
} ecp_type2_spherical_factor_database_t;

typedef struct ecp_type2_spherical_factor_database_array_tag {
    size_t m, n;
    ecp_type2_spherical_factor_database_t *p;
} ecp_type2_spherical_factor_database_array_t;

typedef struct ecp_type2_radial_integral_index_tag {
    int l[2];
    int k;
} ecp_type2_radial_integral_index_t;

typedef struct ecp_type2_radial_integral_entry_tag {
    ecp_type2_radial_integral_index_t i;
    double v;
} ecp_type2_radial_integral_entry_t;

typedef struct ecp_type2_radial_integral_array_tag {
    size_t m, n;
    ecp_type2_radial_integral_entry_t *p;
} ecp_type2_radial_integral_array_t;

typedef struct ecp_type2_radial_integral_database_tag {
    double gr0, gr1, g0, g1, gc, e0, e1; /* The parameters. */
    double tol; /* The error tolerance. */
    size_t n; /* The number of the entries. */
    size_t_array_t h; /* The hash table. [HASH_TABLE_SIZE] */
    size_t_array_t c; /* The next conflicted entry. [n] */
    ecp_type2_radial_integral_array_t a; /* The entries. [n] */
    double_array_t w; /* The work array. */
} ecp_type2_radial_integral_database_t;

void gtoint__ecp_type2_angular_integral_entry__initialize(ecp_type2_angular_integral_entry_t *obj);
void gtoint__ecp_type2_angular_integral_entry__finalize(ecp_type2_angular_integral_entry_t *obj);
bool gtoint__ecp_type2_angular_integral_entry__resize(ecp_type2_angular_integral_entry_t *obj, size_t num);
bool gtoint__ecp_type2_angular_integral_entry__copy(ecp_type2_angular_integral_entry_t *obj, const ecp_type2_angular_integral_entry_t *src);
void gtoint__ecp_type2_angular_integral_entry__move(ecp_type2_angular_integral_entry_t *obj, ecp_type2_angular_integral_entry_t *src);
void gtoint__ecp_type2_angular_integral_entry__compact(ecp_type2_angular_integral_entry_t *obj);

void gtoint__ecp_type2_angular_integral_array__initialize(ecp_type2_angular_integral_array_t *obj);
void gtoint__ecp_type2_angular_integral_array__finalize(ecp_type2_angular_integral_array_t *obj);
bool gtoint__ecp_type2_angular_integral_array__resize(ecp_type2_angular_integral_array_t *obj, size_t num);
bool gtoint__ecp_type2_angular_integral_array__copy(ecp_type2_angular_integral_array_t *obj, const ecp_type2_angular_integral_array_t *src);
void gtoint__ecp_type2_angular_integral_array__move(ecp_type2_angular_integral_array_t *obj, ecp_type2_angular_integral_array_t *src);
void gtoint__ecp_type2_angular_integral_array__compact(ecp_type2_angular_integral_array_t *obj);

void gtoint__ecp_type2_angular_integral_database__initialize(ecp_type2_angular_integral_database_t *obj);
void gtoint__ecp_type2_angular_integral_database__finalize(ecp_type2_angular_integral_database_t *obj);
void gtoint__ecp_type2_angular_integral_database__clear(ecp_type2_angular_integral_database_t *obj);
bool gtoint__ecp_type2_angular_integral_database__fetch(
    ecp_type2_angular_integral_database_t *obj, spherical_harmonics_database_t *sph, const ecp_type2_angular_integral_index_t *index, double_array_t *out
);
bool gtoint__ecp_type2_angular_integral_database__index(
    ecp_type2_angular_integral_database_t *obj, spherical_harmonics_database_t *sph, const ecp_type2_angular_integral_index_t *index, size_t *out
);
bool gtoint__ecp_type2_angular_integral_database__copy(ecp_type2_angular_integral_database_t *obj, const ecp_type2_angular_integral_database_t *src);
void gtoint__ecp_type2_angular_integral_database__move(ecp_type2_angular_integral_database_t *obj, ecp_type2_angular_integral_database_t *src);
void gtoint__ecp_type2_angular_integral_database__compact(ecp_type2_angular_integral_database_t *obj);

void gtoint__ecp_type2_spherical_factor_entry__initialize(ecp_type2_spherical_factor_entry_t *obj);
void gtoint__ecp_type2_spherical_factor_entry__finalize(ecp_type2_spherical_factor_entry_t *obj);
bool gtoint__ecp_type2_spherical_factor_entry__resize(ecp_type2_spherical_factor_entry_t *obj, size_t num);
bool gtoint__ecp_type2_spherical_factor_entry__copy(ecp_type2_spherical_factor_entry_t *obj, const ecp_type2_spherical_factor_entry_t *src);
void gtoint__ecp_type2_spherical_factor_entry__move(ecp_type2_spherical_factor_entry_t *obj, ecp_type2_spherical_factor_entry_t *src);
void gtoint__ecp_type2_spherical_factor_entry__compact(ecp_type2_spherical_factor_entry_t *obj);

void gtoint__ecp_type2_spherical_factor_array__initialize(ecp_type2_spherical_factor_array_t *obj);
void gtoint__ecp_type2_spherical_factor_array__finalize(ecp_type2_spherical_factor_array_t *obj);
bool gtoint__ecp_type2_spherical_factor_array__resize(ecp_type2_spherical_factor_array_t *obj, size_t num);
bool gtoint__ecp_type2_spherical_factor_array__copy(ecp_type2_spherical_factor_array_t *obj, const ecp_type2_spherical_factor_array_t *src);
void gtoint__ecp_type2_spherical_factor_array__move(ecp_type2_spherical_factor_array_t *obj, ecp_type2_spherical_factor_array_t *src);
void gtoint__ecp_type2_spherical_factor_array__compact(ecp_type2_spherical_factor_array_t *obj);

void gtoint__ecp_type2_spherical_factor_database__initialize(ecp_type2_spherical_factor_database_t *obj);
void gtoint__ecp_type2_spherical_factor_database__finalize(ecp_type2_spherical_factor_database_t *obj);
void gtoint__ecp_type2_spherical_factor_database__clear(ecp_type2_spherical_factor_database_t *obj);
void gtoint__ecp_type2_spherical_factor_database__reset(ecp_type2_spherical_factor_database_t *obj, double3_t vec);
bool gtoint__ecp_type2_spherical_factor_database__fetch(
    ecp_type2_spherical_factor_database_t *obj, spherical_harmonics_database_t *sph, int azim, double_array_t *out
);
bool gtoint__ecp_type2_spherical_factor_database__index(
    ecp_type2_spherical_factor_database_t *obj, spherical_harmonics_database_t *sph, int azim, size_t *out
);
bool gtoint__ecp_type2_spherical_factor_database__copy(ecp_type2_spherical_factor_database_t *obj, const ecp_type2_spherical_factor_database_t *src);
void gtoint__ecp_type2_spherical_factor_database__move(ecp_type2_spherical_factor_database_t *obj, ecp_type2_spherical_factor_database_t *src);
void gtoint__ecp_type2_spherical_factor_database__compact(ecp_type2_spherical_factor_database_t *obj);

void gtoint__ecp_type2_spherical_factor_database_array__initialize(ecp_type2_spherical_factor_database_array_t *obj);
void gtoint__ecp_type2_spherical_factor_database_array__finalize(ecp_type2_spherical_factor_database_array_t *obj);
bool gtoint__ecp_type2_spherical_factor_database_array__resize(ecp_type2_spherical_factor_database_array_t *obj, size_t num);
bool gtoint__ecp_type2_spherical_factor_database_array__copy(ecp_type2_spherical_factor_database_array_t *obj, const ecp_type2_spherical_factor_database_array_t *src);
void gtoint__ecp_type2_spherical_factor_database_array__move(ecp_type2_spherical_factor_database_array_t *obj, ecp_type2_spherical_factor_database_array_t *src);
void gtoint__ecp_type2_spherical_factor_database_array__compact(ecp_type2_spherical_factor_database_array_t *obj);

void gtoint__ecp_type2_radial_integral_array__initialize(ecp_type2_radial_integral_array_t *obj);
void gtoint__ecp_type2_radial_integral_array__finalize(ecp_type2_radial_integral_array_t *obj);
bool gtoint__ecp_type2_radial_integral_array__resize(ecp_type2_radial_integral_array_t *obj, size_t num);
bool gtoint__ecp_type2_radial_integral_array__copy(ecp_type2_radial_integral_array_t *obj, const ecp_type2_radial_integral_array_t *src);
void gtoint__ecp_type2_radial_integral_array__move(ecp_type2_radial_integral_array_t *obj, ecp_type2_radial_integral_array_t *src);
void gtoint__ecp_type2_radial_integral_array__compact(ecp_type2_radial_integral_array_t *obj);

void gtoint__ecp_type2_radial_integral_database__initialize(ecp_type2_radial_integral_database_t *obj);
void gtoint__ecp_type2_radial_integral_database__finalize(ecp_type2_radial_integral_database_t *obj);
void gtoint__ecp_type2_radial_integral_database__clear(ecp_type2_radial_integral_database_t *obj);
void gtoint__ecp_type2_radial_integral_database__reset(ecp_type2_radial_integral_database_t *obj, double3_t r0c, double3_t r1c, double g0, double g1, double gc, double e0, double e1);
bool gtoint__ecp_type2_radial_integral_database__fetch(ecp_type2_radial_integral_database_t *obj, const ecp_type2_radial_integral_index_t *index, double *out);
bool gtoint__ecp_type2_radial_integral_database__copy(ecp_type2_radial_integral_database_t *obj, const ecp_type2_radial_integral_database_t *src);
void gtoint__ecp_type2_radial_integral_database__move(ecp_type2_radial_integral_database_t *obj, ecp_type2_radial_integral_database_t *src);
void gtoint__ecp_type2_radial_integral_database__compact(ecp_type2_radial_integral_database_t *obj);

#ifdef __cplusplus
}
#endif

#endif /* !GTOINT_INCLUDED_INTEGRAL_ECP_2_H */
