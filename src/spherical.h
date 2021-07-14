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

#ifndef GTOINT_INCLUDED_SPHERICAL_H
#define GTOINT_INCLUDED_SPHERICAL_H

#include "type.h"
#include "array.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct spherical_harmonics_tag {
    size_t n; /* The number of spherical bases. */
    size_t m; /* The maximum number of the Cartesian bases composing the spherical bases. */
    size_t_array_t l; /* The numbers of the Cartesian bases composing the respective spherical bases. [n] */
    size_t_array_t i; /* The indices of the Cartesian bases composing the respective spherical bases. [n][m] */
    double_array_t c; /* The coefficients of the Cartesian bases composing the respective spherical bases. [n][m] */
    int3_array_t a;   /* The list of the moments. [max(i[:][:])] */
} spherical_harmonics_t;

typedef struct spherical_harmonics_array_tag {
    size_t m, n;
    spherical_harmonics_t *p;
} spherical_harmonics_array_t;

typedef struct spherical_harmonics_database_tag {
    spherical_harmonics_array_t d;
    bool r; /* True if reducing to integers (default). */
} spherical_harmonics_database_t;

void gtoint__spherical_harmonics__initialize(spherical_harmonics_t *obj);
void gtoint__spherical_harmonics__finalize(spherical_harmonics_t *obj);
bool gtoint__spherical_harmonics__compute_reduced(spherical_harmonics_t *obj, int azim);
bool gtoint__spherical_harmonics__compute_normalized(spherical_harmonics_t *obj, int azim);
bool gtoint__spherical_harmonics__makeup_cartesian(spherical_harmonics_t *obj, int azim);
bool gtoint__spherical_harmonics__copy(spherical_harmonics_t *obj, const spherical_harmonics_t *src);
void gtoint__spherical_harmonics__move(spherical_harmonics_t *obj, spherical_harmonics_t *src);
void gtoint__spherical_harmonics__compact(spherical_harmonics_t *obj);

void gtoint__spherical_harmonics_array__initialize(spherical_harmonics_array_t *obj);
void gtoint__spherical_harmonics_array__finalize(spherical_harmonics_array_t *obj);
bool gtoint__spherical_harmonics_array__resize(spherical_harmonics_array_t *obj, size_t num);
bool gtoint__spherical_harmonics_array__copy(spherical_harmonics_array_t *obj, const spherical_harmonics_array_t *src);
void gtoint__spherical_harmonics_array__move(spherical_harmonics_array_t *obj, spherical_harmonics_array_t *src);
void gtoint__spherical_harmonics_array__compact(spherical_harmonics_array_t *obj);

void gtoint__spherical_harmonics_database__initialize(spherical_harmonics_database_t *obj);
void gtoint__spherical_harmonics_database__finalize(spherical_harmonics_database_t *obj);
void gtoint__spherical_harmonics_database__clear(spherical_harmonics_database_t *obj);
void gtoint__spherical_harmonics_database__reset(spherical_harmonics_database_t *obj, bool r);
bool gtoint__spherical_harmonics_database__fetch(spherical_harmonics_database_t *obj, int azim, spherical_harmonics_t *out);
bool gtoint__spherical_harmonics_database__index(spherical_harmonics_database_t *obj, int azim, size_t *out);
bool gtoint__spherical_harmonics_database__copy(spherical_harmonics_database_t *obj, const spherical_harmonics_database_t *src);
void gtoint__spherical_harmonics_database__move(spherical_harmonics_database_t *obj, spherical_harmonics_database_t *src);
void gtoint__spherical_harmonics_database__compact(spherical_harmonics_database_t *obj);

#ifdef __cplusplus
}
#endif

#endif /* !GTOINT_INCLUDED_SPHERICAL_H */
