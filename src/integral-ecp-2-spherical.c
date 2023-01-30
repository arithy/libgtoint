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

#include "integral-ecp-2.h"

#include <stdlib.h>
#include <assert.h>

#ifndef ARRAY_INIT_ALLOC
#define ARRAY_INIT_ALLOC 16
#endif

inline static double power_(double a, int b) {
    if (b < 0) { b = -b; a = 1.0 / a; }
    double v = 1.0;
    while (b != 0) {
        if (b & 1) v *= a;
        b >>= 1; a *= a;
    }
    return v;
}

/*== ecp_type2_spherical_factor_entry_t ==*/

static bool ecp_type2_spherical_factor_entry__compute_(
    ecp_type2_spherical_factor_entry_t *obj, spherical_harmonics_database_t *sph, int azim, double3_t vec
) {
    assert(!sph->r);
    size_t k;
    if (!gtoint__spherical_harmonics_database__index(sph, azim, &k)) return false;
    const spherical_harmonics_t *const s = &(sph->d.p[k]);
    assert(s->n == (size_t)(2 * azim + 1));
    if (!gtoint__ecp_type2_spherical_factor_entry__resize(obj, s->n)) return false;
    obj->l = azim;
    const size_t l = s->m;
    for (size_t i = 0; i < s->n; i++) {
        double v = 0.0;
        for (size_t j = 0; j < s->l.p[i]; j++) {
            const int3_t m = s->a.p[s->i.p[j + l * i]];
            v += s->c.p[j + l * i] * power_(vec.x, m.x) * power_(vec.y, m.y) * power_(vec.z, m.z);
        }
        obj->v.p[i] = v;
    }
    return true;
}

void gtoint__ecp_type2_spherical_factor_entry__initialize(ecp_type2_spherical_factor_entry_t *obj) {
    gtoint__double_array__initialize(&(obj->v));
}

void gtoint__ecp_type2_spherical_factor_entry__finalize(ecp_type2_spherical_factor_entry_t *obj) {
    gtoint__double_array__finalize(&(obj->v));
}

bool gtoint__ecp_type2_spherical_factor_entry__resize(ecp_type2_spherical_factor_entry_t *obj, size_t num) {
    return gtoint__double_array__resize(&(obj->v), num);
}

bool gtoint__ecp_type2_spherical_factor_entry__copy(ecp_type2_spherical_factor_entry_t *obj, const ecp_type2_spherical_factor_entry_t *src) {
    if (!gtoint__double_array__copy(&(obj->v), &(src->v))) return false;
    obj->l = src->l;
    return true;
}

void gtoint__ecp_type2_spherical_factor_entry__move(ecp_type2_spherical_factor_entry_t *obj, ecp_type2_spherical_factor_entry_t *src) {
    gtoint__ecp_type2_spherical_factor_entry__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_spherical_factor_entry__initialize(src);
}

void gtoint__ecp_type2_spherical_factor_entry__compact(ecp_type2_spherical_factor_entry_t *obj) {
    gtoint__double_array__compact(&(obj->v));
}

/*== ecp_type2_spherical_factor_array_t ==*/

inline static void ecp_type2_spherical_factor_array__clear_(ecp_type2_spherical_factor_array_t *obj) {
    obj->p = NULL;
}

inline static bool ecp_type2_spherical_factor_array__realloc_(ecp_type2_spherical_factor_array_t *obj, size_t num) {
    ecp_type2_spherical_factor_entry_t *const p = (ecp_type2_spherical_factor_entry_t *)realloc(obj->p, sizeof(ecp_type2_spherical_factor_entry_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void ecp_type2_spherical_factor_array__free_(ecp_type2_spherical_factor_array_t *obj) {
    free(obj->p);
}

void gtoint__ecp_type2_spherical_factor_array__initialize(ecp_type2_spherical_factor_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    ecp_type2_spherical_factor_array__clear_(obj);
}

void gtoint__ecp_type2_spherical_factor_array__finalize(ecp_type2_spherical_factor_array_t *obj) {
    for (size_t i = 0; i < obj->n; i++) gtoint__ecp_type2_spherical_factor_entry__finalize(&(obj->p[i]));
    ecp_type2_spherical_factor_array__free_(obj);
}

bool gtoint__ecp_type2_spherical_factor_array__resize(ecp_type2_spherical_factor_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m <= 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!ecp_type2_spherical_factor_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    for (size_t i = num; i < obj->n; i++) gtoint__ecp_type2_spherical_factor_entry__finalize(&(obj->p[i]));
    for (size_t i = obj->n; i < num; i++) gtoint__ecp_type2_spherical_factor_entry__initialize(&(obj->p[i]));
    obj->n = num;
    return true;
}

bool gtoint__ecp_type2_spherical_factor_array__copy(ecp_type2_spherical_factor_array_t *obj, const ecp_type2_spherical_factor_array_t *src) {
    if (!gtoint__ecp_type2_spherical_factor_array__resize(obj, src->n)) return false;
    for (size_t i = 0; i < src->n; i++) {
        if (!gtoint__ecp_type2_spherical_factor_entry__copy(&(obj->p[i]), &(src->p[i]))) return false;
    }
    return true;
}

void gtoint__ecp_type2_spherical_factor_array__move(ecp_type2_spherical_factor_array_t *obj, ecp_type2_spherical_factor_array_t *src) {
    gtoint__ecp_type2_spherical_factor_array__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_spherical_factor_array__initialize(src);
}

void gtoint__ecp_type2_spherical_factor_array__compact(ecp_type2_spherical_factor_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n) m <<= 1;
        if (m >= obj->m) return;
        ecp_type2_spherical_factor_array__realloc_(obj, m);
        obj->m = m;
        for (size_t i = obj->n; i < obj->n; i++) gtoint__ecp_type2_spherical_factor_entry__compact(&(obj->p[i]));
    }
    else {
        ecp_type2_spherical_factor_array__free_(obj);
        obj->m = 0;
        ecp_type2_spherical_factor_array__clear_(obj);
    }
}

/*== ecp_type2_spherical_factor_database_t ==*/

void gtoint__ecp_type2_spherical_factor_database__initialize(ecp_type2_spherical_factor_database_t *obj) {
    obj->r = double3__new(0.0, 0.0, 0.0);
    obj->n = 0;
    gtoint__ecp_type2_spherical_factor_array__initialize(&(obj->a));
}

void gtoint__ecp_type2_spherical_factor_database__finalize(ecp_type2_spherical_factor_database_t *obj) {
    gtoint__ecp_type2_spherical_factor_array__finalize(&(obj->a));
}

void gtoint__ecp_type2_spherical_factor_database__clear(ecp_type2_spherical_factor_database_t *obj) {
    obj->n = 0;
    gtoint__ecp_type2_spherical_factor_array__resize(&(obj->a), 0);
}

void gtoint__ecp_type2_spherical_factor_database__reset(ecp_type2_spherical_factor_database_t *obj, double3_t vec) {
    if (obj->r.x == vec.x && obj->r.y == vec.y && obj->r.z == vec.z) return;
    obj->r = vec;
    obj->n = 0;
    gtoint__ecp_type2_spherical_factor_array__resize(&(obj->a), 0);
}

bool gtoint__ecp_type2_spherical_factor_database__fetch(
    ecp_type2_spherical_factor_database_t *obj, spherical_harmonics_database_t *sph, int azim, double_array_t *out
) {
    if (azim < 0) azim = ~azim;
    if ((size_t)azim >= obj->n) {
        if (!gtoint__ecp_type2_spherical_factor_array__resize(&(obj->a), (size_t)azim + 1)) return false;
    }
    if (obj->a.p[azim].v.n <= 0) {
        if (!ecp_type2_spherical_factor_entry__compute_(&(obj->a.p[azim]), sph, azim, obj->r)) return false;
    }
    if (out) {
        if (!gtoint__double_array__copy(out, &(obj->a.p[azim].v))) return false;
    }
    return true;
}

bool gtoint__ecp_type2_spherical_factor_database__index(
    ecp_type2_spherical_factor_database_t *obj, spherical_harmonics_database_t *sph, int azim, size_t *out
) {
    if (azim < 0) azim = ~azim;
    if ((size_t)azim >= obj->n) {
        if (!gtoint__ecp_type2_spherical_factor_array__resize(&(obj->a), (size_t)azim + 1)) return false;
    }
    if (obj->a.p[azim].v.n <= 0) {
        if (!ecp_type2_spherical_factor_entry__compute_(&(obj->a.p[azim]), sph, azim, obj->r)) return false;
    }
    if (out) *out = (size_t)azim;
    return true;
}

bool gtoint__ecp_type2_spherical_factor_database__copy(ecp_type2_spherical_factor_database_t *obj, const ecp_type2_spherical_factor_database_t *src) {
    if (!gtoint__ecp_type2_spherical_factor_array__copy(&(obj->a), &(src->a))) return false;
    obj->n = src->n;
    obj->r = src->r;
    return true;
}

void gtoint__ecp_type2_spherical_factor_database__move(ecp_type2_spherical_factor_database_t *obj, ecp_type2_spherical_factor_database_t *src) {
    gtoint__ecp_type2_spherical_factor_database__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_spherical_factor_database__initialize(src);
}

void gtoint__ecp_type2_spherical_factor_database__compact(ecp_type2_spherical_factor_database_t *obj) {
    gtoint__ecp_type2_spherical_factor_array__compact(&(obj->a));
}
