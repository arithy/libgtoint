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

#include "integral-ecp-2.h"

#include <stdlib.h>
#include <assert.h>

#ifndef ARRAY_INIT_ALLOC
#define ARRAY_INIT_ALLOC 16
#endif

#ifndef HASH_TABLE_SIZE
#define HASH_TABLE_SIZE 0x1000 /* must be a power of 2 */
#endif

#define D_4PI 12.566370614359172953850573533118

inline static double factorial2_(int n) {
#define N 20
    static const double t[N] = {
        1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135,
        645120, 2027025, 10321920, 34459425, 185794560, 654729075
    };
    if (n < 0) return 1.0; /* For -1 and safety. */
    if (n < N) return t[n];
    double v = 1.0;
    for (int i = 2 - (n & 1); i <= n; i += 2) v *= i;
    return v;
#undef N
}

/*== ecp_type2_angular_integral_entry_t ==*/

static bool ecp_type2_angular_integral_entry__compute_(
    ecp_type2_angular_integral_entry_t *obj, spherical_harmonics_database_t *sph, const ecp_type2_angular_integral_index_t *index
) {
    assert(!sph->r);
    size_t k0, k1;
    if (!gtoint__spherical_harmonics_database__index(sph, index->l[0], &k0)) return false;
    if (!gtoint__spherical_harmonics_database__index(sph, index->l[1], &k1)) return false;
    const spherical_harmonics_t *const s0 = &(sph->d.p[k0]);
    const spherical_harmonics_t *const s1 = &(sph->d.p[k1]);
    assert(s0->n == (size_t)(2 * index->l[0] + 1));
    assert(s1->n == (size_t)(2 * index->l[1] + 1));
    if (!gtoint__ecp_type2_angular_integral_entry__resize(obj, s0->n * s1->n)) return false;
    obj->i = *index;
    const int3_t k = obj->i.k;
    const size_t l0 = s0->m;
    const size_t l1 = s1->m;
    size_t i = 0;
    for (size_t i1 = 0; i1 < s1->n; i1++) {
    for (size_t i0 = 0; i0 < s0->n; i0++) {
        double v = 0.0;
        for (size_t j1 = 0; j1 < s1->l.p[i1]; j1++) {
            const int3_t m1 = s1->a.p[s1->i.p[j1 + l1 * i1]];
            double u = 0.0;
            for (size_t j0 = 0; j0 < s0->l.p[i0]; j0++) {
                const int3_t m0 = s0->a.p[s0->i.p[j0 + l0 * i0]];
                const int3_t t = {
                    k.x + m0.x + m1.x,
                    k.y + m0.y + m1.y,
                    k.z + m0.z + m1.z
                };
                if ((t.x & 1) == 0 && (t.y & 1) == 0 && (t.z & 1) == 0) {
                    u += s0->c.p[j0 + l0 * i0] *
                        factorial2_(t.x - 1) * factorial2_(t.y - 1) * factorial2_(t.z - 1) / factorial2_(t.x + t.y + t.z + 1);
                }
            }
            v += s1->c.p[j1 + l1 * i1] * u;
        }
        obj->v.p[i++] = D_4PI * v;
    }
    }
    return true;
}

void gtoint__ecp_type2_angular_integral_entry__initialize(ecp_type2_angular_integral_entry_t *obj) {
    gtoint__double_array__initialize(&(obj->v));
}

void gtoint__ecp_type2_angular_integral_entry__finalize(ecp_type2_angular_integral_entry_t *obj) {
    gtoint__double_array__finalize(&(obj->v));
}

bool gtoint__ecp_type2_angular_integral_entry__resize(ecp_type2_angular_integral_entry_t *obj, size_t num) {
    return gtoint__double_array__resize(&(obj->v), num);
}

bool gtoint__ecp_type2_angular_integral_entry__copy(ecp_type2_angular_integral_entry_t *obj, const ecp_type2_angular_integral_entry_t *src) {
    if (!gtoint__double_array__copy(&(obj->v), &(src->v))) return false;
    obj->i = src->i;
    return true;
}

void gtoint__ecp_type2_angular_integral_entry__move(ecp_type2_angular_integral_entry_t *obj, ecp_type2_angular_integral_entry_t *src) {
    gtoint__ecp_type2_angular_integral_entry__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_angular_integral_entry__initialize(src);
}

void gtoint__ecp_type2_angular_integral_entry__compact(ecp_type2_angular_integral_entry_t *obj) {
    gtoint__double_array__compact(&(obj->v));
}

/*== ecp_type2_angular_integral_array_t ==*/

inline static void ecp_type2_angular_integral_array__clear_(ecp_type2_angular_integral_array_t *obj) {
    obj->p = NULL;
}

inline static bool ecp_type2_angular_integral_array__realloc_(ecp_type2_angular_integral_array_t *obj, size_t num) {
    ecp_type2_angular_integral_entry_t *const p = (ecp_type2_angular_integral_entry_t *)realloc(obj->p, sizeof(ecp_type2_angular_integral_entry_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void ecp_type2_angular_integral_array__free_(ecp_type2_angular_integral_array_t *obj) {
    free(obj->p);
}

void gtoint__ecp_type2_angular_integral_array__initialize(ecp_type2_angular_integral_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    ecp_type2_angular_integral_array__clear_(obj);
}

void gtoint__ecp_type2_angular_integral_array__finalize(ecp_type2_angular_integral_array_t *obj) {
    for (size_t i = 0; i < obj->n; i++) gtoint__ecp_type2_angular_integral_entry__finalize(&(obj->p[i]));
    ecp_type2_angular_integral_array__free_(obj);
}

bool gtoint__ecp_type2_angular_integral_array__resize(ecp_type2_angular_integral_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m <= 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!ecp_type2_angular_integral_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    for (size_t i = num; i < obj->n; i++) gtoint__ecp_type2_angular_integral_entry__finalize(&(obj->p[i]));
    for (size_t i = obj->n; i < num; i++) gtoint__ecp_type2_angular_integral_entry__initialize(&(obj->p[i]));
    obj->n = num;
    return true;
}

bool gtoint__ecp_type2_angular_integral_array__copy(ecp_type2_angular_integral_array_t *obj, const ecp_type2_angular_integral_array_t *src) {
    if (!gtoint__ecp_type2_angular_integral_array__resize(obj, src->n)) return false;
    for (size_t i = 0; i < src->n; i++) {
        if (!gtoint__ecp_type2_angular_integral_entry__copy(&(obj->p[i]), &(src->p[i]))) return false;
    }
    return true;
}

void gtoint__ecp_type2_angular_integral_array__move(ecp_type2_angular_integral_array_t *obj, ecp_type2_angular_integral_array_t *src) {
    gtoint__ecp_type2_angular_integral_array__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_angular_integral_array__initialize(src);
}

void gtoint__ecp_type2_angular_integral_array__compact(ecp_type2_angular_integral_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n) m <<= 1;
        if (m >= obj->m) return;
        ecp_type2_angular_integral_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        ecp_type2_angular_integral_array__free_(obj);
        obj->m = 0;
        ecp_type2_angular_integral_array__clear_(obj);
    }
}

/*== ecp_type2_angular_integral_database_t ==*/

inline static bool equals_(const ecp_type2_angular_integral_index_t *a, const ecp_type2_angular_integral_index_t *b) {
    if (a->k.x != b->k.x) return false;
    if (a->k.y != b->k.y) return false;
    if (a->k.z != b->k.z) return false;
    if (a->l[0] != b->l[0]) return false;
    if (a->l[1] != b->l[1]) return false;
    return true;
}

inline static int hash_(const ecp_type2_angular_integral_index_t *a) {
    return
        a->k.x + 7 * (
        a->k.y + 7 * (
        a->k.z + 7 * (
        a->l[0] + 7 * (
        a->l[1]))));
}

inline static bool ecp_type2_angular_integral_database__find_entry_(const ecp_type2_angular_integral_database_t *obj, const ecp_type2_angular_integral_index_t *index, size_t *out) {
    const size_t h = (size_t)hash_(index) & (obj->h.n - 1);
    for (size_t i = obj->h.p[h]; i != DATABASE_VOID_INDEX; i = obj->c.p[i]) {
        if (equals_(&(obj->a.p[i].i), index)) {
            *out = i;
            return true;
        }
    }
    *out = h;
    return false;
}

void gtoint__ecp_type2_angular_integral_database__initialize(ecp_type2_angular_integral_database_t *obj) {
    obj->n = 0;
    gtoint__size_t_array__initialize(&(obj->h));
    gtoint__size_t_array__initialize(&(obj->c));
    gtoint__ecp_type2_angular_integral_array__initialize(&(obj->a));
}

void gtoint__ecp_type2_angular_integral_database__finalize(ecp_type2_angular_integral_database_t *obj) {
    gtoint__size_t_array__finalize(&(obj->h));
    gtoint__size_t_array__finalize(&(obj->c));
    gtoint__ecp_type2_angular_integral_array__finalize(&(obj->a));
}

void gtoint__ecp_type2_angular_integral_database__clear(ecp_type2_angular_integral_database_t *obj) {
    obj->n = 0;
    gtoint__size_t_array__resize(&(obj->h), 0);
    gtoint__size_t_array__resize(&(obj->c), 0);
    gtoint__ecp_type2_angular_integral_array__resize(&(obj->a), 0);
}

bool gtoint__ecp_type2_angular_integral_database__fetch(
    ecp_type2_angular_integral_database_t *obj, spherical_harmonics_database_t *sph, const ecp_type2_angular_integral_index_t *index, double_array_t *out
) {
    if (obj->h.n <= 0) {
        if (!gtoint__size_t_array__resize(&(obj->h), HASH_TABLE_SIZE)) return false;
        for (size_t i = 0; i < obj->h.n; i++) obj->h.p[i] = DATABASE_VOID_INDEX;
    }
    size_t i;
    if (!ecp_type2_angular_integral_database__find_entry_(obj, index, &i)) {
        const size_t h = i;
        i = obj->n;
        if (!gtoint__size_t_array__resize(&(obj->c), obj->n + 1)) return false;
        if (!gtoint__ecp_type2_angular_integral_array__resize(&(obj->a), obj->n + 1)) return false;
        if (!ecp_type2_angular_integral_entry__compute_(&(obj->a.p[i]), sph, index)) return false;
        obj->n++;
        obj->a.p[i].i = *index;
        obj->c.p[i] = obj->h.p[h];
        obj->h.p[h] = i;
    }
    if (out) {
        if (!gtoint__double_array__copy(out, &(obj->a.p[i].v))) return false;
    }
    return true;
}

bool gtoint__ecp_type2_angular_integral_database__index(
    ecp_type2_angular_integral_database_t *obj, spherical_harmonics_database_t *sph, const ecp_type2_angular_integral_index_t *index, size_t *out
) {
    if (obj->h.n <= 0) {
        if (!gtoint__size_t_array__resize(&(obj->h), HASH_TABLE_SIZE)) return false;
        for (size_t i = 0; i < obj->h.n; i++) obj->h.p[i] = DATABASE_VOID_INDEX;
    }
    size_t i;
    if (!ecp_type2_angular_integral_database__find_entry_(obj, index, &i)) {
        const size_t h = i;
        i = obj->n;
        if (!gtoint__size_t_array__resize(&(obj->c), obj->n + 1)) return false;
        if (!gtoint__ecp_type2_angular_integral_array__resize(&(obj->a), obj->n + 1)) return false;
        if (!ecp_type2_angular_integral_entry__compute_(&(obj->a.p[i]), sph, index)) return false;
        obj->n++;
        obj->a.p[i].i = *index;
        obj->c.p[i] = obj->h.p[h];
        obj->h.p[h] = i;
    }
    if (out) *out = i;
    return true;
}

bool gtoint__ecp_type2_angular_integral_database__copy(ecp_type2_angular_integral_database_t *obj, const ecp_type2_angular_integral_database_t *src) {
    if (!gtoint__size_t_array__copy(&(obj->h), &(src->h))) return false;
    if (!gtoint__size_t_array__copy(&(obj->c), &(src->c))) return false;
    if (!gtoint__ecp_type2_angular_integral_array__copy(&(obj->a), &(src->a))) return false;
    obj->n = src->n;
    return true;
}

void gtoint__ecp_type2_angular_integral_database__move(ecp_type2_angular_integral_database_t *obj, ecp_type2_angular_integral_database_t *src) {
    gtoint__ecp_type2_angular_integral_database__finalize(obj);
    *obj = *src;
    gtoint__ecp_type2_angular_integral_database__initialize(src);
}

void gtoint__ecp_type2_angular_integral_database__compact(ecp_type2_angular_integral_database_t *obj) {
    gtoint__size_t_array__compact(&(obj->h));
    gtoint__size_t_array__compact(&(obj->c));
    gtoint__ecp_type2_angular_integral_array__compact(&(obj->a));
}
