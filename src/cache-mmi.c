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

#include "cache.h"

inline static bool equals_mmi_(const index_t *a, const index_t *b) {
    if (a->mmi.m.x != b->mmi.m.x) return false;
    if (a->mmi.m.y != b->mmi.m.y) return false;
    if (a->mmi.m.z != b->mmi.m.z) return false;
    if (a->mmi.a[0].x != b->mmi.a[0].x) return false;
    if (a->mmi.a[0].y != b->mmi.a[0].y) return false;
    if (a->mmi.a[0].z != b->mmi.a[0].z) return false;
    if (a->mmi.a[1].x != b->mmi.a[1].x) return false;
    if (a->mmi.a[1].y != b->mmi.a[1].y) return false;
    if (a->mmi.a[1].z != b->mmi.a[1].z) return false;
    if (a->mmi.d[0].x != b->mmi.d[0].x) return false;
    if (a->mmi.d[0].y != b->mmi.d[0].y) return false;
    if (a->mmi.d[0].z != b->mmi.d[0].z) return false;
    if (a->mmi.d[1].x != b->mmi.d[1].x) return false;
    if (a->mmi.d[1].y != b->mmi.d[1].y) return false;
    if (a->mmi.d[1].z != b->mmi.d[1].z) return false;
    return true;
}

inline static int hash_mmi_(const index_t *a) {
    return
        a->mmi.m.x + 3 * (
        a->mmi.m.y + 3 * (
        a->mmi.m.z + 3 * (
        a->mmi.a[0].x + 5 * (
        a->mmi.a[0].y + 5 * (
        a->mmi.a[0].z + 5 * (
        a->mmi.a[1].x + 5 * (
        a->mmi.a[1].y + 5 * (
        a->mmi.a[1].z + 5 * (
        a->mmi.d[0].x + 3 * (
        a->mmi.d[0].y + 3 * (
        a->mmi.d[0].z + 3 * (
        a->mmi.d[1].x + 3 * (
        a->mmi.d[1].y + 3 * (
        a->mmi.d[1].z))))))))))))));
}

inline static bool cache__find_entry_mmi_(const cache_t *obj, const index_t *index, size_t *out) {
    const size_t h = (size_t)hash_mmi_(index) & (obj->h.n - 1);
    for (size_t i = obj->h.p[h]; i != CACHE_VOID_INDEX; i = obj->c.p[i]) {
        if (equals_mmi_(&(obj->i.p[i]), index)) {
            *out = i;
            return true;
        }
    }
    *out = h;
    return false;
}

bool gtoint__cache__reference_to_store_mmi(cache_t *obj, const index_t *index, double **value) {
    if (obj->h.n <= 0 && !gtoint__cache__reset(obj, obj->l)) return false;
    size_t i;
    if (!cache__find_entry_mmi_(obj, index, &i)) {
        const size_t h = i;
        i = obj->n;
        if (!gtoint__size_t_array__resize(&(obj->c), obj->n + 1)) return false;
        if (!gtoint__cache_index_array__resize(&(obj->i), obj->n + 1)) return false;
        if (!gtoint__double_array__resize(&(obj->v), obj->l * (obj->n + 1))) return false;
        obj->n++;
        obj->i.p[i] = *index;
        obj->c.p[i] = obj->h.p[h];
        obj->h.p[h] = i;
    }
    if (value) *value = obj->v.p + obj->l * i;
    return true;
}

bool gtoint__cache__reference_to_fetch_mmi(const cache_t *obj, const index_t *index, const double **value) {
    size_t i;
    if (!cache__find_entry_mmi_(obj, index, &i)) return false;
    if (value) *value = obj->v.p + obj->l * i;
    return true;
}
