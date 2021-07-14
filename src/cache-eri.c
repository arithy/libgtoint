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

#include <string.h>

inline static int compare_indices_eri_(const index_t *a, const index_t *b) {
    if (a->eri.m > b->eri.m) return 1;
    if (a->eri.m < b->eri.m) return -1;
    if (a->eri.a[0].x > b->eri.a[0].x) return 1;
    if (a->eri.a[0].x < b->eri.a[0].x) return -1;
    if (a->eri.a[0].y > b->eri.a[0].y) return 1;
    if (a->eri.a[0].y < b->eri.a[0].y) return -1;
    if (a->eri.a[0].z > b->eri.a[0].z) return 1;
    if (a->eri.a[0].z < b->eri.a[0].z) return -1;
    if (a->eri.a[1].x > b->eri.a[1].x) return 1;
    if (a->eri.a[1].x < b->eri.a[1].x) return -1;
    if (a->eri.a[1].y > b->eri.a[1].y) return 1;
    if (a->eri.a[1].y < b->eri.a[1].y) return -1;
    if (a->eri.a[1].z > b->eri.a[1].z) return 1;
    if (a->eri.a[1].z < b->eri.a[1].z) return -1;
    if (a->eri.a[2].x > b->eri.a[2].x) return 1;
    if (a->eri.a[2].x < b->eri.a[2].x) return -1;
    if (a->eri.a[2].y > b->eri.a[2].y) return 1;
    if (a->eri.a[2].y < b->eri.a[2].y) return -1;
    if (a->eri.a[2].z > b->eri.a[2].z) return 1;
    if (a->eri.a[2].z < b->eri.a[2].z) return -1;
    if (a->eri.a[3].x > b->eri.a[3].x) return 1;
    if (a->eri.a[3].x < b->eri.a[3].x) return -1;
    if (a->eri.a[3].y > b->eri.a[3].y) return 1;
    if (a->eri.a[3].y < b->eri.a[3].y) return -1;
    if (a->eri.a[3].z > b->eri.a[3].z) return 1;
    if (a->eri.a[3].z < b->eri.a[3].z) return -1;
    if (a->eri.d[0].x > b->eri.d[0].x) return 1;
    if (a->eri.d[0].x < b->eri.d[0].x) return -1;
    if (a->eri.d[0].y > b->eri.d[0].y) return 1;
    if (a->eri.d[0].y < b->eri.d[0].y) return -1;
    if (a->eri.d[0].z > b->eri.d[0].z) return 1;
    if (a->eri.d[0].z < b->eri.d[0].z) return -1;
    if (a->eri.d[1].x > b->eri.d[1].x) return 1;
    if (a->eri.d[1].x < b->eri.d[1].x) return -1;
    if (a->eri.d[1].y > b->eri.d[1].y) return 1;
    if (a->eri.d[1].y < b->eri.d[1].y) return -1;
    if (a->eri.d[1].z > b->eri.d[1].z) return 1;
    if (a->eri.d[1].z < b->eri.d[1].z) return -1;
    if (a->eri.d[2].x > b->eri.d[2].x) return 1;
    if (a->eri.d[2].x < b->eri.d[2].x) return -1;
    if (a->eri.d[2].y > b->eri.d[2].y) return 1;
    if (a->eri.d[2].y < b->eri.d[2].y) return -1;
    if (a->eri.d[2].z > b->eri.d[2].z) return 1;
    if (a->eri.d[2].z < b->eri.d[2].z) return -1;
    if (a->eri.d[3].x > b->eri.d[3].x) return 1;
    if (a->eri.d[3].x < b->eri.d[3].x) return -1;
    if (a->eri.d[3].y > b->eri.d[3].y) return 1;
    if (a->eri.d[3].y < b->eri.d[3].y) return -1;
    if (a->eri.d[3].z > b->eri.d[3].z) return 1;
    if (a->eri.d[3].z < b->eri.d[3].z) return -1;
    return 0;
}

inline static bool cache__find_entry_eri_(const cache_t *obj, const index_t *index, size_t *out) {
    if (obj->n <= 0) { *out = 0; return false; }
    size_t min = 0;
    size_t max = obj->n - 1;
    for (;;) {
        const size_t mid = min + ((max - min) >> 1);
        const int c = compare_indices_eri_(index, &(obj->i.p[obj->o.p[mid]]));
        if (c > 0) {
            if (min == mid) { *out = min + 1; return false; }
            min = mid + 1;
        }
        else if (c < 0) {
            if (min == mid) { *out = min; return false; }
            max = mid - 1;
        }
        else {
            *out = mid;
            return true;
        }
    }
}

bool gtoint__cache__reference_to_store_eri(cache_t *obj, const index_t *index, double **value) {
    size_t m, k;
    if (cache__find_entry_eri_(obj, index, &m)) {
        k = obj->o.p[m];
    }
    else {
        k = obj->n;
        if (!gtoint__size_t_array__resize(&(obj->o), obj->n + 1)) return false;
        if (!gtoint__cache_index_array__resize(&(obj->i), obj->n + 1)) return false;
        if (!gtoint__double_array__resize(&(obj->v), obj->l * (obj->n + 1))) return false;
        obj->n++;
        memmove(obj->o.p + m + 1, obj->o.p + m, sizeof(size_t) * (obj->o.n - m - 1));
        obj->o.p[m] = k;
        obj->i.p[k] = *index;
    }
    if (value) *value = obj->v.p + obj->l * k;
    return true;
}

bool gtoint__cache__reference_to_fetch_eri(const cache_t *obj, const index_t *index, const double **value) {
    size_t m;
    if (!cache__find_entry_eri_(obj, index, &m)) return false;
    if (value) *value = obj->v.p + obj->l * obj->o.p[m];
    return true;
}
