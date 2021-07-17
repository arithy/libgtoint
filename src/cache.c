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

#include <stdlib.h>
#include <string.h>

#ifndef ARRAY_INIT_ALLOC
#define ARRAY_INIT_ALLOC 16
#endif

#ifndef HASH_TABLE_SIZE
#define HASH_TABLE_SIZE 0x1000 /* must be a power of 2 */
#endif

/*== cache_index_array_t ==*/

inline static void cache_index_array__clear_(cache_index_array_t *obj) {
    obj->p = NULL;
}

inline static bool cache_index_array__realloc_(cache_index_array_t *obj, size_t num) {
    index_t *const p = (index_t *)realloc(obj->p, sizeof(index_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void cache_index_array__free_(cache_index_array_t *obj) {
    free(obj->p);
}

void gtoint__cache_index_array__initialize(cache_index_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    cache_index_array__clear_(obj);
}

void gtoint__cache_index_array__finalize(cache_index_array_t *obj) {
    cache_index_array__free_(obj);
}

bool gtoint__cache_index_array__resize(cache_index_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!cache_index_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__cache_index_array__copy(cache_index_array_t *obj, const cache_index_array_t *src) {
    if (!gtoint__cache_index_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(index_t) * src->n);
    return true;
}

void gtoint__cache_index_array__move(cache_index_array_t *obj, cache_index_array_t *src) {
    gtoint__cache_index_array__finalize(obj);
    *obj = *src;
    gtoint__cache_index_array__initialize(src);
}

void gtoint__cache_index_array__compact(cache_index_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n && m != 0) m <<= 1;
        if (m == 0) m = obj->n; /* in case of shift overflow */
        if (m >= obj->m) return;
        cache_index_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        cache_index_array__free_(obj);
        obj->m = 0;
        cache_index_array__clear_(obj);
    }
}

/*== cache_t ==*/

void gtoint__cache__initialize(cache_t *obj) {
    obj->l = 0;
    obj->n = 0;
    gtoint__size_t_array__initialize(&(obj->h));
    gtoint__size_t_array__initialize(&(obj->c));
    gtoint__cache_index_array__initialize(&(obj->i));
    gtoint__double_array__initialize(&(obj->v));
}

void gtoint__cache__finalize(cache_t *obj) {
    gtoint__size_t_array__finalize(&(obj->h));
    gtoint__size_t_array__finalize(&(obj->c));
    gtoint__cache_index_array__finalize(&(obj->i));
    gtoint__double_array__finalize(&(obj->v));
}

void gtoint__cache__clear(cache_t *obj) {
    obj->l = 0;
    obj->n = 0;
    gtoint__size_t_array__resize(&(obj->h), 0);
    gtoint__size_t_array__resize(&(obj->c), 0);
    gtoint__cache_index_array__resize(&(obj->i), 0);
    gtoint__double_array__resize(&(obj->v), 0);
}

bool gtoint__cache__reset(cache_t *obj, size_t ncc) {
    if (!gtoint__size_t_array__resize(&(obj->h), HASH_TABLE_SIZE)) return false;
    obj->l = ncc;
    obj->n = 0;
    gtoint__size_t_array__resize(&(obj->c), 0);
    gtoint__cache_index_array__resize(&(obj->i), 0);
    gtoint__double_array__resize(&(obj->v), 0);
    for (size_t i = 0; i < obj->h.n; i++) obj->h.p[i] = CACHE_VOID_INDEX;
    return true;
}

bool gtoint__cache__copy(cache_t *obj, const cache_t *src) {
    if (!gtoint__size_t_array__copy(&(obj->h), &(src->h))) return false;
    if (!gtoint__size_t_array__copy(&(obj->c), &(src->c))) return false;
    if (!gtoint__cache_index_array__copy(&(obj->i), &(src->i))) return false;
    if (!gtoint__double_array__copy(&(obj->v), &(src->v))) return false;
    obj->l = src->l;
    obj->n = src->n;
    return true;
}

void gtoint__cache__move(cache_t *obj, cache_t *src) {
    gtoint__cache__finalize(obj);
    *obj = *src;
    gtoint__cache__initialize(src);
}

void gtoint__cache__compact(cache_t *obj) {
    gtoint__size_t_array__compact(&(obj->h));
    gtoint__size_t_array__compact(&(obj->c));
    gtoint__cache_index_array__compact(&(obj->i));
    gtoint__double_array__compact(&(obj->v));
}
