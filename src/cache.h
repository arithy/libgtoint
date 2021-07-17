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

#ifndef GTOINT_INCLUDED_CACHE_H
#define GTOINT_INCLUDED_CACHE_H

#include "type.h"
#include "index.h"
#include "array.h"

#define CACHE_VOID_INDEX (~(~(size_t)0 >> 1))

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cache_index_array_tag {
    size_t m, n;
    index_t *p;
} cache_index_array_t;

typedef struct cache_tag {
    size_t l; /* The number of the contract combinations. */
    size_t n; /* The number of the entries. */
    size_t_array_t h; /* The hash table. [HASH_TABLE_SIZE] */
    size_t_array_t c; /* The next conflicted entry. [n] */
    cache_index_array_t i; /* The integral indices. [n] */
    double_array_t v; /* The integral values. [n][l] */
} cache_t;

void gtoint__cache_index_array__initialize(cache_index_array_t *obj);
void gtoint__cache_index_array__finalize(cache_index_array_t *obj);
bool gtoint__cache_index_array__resize(cache_index_array_t *obj, size_t num);
bool gtoint__cache_index_array__copy(cache_index_array_t *obj, const cache_index_array_t *src);
void gtoint__cache_index_array__move(cache_index_array_t *obj, cache_index_array_t *src);
void gtoint__cache_index_array__compact(cache_index_array_t *obj);

void gtoint__cache__initialize(cache_t *obj);
void gtoint__cache__finalize(cache_t *obj);
void gtoint__cache__clear(cache_t *obj);
bool gtoint__cache__reset(cache_t *obj, size_t ncc);
bool gtoint__cache__reference_to_store_oi(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_store_kei(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_store_nai(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_store_mmi(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_store_eri(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_store_ecp0(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_store_ecp1(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_store_ecp2(cache_t *obj, const index_t *index, double **value);
bool gtoint__cache__reference_to_fetch_oi(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__reference_to_fetch_kei(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__reference_to_fetch_nai(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__reference_to_fetch_mmi(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__reference_to_fetch_eri(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__reference_to_fetch_ecp0(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__reference_to_fetch_ecp1(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__reference_to_fetch_ecp2(const cache_t *obj, const index_t *index, const double **value);
bool gtoint__cache__copy(cache_t *obj, const cache_t *src);
void gtoint__cache__move(cache_t *obj, cache_t *src);
void gtoint__cache__compact(cache_t *obj);

#ifdef __cplusplus
}
#endif

#endif /* !GTOINT_INCLUDED_CACHE_H */
