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

#include "stack.h"

#include <stdlib.h>
#include <string.h>

#ifndef ARRAY_INIT_ALLOC
#define ARRAY_INIT_ALLOC 16
#endif

/*== stack_index_array_t ==*/

inline static void stack_index_array__clear_(stack_index_array_t *obj) {
    obj->p = NULL;
}

inline static bool stack_index_array__realloc_(stack_index_array_t *obj, size_t num) {
    stack_index_t *const p = (stack_index_t *)realloc(obj->p, sizeof(stack_index_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void stack_index_array__free_(stack_index_array_t *obj) {
    free(obj->p);
}

void gtoint__stack_index_array__initialize(stack_index_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    stack_index_array__clear_(obj);
}

void gtoint__stack_index_array__finalize(stack_index_array_t *obj) {
    stack_index_array__free_(obj);
}

bool gtoint__stack_index_array__resize(stack_index_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!stack_index_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__stack_index_array__copy(stack_index_array_t *obj, const stack_index_array_t *src) {
    if (!gtoint__stack_index_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(stack_index_t) * src->n);
    return true;
}

void gtoint__stack_index_array__move(stack_index_array_t *obj, stack_index_array_t *src) {
    gtoint__stack_index_array__finalize(obj);
    *obj = *src;
    gtoint__stack_index_array__initialize(src);
}

void gtoint__stack_index_array__compact(stack_index_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n && m != 0) m <<= 1;
        if (m == 0) m = obj->n; /* in case of shift overflow */
        if (m >= obj->m) return;
        stack_index_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        stack_index_array__free_(obj);
        obj->m = 0;
        stack_index_array__clear_(obj);
    }
}

/*== stack_t ==*/

void gtoint__stack__initialize(stack_t *obj) {
    obj->l = 0;
    obj->n = 0;
    gtoint__stack_index_array__initialize(&(obj->i));
    gtoint__double_array__initialize(&(obj->c));
    gtoint__double_array__initialize(&(obj->v));
}

void gtoint__stack__finalize(stack_t *obj) {
    gtoint__stack_index_array__finalize(&(obj->i));
    gtoint__double_array__finalize(&(obj->c));
    gtoint__double_array__finalize(&(obj->v));
}

void gtoint__stack__reset(stack_t *obj, size_t ncc) {
    obj->l = ncc;
    obj->n = 0;
    gtoint__stack_index_array__resize(&(obj->i), 0);
    gtoint__double_array__resize(&(obj->c), 0);
    gtoint__double_array__resize(&(obj->v), 0);
}

bool gtoint__stack__push(stack_t *obj) {
    const size_t n = obj->n;
    if (!gtoint__stack_index_array__resize(&(obj->i), n + 1)) return false;
    if (!gtoint__double_array__resize(&(obj->c), obj->l * (n + 1))) return false;
    if (!gtoint__double_array__resize(&(obj->v), obj->l * (n + 1))) return false;
    obj->n = n + 1;
    return true;
}

bool gtoint__stack__pop(stack_t *obj) {
    if (obj->n <= 0) return false;
    const size_t n = obj->n - 1;
    gtoint__stack_index_array__resize(&(obj->i), n);
    gtoint__double_array__resize(&(obj->c), obj->l * n);
    gtoint__double_array__resize(&(obj->v), obj->l * n);
    obj->n = n;
    return true;
}

bool gtoint__stack__write(stack_t *obj, const stack_index_t *index, size_t *ient) {
    if (obj->n <= 0) return false;
    const size_t n = obj->n - 1;
    obj->i.p[n] = *index;
    if (ient) *ient = n;
    return true;
}

bool gtoint__stack__read(const stack_t *obj, stack_index_t *index, size_t *ient) {
    if (obj->n <= 0) return false;
    const size_t n = obj->n - 1;
    if (index) *index = obj->i.p[n];
    if (ient) *ient = n;
    return true;
}

bool gtoint__stack__copy(stack_t *obj, const stack_t *src) {
    if (!gtoint__stack_index_array__copy(&(obj->i), &(src->i))) return false;
    if (!gtoint__double_array__copy(&(obj->c), &(src->c))) return false;
    if (!gtoint__double_array__copy(&(obj->v), &(src->v))) return false;
    obj->l = src->l;
    obj->n = src->n;
    return true;
}

void gtoint__stack__move(stack_t *obj, stack_t *src) {
    gtoint__stack__finalize(obj);
    *obj = *src;
    gtoint__stack__initialize(src);
}

void gtoint__stack__compact(stack_t *obj) {
    gtoint__stack_index_array__compact(&(obj->i));
    gtoint__double_array__compact(&(obj->c));
    gtoint__double_array__compact(&(obj->v));
}
