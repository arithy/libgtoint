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

#ifndef GTOINT_INCLUDED_STACK_H
#define GTOINT_INCLUDED_STACK_H

#include "type.h"
#include "index.h"
#include "array.h"

#define STACK_VOID_INDEX (~(~(size_t)0 >> 1))

#ifdef __cplusplus
extern "C" {
#endif

typedef struct stack_index_tag {
    index_t i;
    size_t o; /* The index of the stack entry to which integrals multiplied by coefficients are to be added. STACK_VOID_INDEX means the output buffer. */
    bool b; /* True if this entry is for having the final results of the integrals. */
} stack_index_t;

typedef struct stack_index_array_tag {
    size_t m, n;
    stack_index_t *p;
} stack_index_array_t;

typedef struct stack_tag {
    size_t l; /* The number of the contract combinations. */
    size_t n; /* The number of the stack entries. */
    stack_index_array_t i; /* The integral indices. [n] */
    double_array_t c; /* The coefficients. [n][l] */
    double_array_t v; /* The integral values. [n][l] */
} stack_t;

void gtoint__stack_index_array__initialize(stack_index_array_t *obj);
void gtoint__stack_index_array__finalize(stack_index_array_t *obj);
bool gtoint__stack_index_array__resize(stack_index_array_t *obj, size_t num);
bool gtoint__stack_index_array__copy(stack_index_array_t *obj, const stack_index_array_t *src);
void gtoint__stack_index_array__move(stack_index_array_t *obj, stack_index_array_t *src);
void gtoint__stack_index_array__compact(stack_index_array_t *obj);

void gtoint__stack__initialize(stack_t *obj);
void gtoint__stack__finalize(stack_t *obj);
void gtoint__stack__reset(stack_t *obj, size_t ncc);
bool gtoint__stack__push(stack_t *obj);
bool gtoint__stack__pop(stack_t *obj);
bool gtoint__stack__write(stack_t *obj, const stack_index_t *index, size_t *ient);
bool gtoint__stack__read(const stack_t *obj, stack_index_t *index, size_t *ient);
bool gtoint__stack__copy(stack_t *obj, const stack_t *src);
void gtoint__stack__move(stack_t *obj, stack_t *src);
void gtoint__stack__compact(stack_t *obj);

#ifdef __cplusplus
}
#endif

inline static bool gtoint__stack__is_empty(const stack_t *obj) {
    return (obj->n <= 0);
}

inline static double *gtoint__stack__coefficients(stack_t *obj, size_t ient) {
    return obj->c.p + obj->l * ient;
}

inline static double *gtoint__stack__integrals(stack_t *obj, size_t ient) {
    return obj->v.p + obj->l * ient;
}

#endif /* !GTOINT_INCLUDED_STACK_H */
