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

#ifndef GTOINT_INCLUDED_ARRAY_H
#define GTOINT_INCLUDED_ARRAY_H

#include "type.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct int_array_tag {
    size_t m, n;
    int *p;
} int_array_t;

typedef struct int3_array_tag {
    size_t m, n;
    int3_t *p;
} int3_array_t;

typedef struct size_t_array_tag {
    size_t m, n;
    size_t *p;
} size_t_array_t;

typedef struct double_array_tag {
    size_t m, n;
    double *p;
} double_array_t;

typedef struct double_pointer_array_tag {
    size_t m, n;
    double **p;
} double_pointer_array_t;

void gtoint__int_array__initialize(int_array_t *obj);
void gtoint__int_array__finalize(int_array_t *obj);
bool gtoint__int_array__resize(int_array_t *obj, size_t num);
bool gtoint__int_array__copy(int_array_t *obj, const int_array_t *src);
void gtoint__int_array__move(int_array_t *obj, int_array_t *src);
void gtoint__int_array__compact(int_array_t *obj);

void gtoint__int3_array__initialize(int3_array_t *obj);
void gtoint__int3_array__finalize(int3_array_t *obj);
bool gtoint__int3_array__resize(int3_array_t *obj, size_t num);
bool gtoint__int3_array__copy(int3_array_t *obj, const int3_array_t *src);
void gtoint__int3_array__move(int3_array_t *obj, int3_array_t *src);
void gtoint__int3_array__compact(int3_array_t *obj);

void gtoint__size_t_array__initialize(size_t_array_t *obj);
void gtoint__size_t_array__finalize(size_t_array_t *obj);
bool gtoint__size_t_array__resize(size_t_array_t *obj, size_t num);
bool gtoint__size_t_array__copy(size_t_array_t *obj, const size_t_array_t *src);
void gtoint__size_t_array__move(size_t_array_t *obj, size_t_array_t *src);
void gtoint__size_t_array__compact(size_t_array_t *obj);

void gtoint__double_array__initialize(double_array_t *obj);
void gtoint__double_array__finalize(double_array_t *obj);
bool gtoint__double_array__resize(double_array_t *obj, size_t num);
bool gtoint__double_array__copy(double_array_t *obj, const double_array_t *src);
void gtoint__double_array__move(double_array_t *obj, double_array_t *src);
void gtoint__double_array__compact(double_array_t *obj);

void gtoint__double_pointer_array__initialize(double_pointer_array_t *obj);
void gtoint__double_pointer_array__finalize(double_pointer_array_t *obj);
bool gtoint__double_pointer_array__resize(double_pointer_array_t *obj, size_t num);
bool gtoint__double_pointer_array__copy(double_pointer_array_t *obj, const double_pointer_array_t *src);
void gtoint__double_pointer_array__move(double_pointer_array_t *obj, double_pointer_array_t *src);
void gtoint__double_pointer_array__compact(double_pointer_array_t *obj);

#ifdef __cplusplus
}
#endif

#endif /* !GTOINT_INCLUDED_ARRAY_H */
