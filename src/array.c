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

#include "array.h"

#include <stdlib.h>
#include <string.h>

#ifndef ARRAY_INIT_ALLOC
#define ARRAY_INIT_ALLOC 16
#endif

/*== int_array_t ==*/

inline static void int_array__clear_(int_array_t *obj) {
    obj->p = NULL;
}

inline static bool int_array__realloc_(int_array_t *obj, size_t num) {
    int *const p = (int *)realloc(obj->p, sizeof(int) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void int_array__free_(int_array_t *obj) {
    free(obj->p);
}

void gtoint__int_array__initialize(int_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    int_array__clear_(obj);
}

void gtoint__int_array__finalize(int_array_t *obj) {
    int_array__free_(obj);
}

bool gtoint__int_array__resize(int_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!int_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__int_array__copy(int_array_t *obj, const int_array_t *src) {
    if (!gtoint__int_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(int) * src->n);
    return true;
}

void gtoint__int_array__move(int_array_t *obj, int_array_t *src) {
    gtoint__int_array__finalize(obj);
    *obj = *src;
    gtoint__int_array__initialize(src);
}

void gtoint__int_array__compact(int_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n && m != 0) m <<= 1;
        if (m == 0) m = obj->n; /* in case of shift overflow */
        if (m >= obj->m) return;
        int_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        int_array__free_(obj);
        obj->m = 0;
        int_array__clear_(obj);
    }
}

/*== int3_array_t ==*/

inline static void int3_array__clear_(int3_array_t *obj) {
    obj->p = NULL;
}

inline static bool int3_array__realloc_(int3_array_t *obj, size_t num) {
    int3_t *const p = (int3_t *)realloc(obj->p, sizeof(int3_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void int3_array__free_(int3_array_t *obj) {
    free(obj->p);
}

void gtoint__int3_array__initialize(int3_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    int3_array__clear_(obj);
}

void gtoint__int3_array__finalize(int3_array_t *obj) {
    int3_array__free_(obj);
}

bool gtoint__int3_array__resize(int3_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!int3_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__int3_array__copy(int3_array_t *obj, const int3_array_t *src) {
    if (!gtoint__int3_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(int3_t) * src->n);
    return true;
}

void gtoint__int3_array__move(int3_array_t *obj, int3_array_t *src) {
    gtoint__int3_array__finalize(obj);
    *obj = *src;
    gtoint__int3_array__initialize(src);
}

void gtoint__int3_array__compact(int3_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n && m != 0) m <<= 1;
        if (m == 0) m = obj->n; /* in case of shift overflow */
        if (m >= obj->m) return;
        int3_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        int3_array__free_(obj);
        obj->m = 0;
        int3_array__clear_(obj);
    }
}

/*== size_t_array_t ==*/

inline static void size_t_array__clear_(size_t_array_t *obj) {
    obj->p = NULL;
}

inline static bool size_t_array__realloc_(size_t_array_t *obj, size_t num) {
    size_t *const p = (size_t *)realloc(obj->p, sizeof(size_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void size_t_array__free_(size_t_array_t *obj) {
    free(obj->p);
}

void gtoint__size_t_array__initialize(size_t_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    size_t_array__clear_(obj);
}

void gtoint__size_t_array__finalize(size_t_array_t *obj) {
    size_t_array__free_(obj);
}

bool gtoint__size_t_array__resize(size_t_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!size_t_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__size_t_array__copy(size_t_array_t *obj, const size_t_array_t *src) {
    if (!gtoint__size_t_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(size_t) * src->n);
    return true;
}

void gtoint__size_t_array__move(size_t_array_t *obj, size_t_array_t *src) {
    gtoint__size_t_array__finalize(obj);
    *obj = *src;
    gtoint__size_t_array__initialize(src);
}

void gtoint__size_t_array__compact(size_t_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n && m != 0) m <<= 1;
        if (m == 0) m = obj->n; /* in case of shift overflow */
        if (m >= obj->m) return;
        size_t_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        size_t_array__free_(obj);
        obj->m = 0;
        size_t_array__clear_(obj);
    }
}

/*== double_array_t ==*/

inline static void double_array__clear_(double_array_t *obj) {
    obj->p = NULL;
}

inline static bool double_array__realloc_(double_array_t *obj, size_t num) {
    double *const p = (double *)realloc(obj->p, sizeof(double) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void double_array__free_(double_array_t *obj) {
    free(obj->p);
}

void gtoint__double_array__initialize(double_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    double_array__clear_(obj);
}

void gtoint__double_array__finalize(double_array_t *obj) {
    double_array__free_(obj);
}

bool gtoint__double_array__resize(double_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!double_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__double_array__copy(double_array_t *obj, const double_array_t *src) {
    if (!gtoint__double_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(double) * src->n);
    return true;
}

void gtoint__double_array__move(double_array_t *obj, double_array_t *src) {
    gtoint__double_array__finalize(obj);
    *obj = *src;
    gtoint__double_array__initialize(src);
}

void gtoint__double_array__compact(double_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n && m != 0) m <<= 1;
        if (m == 0) m = obj->n; /* in case of shift overflow */
        if (m >= obj->m) return;
        double_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        double_array__free_(obj);
        obj->m = 0;
        double_array__clear_(obj);
    }
}

/*== double_pointer_array_t ==*/

inline static void double_pointer_array__clear_(double_pointer_array_t *obj) {
    obj->p = NULL;
}

inline static bool double_pointer_array__realloc_(double_pointer_array_t *obj, size_t num) {
    double **const p = (double **)realloc(obj->p, sizeof(double *) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void double_pointer_array__free_(double_pointer_array_t *obj) {
    free(obj->p);
}

void gtoint__double_pointer_array__initialize(double_pointer_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    double_pointer_array__clear_(obj);
}

void gtoint__double_pointer_array__finalize(double_pointer_array_t *obj) {
    double_pointer_array__free_(obj);
}

bool gtoint__double_pointer_array__resize(double_pointer_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!double_pointer_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    obj->n = num;
    return true;
}

bool gtoint__double_pointer_array__copy(double_pointer_array_t *obj, const double_pointer_array_t *src) {
    if (!gtoint__double_pointer_array__resize(obj, src->n)) return false;
    memcpy(obj->p, src->p, sizeof(double *) * src->n);
    return true;
}

void gtoint__double_pointer_array__move(double_pointer_array_t *obj, double_pointer_array_t *src) {
    gtoint__double_pointer_array__finalize(obj);
    *obj = *src;
    gtoint__double_pointer_array__initialize(src);
}

void gtoint__double_pointer_array__compact(double_pointer_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n && m != 0) m <<= 1;
        if (m == 0) m = obj->n; /* in case of shift overflow */
        if (m >= obj->m) return;
        double_pointer_array__realloc_(obj, m);
        obj->m = m;
    }
    else {
        double_pointer_array__free_(obj);
        obj->m = 0;
        double_pointer_array__clear_(obj);
    }
}
