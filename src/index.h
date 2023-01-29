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

#ifndef GTOINT_INCLUDED_INDEX_H
#define GTOINT_INCLUDED_INDEX_H

#include "type.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct index_oi_tag {
    int3_t a[2];
    int3_t d[2];
} index_oi_t;

typedef struct index_kei_tag {
    int3_t a[2];
    int3_t d[2];
    int t;
} index_kei_t;

typedef struct index_nai_tag {
    int3_t a[2];
    int3_t d[3];
    int m;
} index_nai_t;

typedef struct index_mmi_tag {
    int3_t a[2];
    int3_t d[2];
    int3_t m;
} index_mmi_t;

typedef struct index_eri_tag {
    int3_t a[4];
    int3_t d[4];
    int m;
} index_eri_t;

typedef struct index_ecp0_tag {
    int3_t a[2];
    int3_t d[3];
} index_ecp0_t;

typedef struct index_ecp1_tag {
    int3_t a[2];
    int3_t d[3];
    int m;
} index_ecp1_t;

typedef struct index_ecp2_tag {
    int3_t a[2];
    int3_t d[3];
    int3_t k[2];
    int r;
} index_ecp2_t;

typedef struct index_ecp2_w_tag {
    int3_t a[2];
    int3_t d[3], dw[2];
    int3_t k[2];
    int r;
} index_ecp2_w_t;

typedef union index_tag {
    index_oi_t oi;
    index_kei_t kei;
    index_nai_t nai;
    index_mmi_t mmi;
    index_eri_t eri;
    index_ecp0_t ecp0;
    index_ecp1_t ecp1;
    index_ecp2_t ecp2;
    index_ecp2_w_t ecp2_w;
} index_t;

#ifdef __cplusplus
}
#endif

#endif /* !GTOINT_INCLUDED_INDEX_H */
