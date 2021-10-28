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

#include "gtoint-private.h"

#include <stdlib.h>
#include <string.h>

gtoint_error_t gtoint_ecp_shell_create(gtoint_ecp_shell_t *ecp, gtoint_integrator_t itg, int a, int r, int ng, const double *g, const double *c) {
    if (ng < 0) return GTOINT_ERROR_ARGUMENT;
    struct gtoint_ecp_shell_tag *const p = (struct gtoint_ecp_shell_tag *)malloc(sizeof(struct gtoint_ecp_shell_tag));
    if (p == NULL) return GTOINT_ERROR_MEMORY;
    gtoint__double_array__initialize(&(p->prim.expo));
    gtoint__double_array__initialize(&(p->prim.coef));
    if (!gtoint__double_array__resize(&(p->prim.expo), ng)) goto ERROR;
    if (!gtoint__double_array__resize(&(p->prim.coef), ng)) goto ERROR;
    memcpy(p->prim.expo.p, g, sizeof(double) * (size_t)ng);
    memcpy(p->prim.coef.p, c, sizeof(double) * (size_t)ng);
    p->prim.n = ng;
    p->azim = a;
    p->rpow = r;
    *ecp = p;
    return GTOINT_ERROR_OK;

ERROR:;
    gtoint_ecp_shell_destroy(p);
    return GTOINT_ERROR_MEMORY;
}

void gtoint_ecp_shell_destroy(gtoint_ecp_shell_t ecp) {
    if (ecp == GTOINT_NULL) return;
    gtoint__double_array__finalize(&(ecp->prim.expo));
    gtoint__double_array__finalize(&(ecp->prim.coef));
    free(ecp);
}

gtoint_error_t gtoint_ecp_shell_copy(gtoint_ecp_shell_t *ecp, gtoint_ecp_shell_t src) {
    struct gtoint_ecp_shell_tag *const p = (struct gtoint_ecp_shell_tag *)malloc(sizeof(struct gtoint_ecp_shell_tag));
    if (p == NULL) return GTOINT_ERROR_MEMORY;
    gtoint__double_array__initialize(&(p->prim.expo));
    gtoint__double_array__initialize(&(p->prim.coef));
    if (!gtoint__double_array__copy(&(p->prim.expo), &(src->prim.expo))) goto ERROR;
    if (!gtoint__double_array__copy(&(p->prim.coef), &(src->prim.coef))) goto ERROR;
    p->prim.n = src->prim.n;
    p->azim = src->azim;
    p->rpow = src->rpow;
    *ecp = p;
    return GTOINT_ERROR_OK;

ERROR:;
    gtoint_ecp_shell_destroy(p);
    return GTOINT_ERROR_MEMORY;
}

int gtoint_ecp_shell_get_angular_number(gtoint_ecp_shell_t ecp) {
    return ecp->azim;
}
