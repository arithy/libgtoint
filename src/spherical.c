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

#include "spherical.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef ARRAY_INIT_ALLOC
#define ARRAY_INIT_ALLOC 16
#endif

#define D_R2      1.4142135623730950488016887242097  /* sqrt(2) */
#define D_1D2R2PI 0.19947114020071633896997302996719 /* 1 / (2 sqrt(2 PI)) */

inline static int term_index_(int azim, int nx, int ny, int nz) {
    return ((azim + 1 - nx) * (azim - nx) >> 1) + 2 * (azim - nx - ny) - nz;
}

inline static int abs_(int a) {
    return (a >= 0) ? a : -a;
}

inline static double factorial_partial_(int n, int m) {
    double v = 1.0;
    for (int i = m + 1; i <= n; i++) v *= i;
    return v;
}

inline static int binomial_(int n, int k) {
    double v = 1.0;
    for (int i = k + 1; i <= n; i++) v *= i;
    double w = 1.0;
    for (int i = 1; i <= n - k; i++) w *= i;
    return (int)(v / w);
}

inline static int trinomial_(int n, int k, int l) {
    double v = 1.0;
    for (int i = k + 1; i <= n; i++) v *= i;
    double w = 1.0;
    for (int i = 1; i <= l; i++) w *= i;
    for (int i = 1; i <= n - k - l; i++) w *= i;
    return (int)(v / w);
}

inline static int trinomial_incomplete_(int n, int k, int l, int m) {
    double v = 1.0;
    for (int i = k + 1; i <= n; i++) v *= i;
    double w = 1.0;
    for (int i = 1; i <= l; i++) w *= i;
    for (int i = 1; i <= n - k - l - m; i++) w *= i;
    return (int)(v / w);
}

inline static int gcd_(int a, int b) {
    a = abs_(a);
    b = abs_(b);
    while (b != 0) {
        int r = a % b;
        if (r > (b >> 1)) r = b - r;
        a = b;
        b = r;
    }
    return a;
}

inline static void reduce_integers_(int n, int *a) {
    int g = 0;
    for (int i = 0; i < n; i++) g = gcd_(g, a[i]);
    for (int i = 0; i < n; i++) a[i] /= g;
}

/*== spherical_harmonics_t ==*/

void gtoint__spherical_harmonics__initialize(spherical_harmonics_t *obj) {
    obj->n = 0;
    obj->m = 0;
    gtoint__size_t_array__initialize(&(obj->l));
    gtoint__size_t_array__initialize(&(obj->i));
    gtoint__double_array__initialize(&(obj->c));
    gtoint__int3_array__initialize(&(obj->a));
}

void gtoint__spherical_harmonics__finalize(spherical_harmonics_t *obj) {
    gtoint__size_t_array__finalize(&(obj->l));
    gtoint__size_t_array__finalize(&(obj->i));
    gtoint__double_array__finalize(&(obj->c));
    gtoint__int3_array__finalize(&(obj->a));
}

bool gtoint__spherical_harmonics__compute_reduced(spherical_harmonics_t *obj, int azim) {
    static const int r[4] = { 1, 0, -1, 0 };
    bool b = true;
    int_array_t co;
    gtoint__int_array__initialize(&co);
    {
        const int l = (azim >= 0) ? azim : ~azim;
        const int nc = ((l + 1) * (l + 2)) >> 1;
        const int ns = l * 2 + 1;
        if (!(b = gtoint__int_array__resize(&co, (size_t)(nc * ns)))) goto EXIT;
        for (int i = 0; i < nc * ns; i++) co.p[i] = 0;
        { /* m == 0 */
            for (int j = 0; j <= (l >> 1); j++) {
                const int c = ((j & 1) ? -1 : 1) * trinomial_(2 * (l - j), j, l - j);
                for (int ix = 0; ix <= j; ix++) {
                for (int iy = 0; iy <= j - ix; iy++) {
                    const int iz = j - ix - iy;
                    co.p[term_index_(l, 2 * ix, 2 * iy, 2 * iz + l - 2 * j) + nc * l] += c * trinomial_(j, ix, iy);
                }
                }
            }
            reduce_integers_(nc, co.p + nc * l);
        }
        for (int m = 1; m <= l; m++) {
            for (int j = 0; j <= ((l - m) >> 1); j++) {
                const int c = ((j & 1) ? -1 : 1) * trinomial_incomplete_(2 * (l - j), j, l - j, m);
                for (int k = 0; k <= m; k++) {
                    const int d = c * binomial_(m, k);
                    const int dm = d * -r[(m - k + 1) & 3];
                    const int dp = d * r[(m - k) & 3];
                    for (int ix = 0; ix <= j; ix++) {
                    for (int iy = 0; iy <= j - ix; iy++) {
                        const int iz = j - ix - iy;
                        const int t = trinomial_(j, ix, iy);
                        const int i = term_index_(l, 2 * ix + k, 2 * iy + m - k, 2 * iz + l - 2 * j - m);
                        co.p[i + nc * (l - m)] += dm * t;
                        co.p[i + nc * (l + m)] += dp * t;
                    }
                    }
                }
            }
            reduce_integers_(nc, co.p + nc * (l - m));
            reduce_integers_(nc, co.p + nc * (l + m));
        }
        {
            size_t m = 0;
            for (int j = 0; j < ns; j++) {
                size_t n = 0;
                for (int i = 0; i < nc; i++) {
                    if (co.p[i + nc * j] != 0) n++;
                }
                if (m < n) m = n;
            }
            if (!(b = gtoint__size_t_array__resize(&(obj->l), ns))) goto EXIT;
            if (!(b = gtoint__size_t_array__resize(&(obj->i), m * ns))) goto EXIT;
            if (!(b = gtoint__double_array__resize(&(obj->c), m * ns))) goto EXIT;
            if (!(b = gtoint__int3_array__resize(&(obj->a), nc))) goto EXIT;
            obj->n = ns;
            obj->m = m;
            for (int j = 0; j < ns; j++) {
                size_t n = 0;
                for (int i = 0; i < nc; i++) {
                    if (co.p[i + nc * j] != 0) {
                        obj->i.p[n + m * j] = (size_t)i;
                        obj->c.p[n + m * j] = co.p[i + nc * j];
                        n++;
                    }
                }
                obj->l.p[j] = n;
            }
            size_t im = 0;
            for (int ix = l; ix >= 0; ix--) {
            for (int iy = l - ix; iy >= 0; iy--) {
                const int iz = l - ix - iy;
                obj->a.p[im++] = int3__new(ix, iy, iz);
            }
            }
        }
    }
    EXIT:;
    gtoint__int_array__finalize(&co);
    return b;
}

bool gtoint__spherical_harmonics__compute_normalized(spherical_harmonics_t *obj, int azim) {
    static const double r[4] = { 2.0, 0.0, -2.0, 0.0 };
    bool b = true;
    double_array_t co;
    gtoint__double_array__initialize(&co);
    {
        const int l = (azim >= 0) ? azim : ~azim;
        const int nc = ((l + 1) * (l + 2)) >> 1;
        const int ns = l * 2 + 1;
        if (!(b = gtoint__double_array__resize(&co, (size_t)(nc * ns)))) goto EXIT;
        for (int i = 0; i < nc * ns; i++) co.p[i] = 0.0;
        const double g = D_1D2R2PI * exp2(-l) * sqrt(1 + 2 * l);
        { /* m == 0 */
            const double h = g * D_R2;
            for (int j = 0; j <= (l >> 1); j++) {
                const double c = ((j & 1) ? -1 : 1) * trinomial_(2 * (l - j), j, l - j) * h;
                for (int ix = 0; ix <= j; ix++) {
                for (int iy = 0; iy <= j - ix; iy++) {
                    const int iz = j - ix - iy;
                    co.p[term_index_(l, 2 * ix, 2 * iy, 2 * iz + l - 2 * j) + nc * l] += c * trinomial_(j, ix, iy);
                }
                }
            }
        }
        for (int m = 1; m <= l; m++) {
            const double h = g / sqrt(factorial_partial_(l + m, l - m));
            for (int j = 0; j <= ((l - m) >> 1); j++) {
                const double c = ((j & 1) ? -1 : 1) * trinomial_incomplete_(2 * (l - j), j, l - j, m) * h;
                for (int k = 0; k <= m; k++) {
                    const double d = c * binomial_(m, k);
                    const double dm = d * -r[(m - k + 1) & 3];
                    const double dp = d * r[(m - k) & 3];
                    for (int ix = 0; ix <= j; ix++) {
                    for (int iy = 0; iy <= j - ix; iy++) {
                        const int iz = j - ix - iy;
                        const int t = trinomial_(j, ix, iy);
                        const int i = term_index_(l, 2 * ix + k, 2 * iy + m - k, 2 * iz + l - 2 * j - m);
                        co.p[i + nc * (l - m)] += dm * t;
                        co.p[i + nc * (l + m)] += dp * t;
                    }
                    }
                }
            }
        }
        {
            size_t m = 0;
            for (int j = 0; j < ns; j++) {
                size_t n = 0;
                for (int i = 0; i < nc; i++) {
                    if (co.p[i + nc * j] != 0.0) n++;
                }
                if (m < n) m = n;
            }
            if (!(b = gtoint__size_t_array__resize(&(obj->l), ns))) goto EXIT;
            if (!(b = gtoint__size_t_array__resize(&(obj->i), m * ns))) goto EXIT;
            if (!(b = gtoint__double_array__resize(&(obj->c), m * ns))) goto EXIT;
            if (!(b = gtoint__int3_array__resize(&(obj->a), nc))) goto EXIT;
            obj->n = ns;
            obj->m = m;
            for (int j = 0; j < ns; j++) {
                size_t n = 0;
                for (int i = 0; i < nc; i++) {
                    if (co.p[i + nc * j] != 0.0) {
                        obj->i.p[n + m * j] = (size_t)i;
                        obj->c.p[n + m * j] = co.p[i + nc * j];
                        n++;
                    }
                }
                obj->l.p[j] = n;
            }
            size_t im = 0;
            for (int ix = l; ix >= 0; ix--) {
            for (int iy = l - ix; iy >= 0; iy--) {
                const int iz = l - ix - iy;
                obj->a.p[im++] = int3__new(ix, iy, iz);
            }
            }
        }
    }
    EXIT:;
    gtoint__double_array__finalize(&co);
    return b;
}

bool gtoint__spherical_harmonics__makeup_cartesian(spherical_harmonics_t *obj, int azim) {
    const int l = (azim >= 0) ? azim : ~azim;
    const int nc = ((l + 1) * (l + 2)) >> 1;
    if (!gtoint__size_t_array__resize(&(obj->l), nc)) return false;
    if (!gtoint__size_t_array__resize(&(obj->i), nc)) return false;
    if (!gtoint__double_array__resize(&(obj->c), nc)) return false;
    if (!gtoint__int3_array__resize(&(obj->a), nc)) return false;
    obj->n = nc;
    obj->m = 1;
    for (int j = 0; j < nc; j++) {
        obj->l.p[j] = 1;
        obj->i.p[j] = (size_t)j;
        obj->c.p[j] = 1.0;
    }
    size_t im = 0;
    for (int ix = l; ix >= 0; ix--) {
    for (int iy = l - ix; iy >= 0; iy--) {
        const int iz = l - ix - iy;
        obj->a.p[im++] = int3__new(ix, iy, iz);
    }
    }
    return true;
}

bool gtoint__spherical_harmonics__copy(spherical_harmonics_t *obj, const spherical_harmonics_t *src) {
    if (!gtoint__size_t_array__copy(&(obj->l), &(src->l))) return false;
    if (!gtoint__size_t_array__copy(&(obj->i), &(src->i))) return false;
    if (!gtoint__double_array__copy(&(obj->c), &(src->c))) return false;
    if (!gtoint__int3_array__copy(&(obj->a), &(src->a))) return false;
    obj->n = src->n;
    obj->m = src->m;
    return true;
}

void gtoint__spherical_harmonics__move(spherical_harmonics_t *obj, spherical_harmonics_t *src) {
    gtoint__spherical_harmonics__finalize(obj);
    *obj = *src;
    gtoint__spherical_harmonics__initialize(src);
}

void gtoint__spherical_harmonics__compact(spherical_harmonics_t *obj) {
    gtoint__size_t_array__compact(&(obj->l));
    gtoint__size_t_array__compact(&(obj->i));
    gtoint__double_array__compact(&(obj->c));
    gtoint__int3_array__compact(&(obj->a));
}

/*== spherical_harmonics_array_t ==*/

inline static void spherical_harmonics_array__clear_(spherical_harmonics_array_t *obj) {
    obj->p = NULL;
}

inline static bool spherical_harmonics_array__realloc_(spherical_harmonics_array_t *obj, size_t num) {
    spherical_harmonics_t *const p = (spherical_harmonics_t *)realloc(obj->p, sizeof(spherical_harmonics_t) * num);
    if (p == NULL) return false;
    obj->p = p;
    return true;
}

inline static void spherical_harmonics_array__free_(spherical_harmonics_array_t *obj) {
    free(obj->p);
}

void gtoint__spherical_harmonics_array__initialize(spherical_harmonics_array_t *obj) {
    obj->m = 0;
    obj->n = 0;
    spherical_harmonics_array__clear_(obj);
}

void gtoint__spherical_harmonics_array__finalize(spherical_harmonics_array_t *obj) {
    for (size_t i = 0; i < obj->n; i++) gtoint__spherical_harmonics__finalize(&(obj->p[i]));
    spherical_harmonics_array__free_(obj);
}

bool gtoint__spherical_harmonics_array__resize(spherical_harmonics_array_t *obj, size_t num) {
    if (obj->m < num) {
        size_t m = obj->m;
        if (m == 0) m = ARRAY_INIT_ALLOC;
        while (m < num && m != 0) m <<= 1;
        if (m == 0) m = num; /* in case of shift overflow */
        if (!spherical_harmonics_array__realloc_(obj, m)) return false;
        obj->m = m;
    }
    for (size_t i = num; i < obj->n; i++) gtoint__spherical_harmonics__finalize(&(obj->p[i]));
    for (size_t i = obj->n; i < num; i++) gtoint__spherical_harmonics__initialize(&(obj->p[i]));
    obj->n = num;
    return true;
}

bool gtoint__spherical_harmonics_array__copy(spherical_harmonics_array_t *obj, const spherical_harmonics_array_t *src) {
    if (!gtoint__spherical_harmonics_array__resize(obj, src->n)) return false;
    for (size_t i = 0; i < src->n; i++) {
        if (!gtoint__spherical_harmonics__copy(&(obj->p[i]), &(src->p[i]))) return false;
    }
    return true;
}

void gtoint__spherical_harmonics_array__move(spherical_harmonics_array_t *obj, spherical_harmonics_array_t *src) {
    gtoint__spherical_harmonics_array__finalize(obj);
    *obj = *src;
    gtoint__spherical_harmonics_array__initialize(src);
}

void gtoint__spherical_harmonics_array__compact(spherical_harmonics_array_t *obj) {
    if (obj->n > 0) {
        size_t m = ARRAY_INIT_ALLOC;
        while (m < obj->n) m <<= 1;
        if (m >= obj->m) return;
        spherical_harmonics_array__realloc_(obj, m);
        obj->m = m;
        for (size_t i = 0; i < obj->n; i++) gtoint__spherical_harmonics__compact(&(obj->p[i]));
    }
    else {
        spherical_harmonics_array__free_(obj);
        obj->m = 0;
        spherical_harmonics_array__clear_(obj);
    }
}

/*== spherical_harmonics_database_t ==*/

void gtoint__spherical_harmonics_database__initialize(spherical_harmonics_database_t *obj) {
    gtoint__spherical_harmonics_array__initialize(&(obj->d));
    obj->r = true;
}

void gtoint__spherical_harmonics_database__finalize(spherical_harmonics_database_t *obj) {
    gtoint__spherical_harmonics_array__finalize(&(obj->d));
}

void gtoint__spherical_harmonics_database__clear(spherical_harmonics_database_t *obj) {
    gtoint__spherical_harmonics_array__resize(&(obj->d), 0);
}

void gtoint__spherical_harmonics_database__reset(spherical_harmonics_database_t *obj, bool r) {
    if (obj->r == r) return;
    obj->r = r;
    gtoint__spherical_harmonics_array__resize(&(obj->d), 0);
}

bool gtoint__spherical_harmonics_database__fetch(spherical_harmonics_database_t *obj, int azim, spherical_harmonics_t *out) {
    if (azim < 0) azim = ~azim;
    if ((size_t)azim >= obj->d.n) {
        if (!gtoint__spherical_harmonics_array__resize(&(obj->d), (size_t)azim + 1)) return false;
    }
    if (obj->d.p[azim].n <= 0) {
        if (!(obj->r ?
            gtoint__spherical_harmonics__compute_reduced(&(obj->d.p[azim]), azim) :
            gtoint__spherical_harmonics__compute_normalized(&(obj->d.p[azim]), azim)
        )) return false;
    }
    if (out) {
        if (!gtoint__spherical_harmonics__copy(out, &(obj->d.p[azim]))) return false;
    }
    return true;
}

bool gtoint__spherical_harmonics_database__index(spherical_harmonics_database_t *obj, int azim, size_t *out) {
    if (azim < 0) azim = ~azim;
    if ((size_t)azim >= obj->d.n) {
        if (!gtoint__spherical_harmonics_array__resize(&(obj->d), (size_t)azim + 1)) return false;
    }
    if (obj->d.p[azim].n <= 0) {
        if (!(obj->r ?
            gtoint__spherical_harmonics__compute_reduced(&(obj->d.p[azim]), azim) :
            gtoint__spherical_harmonics__compute_normalized(&(obj->d.p[azim]), azim)
        )) return false;
    }
    if (out) *out = (size_t)azim;
    return true;
}

bool gtoint__spherical_harmonics_database__copy(spherical_harmonics_database_t *obj, const spherical_harmonics_database_t *src) {
    if (!gtoint__spherical_harmonics_array__copy(&(obj->d), &(src->d))) return false;
    obj->r = src->r;
    return true;
}

void gtoint__spherical_harmonics_database__move(spherical_harmonics_database_t *obj, spherical_harmonics_database_t *src) {
    gtoint__spherical_harmonics_database__finalize(obj);
    *obj = *src;
    gtoint__spherical_harmonics_database__initialize(src);
}

void gtoint__spherical_harmonics_database__compact(spherical_harmonics_database_t *obj) {
    gtoint__spherical_harmonics_array__compact(&(obj->d));
}
