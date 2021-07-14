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

#include "spherical.h"

#include <stdio.h>

#define NOTHING { 0, 0, 0, 0 }

int main(int argc, char **argv) {
    static const double ref_s[1][1][4] = {
        { { 0, 0, 0, 1 } }
    };
    static const double ref_p[3][3][4] = {
        { { 0, 1, 0, 1 }, NOTHING, },
        { { 0, 0, 1, 1 }, NOTHING, },
        { { 1, 0, 0, 1 }, NOTHING, }
    };
    static const double ref_d[5][6][4] = {
        { { 1, 1, 0, 1 }, NOTHING, },
        { { 0, 1, 1, 1 }, NOTHING, },
        { { 2, 0, 0, -1 }, { 0, 2, 0, -1 }, { 0, 0, 2, 2 }, NOTHING, },
        { { 1, 0, 1, 1 }, NOTHING, },
        { { 2, 0, 0, 1 }, { 0, 2, 0, -1 }, NOTHING }
    };
    static const double ref_f[7][10][4] = {
        { { 2, 1, 0, 3 }, { 0, 3, 0, -1 }, NOTHING, },
        { { 1, 1, 1, 1 }, NOTHING, },
        { { 0, 1, 2, 4 }, { 2, 1, 0, -1 }, { 0, 3, 0, -1 }, NOTHING, },
        { { 0, 0, 3, 2 }, { 2, 0, 1, -3 }, { 0, 2, 1, -3 }, NOTHING, },
        { { 1, 0, 2, 4 }, { 3, 0, 0, -1 }, { 1, 2, 0, -1 }, NOTHING, },
        { { 2, 0, 1, 1 }, { 0, 2, 1, -1 }, NOTHING, },
        { { 3, 0, 0, 1 }, { 1, 2, 0, -3 }, NOTHING, }
    };
    static const double ref_g[9][15][4] = {
        { { 3, 1, 0, 1 }, { 1, 3, 0, -1 }, NOTHING, },
        { { 2, 1, 1, 3 }, { 0, 3, 1, -1 }, NOTHING, },
        { { 1, 1, 2, 6 }, { 3, 1, 0, -1 }, { 1, 3, 0, -1 }, NOTHING, },
        { { 0, 1, 3, 4 }, { 2, 1, 1, -3 }, { 0, 3, 1, -3 }, NOTHING, },
        { { 4, 0, 0, 3 }, { 2, 2, 0, 6 }, { 0, 4, 0, 3 }, { 2, 0, 2, -24 }, { 0, 2, 2, -24 }, { 0, 0, 4, 8 }, NOTHING, },
        { { 3, 0, 1, -3 }, { 1, 2, 1, -3 }, { 1, 0, 3, 4 }, NOTHING, },
        { { 4, 0, 0, -1 }, { 0, 4, 0, 1 }, { 2, 0, 2, 6 }, { 0, 2, 2, -6 }, NOTHING, },
        { { 3, 0, 1, 1 }, { 1, 2, 1, -3 }, NOTHING, },
        { { 4, 0, 0, 1 }, { 2, 2, 0, -6 }, { 0, 4, 0, 1 }, NOTHING, }
    };
    static const double (*const ref[5])[4] = {
        ref_s[0], ref_p[0], ref_d[0], ref_f[0], ref_g[0]
    };
    int r = 0;
    bool o = true;
    spherical_harmonics_t c;
    spherical_harmonics_database_t d;
    gtoint__spherical_harmonics__initialize(&c);
    gtoint__spherical_harmonics_database__initialize(&d);
    for (int l = 0; l <= 4; l++) {
        if (!(o = gtoint__spherical_harmonics_database__fetch(&d, l, &c))) goto EXIT;
        if (c.n != (size_t)(l * 2 + 1)) {
            printf("%d: number of spherical bases: %zu != %d\n", l, c.n, l * 2 + 1);
            r = 1;
            continue;
        }
        const int n = ((l + 1) * (l + 2)) >> 1;
        for (int m = -l; m <= l; m++) {
            int k;
            for (k = 0; k < n; k++) {
                if (ref[l][k + n * (l + m)][3] == 0) break;
            }
            if (c.l.p[l + m] != (size_t)k) {
                printf("%d,%2d: number of Cartesian bases: %zu != %d\n", l, m, c.l.p[l + m], k);
                r = 1;
                continue;
            }
            printf("%d,%2d: ", l, m);
            int h = 0;
            for (int ix = l; ix >= 0; ix--) {
            for (int iy = l - ix; iy >= 0; iy--) {
                const int iz = l - ix - iy;
                int i, j;
                for (j = n * (l + m), i = 0; i < n; i++, j++) {
                    if (ref[l][j][3] == 0) { i = n; break; }
                    if (ref[l][j][0] == ix && ref[l][j][1] == iy && ref[l][j][2] == iz) break;
                }
                const double a = (i < n) ? ref[l][j][3] : 0.0;
                for (j = (int)c.m * (l + m), i = 0; i < k; i++, j++) {
                    if (c.i.p[j] == (size_t)h) break;
                }
                const double b = (i < k) ? c.c.p[j] : 0.0;
                if (a == b) {
                    printf(".");
                }
                else {
                    printf("(%d,%d,%d:%g!=%g)", ix, iy, iz, b, a);
                    r = 1;
                }
                h++;
            }
            }
            printf("\n");
        }
    }
    EXIT:;
    gtoint__spherical_harmonics__finalize(&c);
    gtoint__spherical_harmonics_database__finalize(&d);
    return (r == 0 && o) ? 0 : 1;
}
