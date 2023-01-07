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

#include "gtoint.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define D_BOHR_A 1.889726124559

static const double ref[] = {
    -0.324380e+00,
     0.000000e+00, -0.533121e+00,
     0.000000e+00,  0.000000e+00, -0.533121e+00,
     0.000000e+00,  0.000000e+00,  0.000000e+00, -0.533121e+00,
    -0.237221e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.164697e-01,
     0.000000e+00, -0.124803e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.148220e-04,
     0.000000e+00,  0.000000e+00, -0.124803e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.148220e-04,
     0.000000e+00,  0.000000e+00,  0.000000e+00, -0.124803e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.148220e-04,
    -0.105729e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.171933e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.186968e+01,
     0.000000e+00, -0.133215e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.284071e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.342380e-01,
     0.000000e+00,  0.000000e+00, -0.133215e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.284071e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.342380e-01,
     0.000000e+00,  0.000000e+00,  0.000000e+00, -0.133215e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.284071e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.342380e-01,
    -0.181662e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.263092e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.285034e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.434683e-01,
     0.000000e+00, -0.590424e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.134864e-04,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.151519e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.671220e-04,
     0.000000e+00,  0.000000e+00, -0.590424e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.134864e-04,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.151519e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.671220e-04,
     0.000000e+00,  0.000000e+00,  0.000000e+00, -0.590424e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.134864e-04,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.151519e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.671220e-04,
    -0.337600e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.100996e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.974990e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.150266e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.143614e+01,
    -0.337600e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.100996e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.974990e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.150266e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.478713e+00, -0.143614e+01,
    -0.337600e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.100996e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.974990e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.150266e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.478713e+00, -0.478713e+00, -0.143614e+01,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.143614e+01,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.143614e+01,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.143614e+01,
    -0.120611e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.232589e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.196290e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.307004e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.341963e-01, -0.113988e-01, -0.113988e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.111306e-02,
    -0.120611e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.232589e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.196290e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.307004e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.113988e-01, -0.341963e-01, -0.113988e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.371021e-03, -0.111306e-02,
    -0.120611e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.232589e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.196290e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.307004e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.113988e-01, -0.113988e-01, -0.341963e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.371021e-03, -0.371021e-03, -0.111306e-02,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.341963e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.111306e-02,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.341963e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.111306e-02,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.341963e-01,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.111306e-02,
    -0.128831e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.245570e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.205903e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.322295e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.360370e-02, -0.120123e-02, -0.120123e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.119135e-03, -0.397115e-04, -0.397115e-04,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.127645e-04,
    -0.128831e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.245570e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.205903e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.322295e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.120123e-02, -0.360370e-02, -0.120123e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.397115e-04, -0.119135e-03, -0.397115e-04,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.425485e-05, -0.127645e-04,
    -0.128831e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.245570e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.205903e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.322295e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.120123e-02, -0.120123e-02, -0.360370e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.397115e-04, -0.397115e-04, -0.119135e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.425485e-05, -0.425485e-05, -0.127645e-04,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.360370e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.119135e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.127645e-04,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.360370e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.119135e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.127645e-04,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.360370e-02,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.119135e-03,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.127645e-04,
     0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,
    -0.223128e-02,  0.269283e-02, -0.269283e-02,  0.000000e+00, -0.280167e-02,  0.637382e-05, -0.637382e-05,  0.000000e+00, -0.301934e-01,  0.690511e-03, -0.690511e-03,  0.000000e+00, -0.460675e-02,  0.306060e-04, -0.306060e-04,  0.000000e+00, -0.163209e-01, -0.163209e-01, -0.158377e-01,  0.836923e-03,  0.000000e+00,  0.000000e+00, -0.342286e-03, -0.342286e-03, -0.326189e-03,  0.278810e-04,  0.000000e+00,  0.000000e+00, -0.359840e-04, -0.359840e-04, -0.342584e-04,  0.298882e-05,  0.000000e+00,  0.000000e+00,  0.000000e+00, -0.517432e-03
};

static const int ord[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 20, 17, 21, 18, 22, 25, 26, 23, 27, 24, 28, 31, 32, 29, 33, 30, 34, 35
};

inline static int symmetric_index(int i, int j) {
    return (i <= j) ? ((j * (j + 1)) >> 1) + i : ((i * (i + 1)) >> 1) + j;
}

int main(int argc, char **argv) {
    int r = 0;
    gtoint_error_t e = GTOINT_ERROR_OK;
    gtoint_integrator_t itg = 0;
    gtoint_basis_shell_t bas_he[2] = { 0 }, bas_fe[7] = { 0 };
    gtoint_ecp_shell_t ecp_fe[1] = { 0 };
    double *val = NULL;
    if ((e = gtoint_integrator_create(&itg)) != GTOINT_ERROR_OK) {
        printf("ERROR: %d: gtoint_integrator_create() -> %d\n", __LINE__, e);
        goto EXIT;
    }
    {
        static const int a[1] = { 0 };
        static const double g[2] = { 13.6267000, 1.9993500 };
        static const double c[2] = { 0.1752300, 0.8934830 };
        if ((e = gtoint_basis_shell_create(&(bas_he[0]), itg, 1.0, 1, a, 2, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 0 };
        static const double g[1] = { 0.3829930 };
        static const double c[1] = { 1.0000000 };
        if ((e = gtoint_basis_shell_create(&(bas_he[1]), itg, 1.0, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[2] = { 0, 1 };
        static const double g[4] = { 70.28000, 6.061000, 4.134000, 1.421000 };
        static const double c[8] = {
            -0.002611, -0.692435, 0.362530, 1.140645,
            -0.007940, -0.290151, 0.591028, 0.719448
        };
        if ((e = gtoint_basis_shell_create(&(bas_fe[0]), itg, 1.0, 2, a, 4, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[2] = { 0, 1 };
        static const double g[2] = { 1.978000, 0.121300 };
        static const double c[4] = {
            -0.098172, 1.026957,
            -0.033731, 1.004462
        };
        if ((e = gtoint_basis_shell_create(&(bas_fe[1]), itg, 1.0, 2, a, 2, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[2] = { 0, 1 };
        static const double g[1] = { 0.512100 };
        static const double c[2] = {
            1.000000,
            1.000000
        };
        if ((e = gtoint_basis_shell_create(&(bas_fe[2]), itg, 1.0, 2, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[2] = { 0, 1 };
        static const double g[1] = { 0.041000 };
        static const double c[2] = {
            1.000000,
            1.000000
        };
        if ((e = gtoint_basis_shell_create(&(bas_fe[3]), itg, 1.0, 2, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 2 };
        static const double g[4] = { 47.10000, 13.12000, 4.478000, 1.581000 };
        static const double c[4] = { 0.026608, 0.152010, 0.413827, 0.605542 };
        if ((e = gtoint_basis_shell_create(&(bas_fe[4]), itg, 1.0, 1, a, 4, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 2 };
        static const double g[1] = { 0.510000 };
        static const double c[1] = { 1.000000 };
        if ((e = gtoint_basis_shell_create(&(bas_fe[5]), itg, 1.0, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 2 };
        static const double g[1] = { 0.138200 };
        static const double c[1] = { 1.000000 };
        if ((e = gtoint_basis_shell_create(&(bas_fe[6]), itg, 1.0, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const double g[1] = { 17.6191000 };
        static const double c[1] = { -3.8942300 };
        if ((e = gtoint_ecp_shell_create(&(ecp_fe[0]), itg, -1, 0, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_ecp_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const gtoint_double3_t p[2] = {
            {  1.6 * D_BOHR_A,  0.4 * D_BOHR_A, 1.0 * D_BOHR_A },
            {  0.4 * D_BOHR_A,  1.6 * D_BOHR_A, 1.0 * D_BOHR_A }
        };
        static const gtoint_int3_t d = { 0, 0, 0 };
        const gtoint_basis_shell_t bas[9] = {
            bas_fe[0], bas_fe[1], bas_fe[2], bas_fe[3], bas_fe[4], bas_fe[5], bas_fe[6],
            bas_he[0], bas_he[1]
        };
        const gtoint_double3_t *const pos_b[9] = {
            &(p[0]), &(p[0]), &(p[0]), &(p[0]), &(p[0]), &(p[0]), &(p[0]),
            &(p[1]), &(p[1])
        };
        const gtoint_ecp_shell_t ecp[1] = {
            ecp_fe[0]
        };
        const gtoint_double3_t *const pos_e[1] = {
            &(p[0])
        };
        for (int i1 = 0, k1 = 0; k1 < 9; k1++) {
            const int n1 = gtoint_basis_shell_get_count(bas[k1]);
        for (int i0 = 0, k0 = 0; k0 <= k1; k0++) {
            const int n0 = gtoint_basis_shell_get_count(bas[k0]);
            {
                double *const t = (double *)realloc(val, sizeof(double) * n0 * n1 * 2);
                if (t == NULL) {
                    printf("ERROR: %d: realloc()\n", __LINE__);
                    goto EXIT;
                }
                val = t;
            }
            for (int i = 0; i < n0 * n1; i++) val[i] = 0.0;
            for (int h = 0; h < 1; h++) {
                if ((e = gtoint_compute_ecp_integrals(
                    itg,
                    pos_b[k0], bas[k0],
                    pos_b[k1], bas[k1],
                    pos_e[h], ecp[h],
                    1, &d, &d, &d,
                    val + n0 * n1
                )) != GTOINT_ERROR_OK) {
                    printf("ERROR: %d: gtoint_compute_ecp_integrals() -> %d\n", __LINE__, e);
                    goto EXIT;
                }
                for (int i = 0; i < n0 * n1; i++) val[i] += val[n0 * n1 + i];
            }
            for (int j1 = 0; j1 < n1; j1++) {
            for (int j0 = 0; j0 < n0; j0++) {
                const double u = ref[symmetric_index(ord[i0 + j0], ord[i1 + j1])];
                const double v = val[j0 + n0 * j1];
                if (!(fabs(v - u) <= 1e-5 * fmax(fabs(u), 1.0))) { /* To catch NaN, negated condition is used. */
                    printf("(%d|%d): %9.6f != %9.6f\n", i0 + j0, i1 + j1, v, u);
                    r = 1;
                }
            }
            }
            i0 += n0;
        }
            i1 += n1;
        }
        printf("Result: %s\n", r ? "FAIL" : "PASS");
    }
EXIT:;
    free(val);
    gtoint_ecp_shell_destroy(ecp_fe[0]);
    gtoint_basis_shell_destroy(bas_fe[6]);
    gtoint_basis_shell_destroy(bas_fe[5]);
    gtoint_basis_shell_destroy(bas_fe[4]);
    gtoint_basis_shell_destroy(bas_fe[3]);
    gtoint_basis_shell_destroy(bas_fe[2]);
    gtoint_basis_shell_destroy(bas_fe[1]);
    gtoint_basis_shell_destroy(bas_fe[0]);
    gtoint_basis_shell_destroy(bas_he[1]);
    gtoint_basis_shell_destroy(bas_he[0]);
    gtoint_integrator_destroy(itg);
    return (r == 0 && e == GTOINT_ERROR_OK) ? 0 : 1;
}
