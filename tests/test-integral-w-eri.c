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

#include "gtoint.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct input_tag {
    int s[4];
    int p[4];
    gtoint_int3_t d[4];
} input_t;

int main(int argc, char **argv) {
    int r = 0;
    gtoint_error_t e = GTOINT_ERROR_OK;
    gtoint_integrator_t itg = 0;
    gtoint_basis_shell_t bas[6] = { 0 };
    double *val = NULL;
    if ((e = gtoint_integrator_create(&itg)) != GTOINT_ERROR_OK) {
        printf("ERROR: %d: gtoint_integrator_create() -> %d\n", __LINE__, e);
        goto EXIT;
    }
    {
        static const int a[1] = { 0 };
        static const double g[1] = { 1.9 };
        static const double c[1] = { 0.9 };
        static const double f = 1.1;
        if ((e = gtoint_basis_shell_create(&(bas[0]), itg, f, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 0 };
        static const double g[1] = { 1.8 };
        static const double c[1] = { 0.8 };
        static const double f = 1.2;
        if ((e = gtoint_basis_shell_create(&(bas[1]), itg, f, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 1 };
        static const double g[1] = { 1.7 };
        static const double c[1] = { 0.7 };
        static const double f = 1.1;
        if ((e = gtoint_basis_shell_create(&(bas[2]), itg, f, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 1 };
        static const double g[1] = { 1.6 };
        static const double c[1] = { 0.6 };
        static const double f = 1.2;
        if ((e = gtoint_basis_shell_create(&(bas[3]), itg, f, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 2 };
        static const double g[1] = { 1.5 };
        static const double c[1] = { 0.5 };
        static const double f = 1.1;
        if ((e = gtoint_basis_shell_create(&(bas[4]), itg, f, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const int a[1] = { 2 };
        static const double g[1] = { 1.4 };
        static const double c[1] = { 0.4 };
        static const double f = 1.2;
        if ((e = gtoint_basis_shell_create(&(bas[5]), itg, f, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: %d: gtoint_basis_shell_create() -> %d\n", __LINE__, e);
            goto EXIT;
        }
    }
    {
        static const double gw[4] = { 0.6, 0.2, 0.3, 0.9 };
        static const gtoint_double3_t pw[4] = {
            { 0.2, 0.4, 0.5 },
            { 0.3, 0.5, 0.6 },
            { 0.1, 0.3, 0.4 },
            { 0.3, 0.4, 0.2 }
        };
        static const gtoint_double3_t pos[] = {
            { 0.1, 0.2, 0.3 },
            { 0.8, 0.6, 0.4 },
            { 0.5, 0.3, 0.1 },
            { 0.7, 0.4, 0.2 }
        };
        static const input_t in[] = {
            { { 0, 0, 0, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 1, 0, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 0, 2 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 0, 4, 4 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 1, 0, 0, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 0, 2, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 2, 2, 0, 2 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 2, 3, 0, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 2, 0, 4, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 3, 2, 0, 2 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 4, 0, 2, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 4, 2, 2, 0 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 1 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 1, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 1, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 1, 0, 1 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 1, 1 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 1 }, { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 1 }, { 0, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 1, 0, 1 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 1, 1 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 1, 1, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 1, 0, 1 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 1, 1 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 1, 1, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 1, 0, 1 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 1, 1 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 2, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 2, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 2 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 2, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 2, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 2 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 2, 0, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 2, 0 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 2 }, { 0, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 2, 0, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 2, 0 } } },
            { { 0, 2, 3, 1 }, { 0, 1, 2, 3 }, { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 2 } } }
        };
        static const double ref[] = { /* Using results of numerical integration (not exact). */
             0.349830,  0.341226,  0.252904,  0.120905,  0.035255,  0.039438,  0.067540,  0.005859, -0.000362, -0.001304,
             0.040890,  0.148042, -0.019604, -0.036571,  0.064765,  0.008723,  0.071170,  0.005438,  0.125710, -0.005316,
            -0.014406, -0.030154,  0.006941, -0.001731,  0.000735,  0.130811, -0.001118, -0.019092, -0.037070,  0.043678,
             0.000033, -0.000313,  0.126770, -0.000410,  0.049209,  0.000079, -0.000526, -0.000050,  0.000663,  0.127016,
            -0.006253,  0.043216, -0.000239, -0.001100,  0.045372,  0.006403,  0.136935,  0.342820,  0.012900,  0.063874,
             0.104650,  0.063522,  0.145187,  0.027154,  0.086305, -0.071487,  0.013360,  0.022955,  0.019058, -0.135057,
             0.010980, -0.020349,  0.006667,  0.066587, -0.017842,  0.005539,  0.005317, -0.007324, -0.027577, -0.004198,
             0.001132, -0.042590, -0.002504, -0.003019, -0.025402,  0.040954,  0.023023, -0.011804, -0.132485, -0.225131,
            -0.048674, -0.213513,  0.125646, -0.026250, -0.061624, -0.034919,  0.245233,  0.111291,  0.066474,  0.016845,
             0.009890,  0.002180,  0.000510,  0.017330,  0.001827,  0.000763,  0.117936,  0.079553,  0.017494,  0.019338,
             0.025813,  0.010292,  0.129472,  0.076545,  0.035161,  0.069739,  0.151013,  0.028118,  0.091748, -0.068976,
             0.013933,  0.024574,  0.020128, -0.136374,  0.013669, -0.021038,  0.006763,  0.068106, -0.019497,  0.005673,
             0.005795, -0.007276, -0.028171, -0.004281,  0.001059, -0.044726, -0.002900, -0.002958, -0.026684,  0.040190,
             0.022575, -0.012006,  0.063846,  0.032182,  0.005187,  0.010308,  0.002186,  0.007211,  0.061376,  0.079565,
             0.006959,  0.062620,  0.006926,  0.024641,  0.095509,  0.059361,  0.052405,  0.055548,  0.029440,  0.039429,
             0.018396,  0.000289, -0.000991,  0.011479, -0.000425,  0.006706, -0.035834,  0.017277, -0.002595,  0.002762,
             0.001469, -0.003777, -0.012658, -0.005916,  0.034334, -0.000659,  0.006836,  0.002135, -0.010959, -0.026502,
            -0.001104, -0.063091, -0.006813, -0.026007,  0.005587,  0.022905, -0.001154,  0.042358,  0.002902, -0.000585,
            -0.010109, -0.015262,  0.042577, -0.012146,  0.047998,  0.002767, -0.018552, -0.016804, -0.018389, -0.057129,
            -0.029073, -0.041888, -0.053495,  0.034098, -0.029362,  0.012424,  0.015454, -0.023342,  0.022691,  0.014346,
             0.056352,  0.012187,  0.033139,  0.045322, -0.008019, -0.057929, -0.019897, -0.039117,  0.006142, -0.019638,
            -0.067820, -0.106419,  0.035625,  0.019470,  0.013753, -0.001061, -0.137221,  0.029144, -0.028719, -0.103255,
             0.054458,  0.018909, -0.004295, -0.003989,  0.023654, -0.013172, -0.003761,  0.072190, -0.091384, -0.053488,
             0.096268, -0.126877, -0.019611, -0.013869, -0.018838,  0.000279,  0.019267, -0.038718,  0.107075, -0.033047,
            -0.019611, -0.024718,  0.001102,  0.000279, -0.176840, -0.004166,  0.107075, -0.160750, -0.019688, -0.013869,
             0.001102, -0.034015,  0.019267, -0.004166, -0.139946, -0.033047, -0.019688, -0.246319,  0.808880,  0.456133,
             0.134407,  0.008899, -0.003156,  0.000520,  0.019511, -0.004193,  0.000044,  0.008899, -0.003156,  0.000520,
             0.784405,  0.473564,  0.132127, -0.035213,  0.006381,  0.007219,  0.019511, -0.004193,  0.000044, -0.035213,
             0.006381,  0.007219,  0.757482,  0.423554,  0.174605, -0.466206, -0.231188, -0.062287,  0.092665, -0.005443,
             0.013059,  0.147574,  0.067482, -0.023549, -0.003075,  0.002113, -0.000095, -0.425488, -0.226025, -0.060410,
             0.023588,  0.043633, -0.003104,  0.009902,  0.001222, -0.000096, -0.000225,  0.006132,  0.015273, -0.444904,
            -0.247157, -0.029925, -0.008128,  0.043342, -0.016579, -0.086014,  0.042117, -0.045394, -0.056084,  0.080537,
             0.031043, -0.004364, -0.007062,  0.092949, -0.003241, -0.005422,  0.107171, -0.058678, -0.084857,  0.141370,
            -0.003268,  0.004439,  0.020154, -0.019244,  0.006755,  0.122709, -0.078187,  0.035141,  0.081651, -0.107315,
             0.000637, -0.011165, -0.056843, -0.047842,  0.041916, -0.030094, -0.085638, -0.028605, -0.008091,  0.004201,
             0.018433,  0.001628,  0.006213, -0.113036, -0.041057,  0.080351, -0.151583,  0.004201, -0.008582, -0.019898,
             0.006213, -0.019770, -0.008671,  0.080351, -0.130936, -0.084761,  0.115041, -0.063278,  0.013094, -0.020306,
             0.263608,  0.005763,  0.035471,  0.240539,  0.028803,  0.044100,  0.013094, -0.078298,  0.008735,  0.005763,
             0.201529, -0.021010,  0.028803,  0.371020,  0.013094,  0.001938,  0.028946,  0.005763,  0.054971, -0.013372,
             0.028803, -0.051542,  0.218412, -0.704131,  1.060830, -0.119238,  0.013626, -0.007555,  0.000916, -0.018590,
            -0.009473, -0.000043, -0.201498, -0.119238,  1.397355, -0.004570,  0.000916, -0.010541,  0.011914, -0.000043,
            -0.012853, -0.004570,  0.000916, -0.010541, -0.195629, -0.129410,  1.368720, -0.010171,  0.003251, -0.074737,
             0.234576, -0.015560,  0.042574, -0.088737,  0.032865,  0.003561,  0.006995, -0.001166,  0.000111,  0.379747,
             0.222364, -0.077602,  0.006995, -0.001166,  0.000111, -0.078342,  0.032701,  0.001754,  0.006995, -0.001166,
             0.000111,  0.374810,  0.232063, -0.074776,  0.218566, -0.010376,  0.054778, -0.136144, -0.003002, -0.019223,
            -0.919049, -0.454246, -0.132246,  0.043973, -0.008368, -0.007890, -0.219791, -0.108923,  0.031767,  0.043973,
            -0.008368, -0.007890, -0.871207, -0.402234, -0.179456, -0.002571,  0.000457, -0.000063, -0.190041, -0.104012,
             0.029181, -0.144062, -0.232724, -0.028735,  0.090961,  0.116070,  0.010924,  0.466983,  0.210801,  0.059567,
            -0.026977, -0.044300,  0.002625, -0.003950,  0.000175,  0.067618,  0.004491, -0.006401, -0.015602,  0.494557,
             0.234007,  0.025703, -0.000045, -0.000664,  0.000012,  0.007713,  0.000793,  0.064167,  0.086807,  0.123696,
             0.004563, -0.062339, -0.129499, -0.058688,  0.265697, -0.003773,  0.009163,  0.423030,  0.049095, -0.017032,
            -0.076027,  0.062962,  0.002006,  0.036832,  0.156818,  0.000100,  0.271630,  0.351707, -0.057629, -0.093839,
             0.028026,  0.018186,  0.214870, -0.020501,  0.067015,  0.344260,  0.200691,  0.091753, -0.010672,  0.115041,
             0.044100,  0.446618, -0.020306,  0.008735,  0.730074,  0.035471, -0.021010, -0.063278,  0.068976,  0.001938,
             0.263608,  0.028412,  0.054971,  0.240539,  0.536313, -0.051542, -0.078298,  0.028946,  0.008783,  0.201529,
            -0.013372,  0.110531,  0.371020,  0.218412, -0.189581, -0.275438,  0.096476,  0.008512,  0.234576, -0.015560,
             0.042574,  0.379747,  0.222364, -0.077602, -0.088737,  0.032865,  0.003561,  0.700925, -0.038900,  0.124275,
             0.374810,  0.232063, -0.074776, -0.078342,  0.032701,  0.001754,  0.218566, -0.010376,  0.054778,  1.117116,
             0.647957, -0.181577,  0.897736,  0.481488,  0.124128,  0.146511, -0.007969,  0.032059,  0.226063,  0.172489,
            -0.057579, -0.099857,  0.034045,  0.004497,  0.417679,  0.213867,  0.065234,  0.403004,  0.223793, -0.081993,
            -0.100358,  0.035702,  0.001779,  0.248812, -0.019477,  0.048196,  0.422473,  0.237281,  0.064291
        };
        for (int k = 0, i = 0; i < (int)(sizeof(in) / sizeof(in[0])); i++) {
            const int n0 = gtoint_basis_shell_get_count(bas[in[i].s[0]]);
            const int n1 = gtoint_basis_shell_get_count(bas[in[i].s[1]]);
            const int n2 = gtoint_basis_shell_get_count(bas[in[i].s[2]]);
            const int n3 = gtoint_basis_shell_get_count(bas[in[i].s[3]]);
            {
                double *const t = (double *)realloc(val, sizeof(double) * n0 * n1 * n2 * n3);
                if (t == NULL) {
                    printf("ERROR: %d: realloc()\n", __LINE__);
                    goto EXIT;
                }
                val = t;
            }
            if ((e = gtoint_compute_weighted_electron_repulsion_integrals(
                itg,
                &(pos[in[i].p[0]]), bas[in[i].s[0]], &(pw[0]), gw[0],
                &(pos[in[i].p[1]]), bas[in[i].s[1]], &(pw[1]), gw[1],
                &(pos[in[i].p[2]]), bas[in[i].s[2]], &(pw[2]), gw[2],
                &(pos[in[i].p[3]]), bas[in[i].s[3]], &(pw[3]), gw[3],
                1, &(in[i].d[0]), &(in[i].d[1]), &(in[i].d[2]), &(in[i].d[3]),
                val
            )) != GTOINT_ERROR_OK) {
                printf("ERROR: %d: gtoint_compute_weighted_electron_repulsion_integrals() -> %d\n", __LINE__, e);
                goto EXIT;
            }
            int l = 0;
            for (int i3 = 0; i3 < n3; i3++) {
            for (int i2 = 0; i2 < n2; i2++) {
            for (int i1 = 0; i1 < n1; i1++) {
            for (int i0 = 0; i0 < n0; i0++) {
                const double u = ref[k++];
                const double v = val[l++];
                if (!(fabs(v - u) <= 1e-2 * fmax(fabs(u), 1.0))) { /* To catch NaN, negated condition is used. */
                    printf(
                        "%3d: {%d[%d]%d[%d]%d[%d]%d[%d]}{%d %d %d %d}{{%d %d %d}{%d %d %d}{%d %d %d}{%d %d %d}}: %9.6f != %9.6f\n",
                        i,
                        in[i].s[0], i0, in[i].s[1], i1, in[i].s[2], i2, in[i].s[3], i3,
                        in[i].p[0], in[i].p[1], in[i].p[2], in[i].p[3],
                        in[i].d[0].x, in[i].d[0].y, in[i].d[0].z, in[i].d[1].x, in[i].d[1].y, in[i].d[1].z,
                        in[i].d[2].x, in[i].d[2].y, in[i].d[2].z, in[i].d[3].x, in[i].d[3].y, in[i].d[3].z,
                        v, u
                    );
                    r = 1;
                }
            }
            }
            }
            }
        }
        printf("Result: %s\n", r ? "FAIL" : "PASS");
    }
EXIT:;
    free(val);
    gtoint_basis_shell_destroy(bas[5]);
    gtoint_basis_shell_destroy(bas[4]);
    gtoint_basis_shell_destroy(bas[3]);
    gtoint_basis_shell_destroy(bas[2]);
    gtoint_basis_shell_destroy(bas[1]);
    gtoint_basis_shell_destroy(bas[0]);
    gtoint_integrator_destroy(itg);
    return (r == 0 && e == GTOINT_ERROR_OK) ? 0 : 1;
}
