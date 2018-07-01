
/*!
 * File: Polynomial.c
 *
 * Root finding for upto quintic polynomials.
 *
 * Copyright © 2017 Harikrishnan Gopalakrishnan. All rights reserved.
 *
 * Compile with -fvisibility=hidden.
 */


#include "Polynomial5.h"

#define EXPORT __attribute__((visibility("default")))

#define __MATH_LONG_DOUBLE_CONSTANTS


#include <stdio.h>
#include <math.h>
#include <float.h>
#include <simd/simd.h>


// Tolerances for comparing floating point equalities in different precisions
#define EPS_M2            DBL_EPSILON
#define EPS_M2_L          DBL_EPSILON
#define EPS_M2_LD         LDBL_EPSILON
#define EPS_M2_LD2        (LDBL_EPSILON/2.0)
#define EPS_LINEAR        1.0e-12


/**
 * Insertion sort, for sorting roots and root bounds.
 * (number of elements are less that 6)
 */
static inline void isort(double *arr, int n)
{
    for (int i = 1; i < n; i++) {
        double elem = arr[i];
        int j = i-1;
        while (j >= 0 && arr[j] > elem) {
            arr[j+1] = arr[j];
            j = j-1;
        }
        arr[j+1] = elem;
    }
}


/**
 * Stable quadratic discriminant evaluation.
 * cf. "Solve Quadratic Equation -Floating-Point Computation Without
 * Extra-Precise Arithmetic" - Kahan W
 */
static __m128d __attribute__((const))
_DD_split(const double X)
{
    double bigX = X * 134217729;  // = X*(2^27 + 1)
    double Y = (X - bigX);
    double Xh = Y + bigX;         // Don’t optimize  Y  away!
    double Xt = X - Xh;
    return (__m128d) { Xh, Xt };
}

/**
 * Evaluate the discriminant of a quadratic `ax^2 - 2b + c = 0`.
 * D = b^2 - a*c.
 */
EXPORT
double poly_discriminant(double a, double b, double c)
{
    __m128d ad = _DD_split(a);
    __m128d bd = _DD_split(b);
    __m128d cd = _DD_split(c);

    double p = b*b;
    __m128d bd2 = _mm_mul_pd(bd, bd);
    double dp = ((bd2[0] - p) + 2*bd[0]*bd[1]) + bd2[1];

    double q = a*c;
    __m128d adcd = _mm_mul_pd(ad, cd);
    __m128d adxcd = _mm_mul_pd(ad, _mm_shuffle_pd(cd, cd, 1));
    double dq = ((adcd[0] - q) + (adxcd[0] + adxcd[1])) + adcd[1];

    double d = (p-q) + (dp-dq);  // Don’t omit parentheses!
    return d;
}


/* Internal method that accepts long doubles */
static int poly_quadroots_LD(long double a, long double b, long double c,
                             double minB, double maxB, double *x1, double *x2)
{
    *x1 = INFINITY;
    *x2 = INFINITY;
    int rCount = 0;
    long double B = b;
    b *= -0.5;
    // Calculate discriminant with original coefficients
    long double D = b * b - a * c;
    if (fabsl(D) <= 1e-8 || fabsl(D) > 1e8) {
        // Normalise coefficients. Normalisation is done by scaling coefficients
        // with a *power of 2*, so that all the bits in the mantissa remain unchanged.
        long double s = fmaxl(fabsl(a), fmaxl(fabsl(b), fabsl(c)));
        s = (s == 0) ? LDBL_EPSILON : s;
        int sc = 0;
        frexpl(s, &sc);
        sc = 1 - sc;
        a = scalbnl(a, sc);
        b = scalbnl(b, sc);
        B = scalbnl(B, sc);
        c = scalbnl(c, sc);
        // Recalculate discriminant
        D = b * b - a * c;
    }
    // The following test tests the chances of overflow, if yes, evaluate the
    // discriminant using a split-double routine.
    long double E = b * b + a * c;
    if (3.0L * fabsl(D) < E) {
		D = (long double) poly_discriminant((double) a, (double) b, (double) c);
    }

    if (fabsl(a) < EPS_M2) {
        if (fabsl(B) < EPS_M2) {
            return 0;
        }
        *x1 = (double)(-c / B);
    } else {
        // No real roots if D < 0
        if (D >= -EPS_M2) {
            long double Q = D < 0 ? 0 : sqrtl(D);
            long double R = b + (b < 0 ? -1 : 1) * Q;
            if (R == 0.0) {
                *x1 = (double) (c / a);
                *x2 = (double) (-(*x1));
            } else {
                *x1 = (double) (R / a);
                *x2 = (double) (c / R);
            }
        }
    }
    // Set the output values: x1
    if (isfinite(*x1) && ((*x1) >= minB && (*x1) <= maxB)) {
        ++rCount;
    }
    // Set the output values: x2
    if (isfinite(*x2) && ((*x2) >= minB && (*x2) <= maxB) && (rCount == 0 || (*x1) != (*x2))) {
        if (rCount == 0) {
            *x1 = *x2;
        }
        ++rCount;
        if (rCount == 2 && (*x1) > (*x2)) {
            double tmp = *x1;
            *x1 = *x2;
            *x2 = tmp;
        }
    }
    return rCount;
}


/**
 * Solve a quadratic equation, using numerically stable methods,
 * given an equation of the form `ax^2 + bx + c = 0`.
 *
 * - parameters:
 * - ca: Coefficient of x^2
 * - cb: Coefficient of x
 * - cc: Coefficient of 1
 * - minB: The minimum bound of valid roots; `roots < minB` are omited
 * - maxB: The maximum bound of valid roots; `roots > maxB` are omited
 * - x1, x2: Pointers to set the roots of the polynomial.
 * - returns: The number of real roots of the quadratic.
 */
EXPORT
int poly_quadroots(double a, double b, double c, double minB, double maxB, double *x1, double *x2)
{
    return poly_quadroots_LD((long double)a, (long double)b, (long double)c, minB, maxB, x1, x2);
}


/* Internal method that accepts long doubles */
static int poly_cubicroots_LD(long double a, long double b, long double c, long double d,
                              double minB, double maxB, double *x1, double *x2, double *x3)
{
    long double x = INFINITY;
    long double b1;
	long double c2;
    // Normalise coefficients a la *Jenkins & Traub's RPOLY*
    long double s = fmaxl(fabsl(a), fmaxl(fabsl(b), fmaxl(fabsl(c), fabsl(d))));
    if ((s < 1e-7 || s > 1e7) && s != 0) {
        // Scale the coefficients by a multiple of the exponent of ||coeffs.||
        // so that all the bits in the mantissa are preserved.
        int sc = 0;
        frexpl(s, &sc);
        sc = 1 - sc;
        a = scalbnl(a, sc);
        b = scalbnl(b, sc);
        c = scalbnl(c, sc);
        d = scalbnl(d, sc);
    }
    // If a or d is zero, we only need to solve a quadratic, so we set
    // the coefficients appropriately.
    if (a == 0L) {
        a = b;
        b1 = c;
        c2 = d;
    } else if (d == 0L) {
        b1 = b;
        c2 = c;
        x = 0L;
    } else {
        // Iteratively find the leftmost root.
        const long double ec = nextafterl(1.0L, 2.0L);
        x = -(b / a) / 3.0L;
        long double ax = a * x;
        b1 = ax + b;
        c2 = b1 * x + c;
        long double qd = (ax + b1) * x + c2;
        long double q = c2 * x + d;
        long double t = q / a;
        long double r = cbrtl(fabsl(t));
        const long double s = t < 0 ? -1.0L : 1.0L;
        t = -qd / a;
        // λ ~= 1.324717957244746025960909 in λ^3 = λ + 1
        // See Kahan's notes on why 1.324718*... works.
        if (t > 0) {
            r = 1.324717957244746025 * fmaxl(r, sqrtl(t));
        }
        long double x0 = x - s * r;
        if (x0 != x) {
            do {
                x = x0;
                ax = a * x;
                b1 = ax + b;
                c2 = b1 * x + c;
                qd = (ax + b1) * x + c2;
                q = c2 * x + d;
                // Newton's; divide by ec to avoid x0 crossing over a root
                x0 = qd == 0L ? x : x - q / qd / ec;
            } while ((x0 * s) > (x * s));
            // If we have one root converged
            // Caluculate the coefficients for the 'deflated' quadratic
            if (fabsl(a) * x * x > fabsl(d / x)) {
                c2 = -d / x;
                b1 = (c2 - c) / x;
            }
        }
    }
    // Solve the 'deflated' quadratic
    int cntQuad = poly_quadroots_LD(a, b1, c2, minB, maxB, x1, x2);
    double xd = (double) x;
    // Set the output values
    if ((cntQuad > 0 && (*x1) == xd) || (cntQuad > 1 && (*x2) == xd)) {
        xd = INFINITY;
    }
    if (isfinite(xd) && xd >= minB && xd <= maxB) {
        if (cntQuad == 0) {
            *x1 = xd;
        } else if (cntQuad == 1) {
            *x2 = xd;
        } else if (cntQuad == 2) {
            *x3 = xd;
        }
        ++cntQuad;
    }
    return cntQuad;
}


/**
 * Solve a cubic equation, using numerically stable methods,
 * given an equation of the form `ax^3 + bx^2 + cx + d = 0`.
 * This algorithm avoids the problematic trigonometric/inverse
 * trigonometric calculations required by the "Italians'" formula.
 *
 * cf. Kahan W. - "To Solve a Real Cubic Equation"
 * http://www.cs.berkeley.edu/~wkahan/Math128/Cubic.pdf
 * The paper contains inferences on accuracy of cubic
 * zero-finding methods. Also testing methods for robustness.
 *
 * - parameters:
 * - ac: Coefficient of x^3
 * - bc: Coefficient of x^2
 * - cc: Coefficient of x
 * - dc: Coefficient of 1
 * - minB: The minimum bound of valid roots; `roots < minB` are omited
 * - maxB: The maximum bound of valid roots; `roots > maxB` are omited
 * - x1, x2, x3: Pointers to set the roots.
 * - returns: The number of roots of the cubic.
 */
EXPORT
int poly_cubicroots(double ac , double bc , double cc , double dc,
                    double minB, double maxB, double *x1, double *x2, double *x3)
{
    return poly_cubicroots_LD((long double)ac, (long double)bc, (long double)cc, (long double)dc,
                              minB, maxB, x1, x2, x3);
}


/* Internal method that accepts long doubles */
static int poly_quartroots_LD(long double an, long double bn, long double cn, long double dn,
                              long double en, double minB, double maxB,
                              double *x1, double *x2, double *x3, double *x4,
                              int should_polish_roots)
{
    *x1 = INFINITY;
    *x2 = INFINITY;
    *x3 = INFINITY;
    *x4 = INFINITY;
    int rCount = 0;

#define APPEND_ROOT4(x) do {                    \
        if (rCount == 0) { *x1 = (x); }         \
        else if (rCount == 1) { *x2 = (x); }    \
        else if (rCount == 2) { *x3 = (x); }    \
        else { *x4 = (x); }                     \
        ++rCount;                               \
    } while(0)

    // Normalise coefficients a la *Jankins & Traub's RPOLY*
    long double s = fmaxl(fabsl(an),
                          fmaxl(fabsl(bn), fmaxl(fabsl(cn), fmaxl(fabsl(dn), fabsl(en)))));
    if ((s < 1e-7 || s > 1e7) && s != 0) {
        // Scale the coefficients by a multiple of the exponent of ||coeffs.||
        // so that all the bits in the mantissa are preserved
        int sc = 0;
        frexpl(s, &sc);
        sc = 1 - sc;
        an = scalbnl(an, sc);
        bn = scalbnl(bn, sc);
        cn = scalbnl(cn, sc);
        dn = scalbnl(dn, sc);
        en = scalbnl(en, sc);
    }
    if (fabsl(an) < EPS_M2_L) {
        // the polynomial is actually a cubic
        return poly_cubicroots(bn, cn, dn, en, minB, maxB, x1, x2, x3);
    }
    // Original quartic `P(x) = ac x^4 + bc x^3 + cc x^2 + dc x + ec = 0`
    // Convert to monic form `P(x) → x^4 + a x^3 + b x^2 + c x + d = 0`
    long double a = bn/an;
    long double b = cn/an;
    long double c = dn/an;
    long double d = en/an;
    // Powers of a
    long double a2 = a*a;
    long double a3 = a2*a;
    long double a4 = a3*a;
    double a_4d = (double)(a/4.0L);
    // Substitute `x → y-a/4`  =>  `Q(y) → y^4 + e y^2 + f y + g = 0`
    long double e = b - 3.0L*a2/8.0L;
    long double f = c + a3/8.0L - a*b/2.0L;
    long double g = d - 3.0L*a4/256.0L + a2*b/16.0L - a*c/4.0L;
    if (fabsl(g) <= EPS_M2_L) {
        // One root @ 0
        if (-a_4d >= minB && -a_4d <= maxB) {
            APPEND_ROOT4(-a_4d);
        }
        if (fabsl(f) <= EPS_M2_L) {
            // root @ 0
            if (fabsl(e) <= EPS_M2_L) {
                // root @ 0
            } else if (e < 0) {
                // `g~=0, f~=0`  ⇒  solve `y^2 + e = 0`
                double qr1 = (double)sqrtl(-e);
                double r1 = qr1 - a_4d;
                double r2 = -qr1 - a_4d;
                if (r1 >= minB && r1 <= maxB)  {
                    APPEND_ROOT4(r1);
                }
                if (r2 >= minB && r2 <= maxB)  {
                    APPEND_ROOT4(r2);
                }
            }
        } else {
            // `g ~= 0`  ⇒  solve `y^3 + ey + f = 0`
            double qr1 = 0;
            double qr2 = 0;
            double cr3 = 0;
            int ysc = poly_cubicroots_LD(1.0L, 0.0L, e, f, -DBL_MAX, DBL_MAX, &qr1, &qr2, &cr3);
            if (ysc > 0) {
                double r = qr1 - a_4d;
                // Since we already added a root at zero, avoid duplicate roots
                if (r >= minB && r <= maxB && fabs(qr1) > EPS_LINEAR) {
                    APPEND_ROOT4(r);
                }
            }
            if (ysc > 1) {
                double r = qr2 - a_4d;
                if (r >= minB && r <= maxB && fabs(qr2) > EPS_LINEAR) {
                    APPEND_ROOT4(r);
                }
            }
            if (ysc > 2) {
                double r = cr3 - a_4d;
                if (r >= minB && r <= maxB && fabs(cr3) > EPS_LINEAR) {
                    APPEND_ROOT4(r);
                }
            }
        }
    } else if (fabsl(f) <= EPS_M2_L) {
        // `f ~= 0`  ⇒  `y^4 + e y^2 + g`  is a quadratic in y^2
        double qr1 = 0;
        double qr2 = 0;
        int cr = poly_quadroots_LD(1.0L, e, g, -DBL_MAX, DBL_MAX, &qr1, &qr2);
        if (cr > 0 && qr1 >= 0) {
            double r = sqrt(qr1);
            double r1 = r - a_4d;
            double r2 = -r - a_4d;
            if (r1 >= minB && r1 <= maxB) {
                APPEND_ROOT4(r1);
            }
            if (r2 >= minB && r2 <= maxB) {
                APPEND_ROOT4(r2);
            }
        }
        if (cr > 1 && qr2 >= 0) {
            double r = sqrt(qr2);
            double r3 = r - a_4d;
            double r4 = -r - a_4d;
            if (r3 >= minB && r3 <= maxB) {
                APPEND_ROOT4(r3);
            }
            if (r4 >= minB && r4 <= maxB) {
                APPEND_ROOT4(r4);
            }
        }
    } else {
        // Coefficients of cubic in h^2
        long double b1 = 2.0L*e;
        long double c1 = e*e - 4.0L*g;
        long double d1 = -f*f;
        double qr1 = 0;
        double qr2 = 0;
        double cr3 = 0;
        int ysc = poly_cubicroots_LD(1.0L, b1, c1, d1, -DBL_MAX, DBL_MAX, &qr1, &qr2, &cr3);
        if (ysc == 0) {
            return rCount;
        }
        double h2 = qr1;
        if (ysc > 1) {
            h2 = fmax(h2, qr2);
        }
        if (ysc > 2) {
            h2 = fmax(h2, cr3);
        }
        long double h = sqrtl((long double)h2);
        long double j = (e + ((long double)h2) - f/h)/2.0L;
        qr1 = 0;
        qr2 = 0;
        double r1 = INFINITY;
        double r2 = INFINITY;
        double r3 = INFINITY;
        int cr = poly_quadroots_LD(1.0L, h, j, -DBL_MAX, DBL_MAX, &qr1, &qr2);
        if (cr > 0) {
            r1 = qr1 - a_4d;
            if (r1 >= minB && r1 <= maxB) {
                APPEND_ROOT4(r1);
            }
        }
        if (cr > 1) {
            r2 = qr2 - a_4d;
            if (r2 >= minB && r2 <= maxB && (fabs(r2 - r1) > EPS_M2)) {
                APPEND_ROOT4(r2);
            }
        }
        qr1 = 0;
        qr2 = 0;
        cr = poly_quadroots_LD(1.0L, -h, g/j, -DBL_MAX, DBL_MAX, &qr1, &qr2);
        if (cr > 0) {
            r3 = qr1 - a_4d;
            if (r3 >= minB && r3 <= maxB &&
                ((fabs(r3 - r1) > EPS_M2) && (fabs(r3 - r2) > EPS_M2))) {
                APPEND_ROOT4(r3);
            }
        }
        if (cr > 1) {
            double r4 = qr2 - a_4d;
            if (r4 >= minB && r4 <= maxB &&
                ((fabs(r4 - r1) > EPS_M2) && (fabs(r4 - r2) > EPS_M2) &&
                 (fabs(r4 - r3) > EPS_M2))) {
                APPEND_ROOT4(r4);
            }
        }
    }
    if (should_polish_roots && rCount > 0 && rCount <= 4) {
        // Evaluate the original polynomial and try to improve accuracy with a few
        // rounds of Newton-Raphson iterations
        int i;
        long double roots[4] = {*x1, *x2, *x3, *x4};
        long double zeros[4] = {0, 0, 0, 0};
        char should_polish_any_roots = 0;
        for (i = 0; i < rCount; ++i) {
            long double t = roots[i];
            zeros[i] = (((an * t + bn) * t + cn) * t + dn) * t + en;
            should_polish_any_roots = fabsl(zeros[i]) > EPS_M2_L;
        }
        if (should_polish_any_roots) {
            long double dan = an*4;
            long double dbn = bn*3;
            long double dcn = cn*2;
            long double ddn = dn;
            for (i = 0; i < rCount; ++i) {
                long double t = roots[i];
                long double z = zeros[i];
                long double dz = ((dan * t + dbn) * t + dcn) * t + ddn;
                if (dz == 0L) { continue; }
                t -= z / dz;
                z = (((an * t + bn) * t + cn) * t + dn) * t + en;
                if (fabsl(z) < fabsl(zeros[i]) && fabsl(z) > EPS_M2_L) {
                    roots[i] = t;
                    // 2nd iteration
                    zeros[i] = z;
                    dz = ((dan * t + dbn) * t + dcn) * t + ddn;
                    if (dz == 0L) { continue; }
                    t -= z / dz;
                    z = (((an * t + bn) * t + cn) * t + dn) * t + en;
                    if (fabsl(z) < fabsl(zeros[i]) && fabsl(z) > EPS_M2_L) {
                        roots[i] = t;
                        // 3rd iteration
                        zeros[i] = z;
                        dz = ((dan * t + dbn) * t + dcn) * t + ddn;
                        if (dz == 0L) { continue; }
                        t -= z / dz;
                        z = (((an * t + bn) * t + cn) * t + dn) * t + en;
                        if (fabsl(z) < fabsl(zeros[i])) {
                            roots[i] = t;
                        }
                    }
                }
            }
            *x1 = roots[0];
            *x2 = roots[1];
            *x3 = roots[2];
            *x4 = roots[3];
        }
    }
    return rCount;
#undef APPEND_ROOT4
}


/**
 * Solve a quartic equation, using numerically stable methods,
 * given an equation of the form `ax^4 + bx^3 + cx^2 + dx + e = 0`.
 * This algorithm uses the robust cubic and quartic solvers to avoid
 * roundoff errors as much as possible.
 *
 * - parameters:
 * - ac: Coefficient of x^4
 * - bc: Coefficient of x^3
 * - cc: Coefficient of x^2
 * - dc: Coefficient of x
 * - ec: Coefficient of 1
 * - minB: The minimum bound of valid roots; `roots < minB` are omited
 * - maxB: The maximum bound of valid roots; `roots > maxB` are omited
 * - x1, x2, x3, x4: Pointers to set the roots.
 * - returns: The number of real roots of the quartic.
 */
EXPORT
int poly_quartroots(double an, double bn, double cn, double dn, double en,
                    double minB, double maxB,
                    double *x1, double *x2, double *x3, double *x4)
{
    return poly_quartroots_LD((long double)an, (long double)bn, (long double)cn, (long double)dn,
                              (long double)en,  minB,  maxB, x1,  x2,  x3,  x4, 1);
}


/**
 *  Refines a root (root search) in an interval known to contain exactly one root
 *  of a quintic polynomial.
 *
 *  cf. "Bisected Direct Quadratic Regula Falsi"
 *  - Robert G. Gottlieb and Blair F. Thompson
 */
static double refine_root(long double a, long double b, long double c, long double d,
                          long double e, long double f, double minB, double maxB)
{
    long double x1 = (long double)minB;
    long double x2 = (long double)maxB;
    long double y = ((((a * x1 + b) * x1 + c) * x1 + d) * x1 + e) * x1 + f;
    long double y2 = ((((a * x2 + b) * x2 + c) * x2 + d) * x2 + e) * x2 + f;
    long double tol = EPS_M2_LD2 * 2;
    long double xLast = LDBL_MAX;
    long double xup = x1;
    long double xdn = x2;
    long double yup = y;
    long double ydn = y2;
    long double delta, xm, ym, a_d, a_, b_, xden, x;

    if (y < 0) {
        xup = x2;
        xdn = x1;
        yup = y2;
        ydn = y;
    }
    while (fabsl(y) > tol) {
        delta = (xup - xdn) / 2.0;
        xm = (xup + xdn) / 2.0;
        ym = ((((a * xm + b) * xm + c) * xm + d) * xm + e) * xm + f;
        a_d = (2*delta*delta);
        a_ = (yup + ydn - 2*ym) / a_d;
        b_ = (yup - ydn) / (2*delta);
        xden = sqrtl(1 - 4*a_*ym/(b_*b_));
        x = xm - 2*ym / (b_ * (1 + xden));
        if (x == xLast) {
            break;
        }
        xLast = x;
        y = ((((a * x + b) * x + c) * x + d) * x + e) * x + f;
        if (y > 0) {
            yup = y;
            xup = x;
            if (ym < 0) {
                ydn = ym;
                xdn = xm;
            }
        } else {
            ydn = y;
            xdn = x;
            if (ym > 0) {
                yup = ym;
                xup = xm;
            }
        }
    }
    return (double)xLast;
}


/**
 * Evaluate a quintic polynomial with error bounds.
 * The quintic is given as `f(x) = ax^5 + bx^4 + cx^3 + dx^2 + ex + f = 0`
 *
 * - parameters:
 * - a..f: Coefficients of the quintic
 * - at: The polynomial will be evaluated at x = at
 * - error: Pointer to set the error bound.
 * - returns: The value of the polynomial at 'x = at'.
 */
static long double eval_quint(long double a, long double b, long double c, long double d,
                              long double e, long double f, double at, long double *error)
{

    long double t, val, e1, e2, e3, e4;
    t = (long double)at;
    val = ((((a * t + b) * t + c) * t + d) * t + e) * t + f;
    e1 = (fabsl(b) + t*fabsl(a))*EPS_M2_LD2;
    e2 = (fabsl(c) + t*fabsl(b))*EPS_M2_LD2 + t * e1;
    e3 = (fabsl(d) + t*fabsl(c))*EPS_M2_LD2 + t * e2;
    e4 = (fabsl(e) + t*fabsl(d))*EPS_M2_LD2 + t * e3;
    *error = (fabsl(f) + t*fabsl(e))*EPS_M2_LD2 + t * e4;
    return val;
}


/**
 * Solve a quintic equation, using numerically stable methods,
 * given an equation of the form `a x^5 + b x^4 + c x^3 + d x^2 + e x + f = 0`.
 * This algorithm uses the robust cubic and quartic solvers to avoid
 * roundoff errors as much as possible.
 *
 * - parameters:
 * - a0: Coefficient of x^5
 * - b0: Coefficient of x^4
 * - c0: Coefficient of x^3
 * - d0: Coefficient of x^2
 * - e0: Coefficient of x
 * - f0: Coefficient of 1
 * - minB: The minimum bound of valid roots; `roots < minB` are omited
 * - maxB: The maximum bound of valid roots; `roots > maxB` are omited
 * - x1..x5: Pointers to set the roots of the polynomial.
 * - returns: The number of real roots of the quintic.
 */
EXPORT
int poly_quintroots(double a0, double b0, double c0, double d0, double e0, double f0,
                    double minB, double maxB,
                    double *x1, double *x2, double *x3, double *x4, double *x5)
{
    int rCount = 0;

#define APPEND_ROOT5(x) do {                    \
        if (rCount == 0) { *x1 = (x); }         \
        else if (rCount == 1) { *x2 = (x); }    \
        else if (rCount == 2) { *x3 = (x); }    \
        else if (rCount == 3) { *x4 = (x); }    \
        else { *x5 = (x); }                     \
        ++rCount;                               \
    } while(0)

    long double a0l = (long double)a0;
    long double b0l = (long double)b0;
    long double c0l = (long double)c0;
    long double d0l = (long double)d0;
    long double e0l = (long double)e0;
    long double f0l = (long double)f0;
    // Normalise coefficients a la *Jankins & Traub's RPOLY*
    long double s = fmaxl(fabsl(a0l),
                          fmaxl(fabsl(b0l),
                                fmaxl(fabsl(c0l),
                                      fmaxl(fabsl(d0l), fmaxl(fabsl(e0l), fabsl(f0l))))));
    if ((s < 1e-7 || s > 1e7) && s != 0) {
        // Scale the coefficients by a multiple of the exponent of ||coeffs.||
        // so that all the bits in the mantissa are preserved
        int sc = 0;
        frexpl(s, &sc);
        sc = 1 - sc;
        a0l = scalbnl(a0l, sc);
        b0l = scalbnl(b0l, sc);
        c0l = scalbnl(c0l, sc);
        d0l = scalbnl(d0l, sc);
        e0l = scalbnl(e0l, sc);
        f0l = scalbnl(f0l, sc);
    }
    // Check for obvious cases, f0 == 0, a0 = 0, etc.
    if (fabs(a0) < EPS_M2) {
        return poly_quartroots_LD(b0l, c0l, d0l, e0l, f0l, minB, maxB, x1, x2, x3, x4, 1);
    }
    // Solve the quintic
    // Convert to monic form `P(x) → x^5 + ac x^4 + bc x^3 + cc x^2 + dc x + ec = 0`
    long double ac = b0l / a0l;
    long double bc = c0l / a0l;
    long double cc = d0l / a0l;
    long double dc = e0l / a0l;
    long double ec = f0l / a0l;
    // Convert to a *depressed* quintic
    // `x → y-ac/5` ⇒ `Q(y) → y^5 + a y^3 + b y^2 + c y + d`
    double a_5d = (double)(ac / 5.0);
    long double ac2 = ac * ac;
    long double ac3 = ac2 * ac;
    long double cc25 = 25 * cc;
    long double acbc = ac * bc;
    long double c1 = 3*ac3 - 15*acbc + 2*cc25;
    long double a = -ac2 * 0.4 + bc;               // y^3
    long double b = (c1 + ac3 - cc25) * 0.04;      // y^2
    long double c = dc - c1*ac * 0.008;            // y^1
    long double d = ac*(ac*(4*ac3 - 25*acbc + 5*cc25) - 625*dc) * 0.00032 + ec;  // y^0
    if (fabsl(d) <= EPS_M2_LD2) {
        // The depressed quintic is actually a quartic polynomial
        double nMinB = minB + a_5d;
        double nMaxB = maxB + a_5d;
        rCount = poly_quartroots_LD(1.0L, 0.0L, a, b, c, nMinB, nMaxB,
                                    x1, x2, x3, x4, 1);
        *x1 -= a_5d;
        *x2 -= a_5d;
        *x3 -= a_5d;
        *x4 -= a_5d;
        if (-a_5d >= minB && -a_5d <= maxB) {
            APPEND_ROOT5(-a_5d);
        }
        return rCount;
    }
    // Calculate the *isolator polynomials* (a pair of quadratics)
    // cf. "Isolator Polynomials" - Sederberg T., Chang G.
    long double a2 = a * a;
    long double b2 = b * b;
    // First isolator polynomial A(x)
    long double siAx2 = 12*a2*a + 45*b2 - 40*a*c;
    long double siAx = 8*a2*b + 60*b*c - 50*a*d;
    long double siA1 = 4*a2*c + 75*b*d;
    // second isolator polynomial B(x)
    long double siBx2 = 10*a;
    long double siBx = -15*b;
    long double siB1 = 4*a2;
    double qr1 = 0;
    double qr2 = 0;
    double qr3 = 0;
    double qr4 = 0;
    int siArc = poly_quadroots_LD(siAx2, siAx, siA1, -DBL_MAX, DBL_MAX, &qr1, &qr2);
    int siBrc = poly_quadroots_LD(siBx2, siBx, siB1, -DBL_MAX, DBL_MAX, &qr3, &qr4);
    // Calculcate root bounds of the depressed poly (and original quintic by adding |ac/5|)
    // Bounds due to Lagurre, (cf. Samuelson's inequality)
    // Our roots have an arithmetic mean of 0, since y^4's coefficient = 0
    double rBound = ceil(sqrt((double)(10.0/4.0 * fabsl(a))) * 4.0/5.0 + fabs(a_5d));
    // Sort the bounds of the roots
    double rBounds[6] = {-rBound, rBound, rBound, rBound, rBound, rBound};
    int rbCount = 1;
    if (siArc > 0) {
        rBounds[rbCount++] = qr1-a_5d;
    }
    if (siBrc > 0) {
        rBounds[rbCount++] = qr3-a_5d;
    }
    if (siArc > 1) {
        rBounds[rbCount++] = qr2-a_5d;
    }
    if (siBrc > 1) {
        rBounds[rbCount++] = qr4-a_5d;
    }
    isort(rBounds, rbCount);
    rBounds[rbCount++] = rBound;
    // Find a root interval from the roots of isolator polynomials for
    // which the polynomial changes sign at the limits.
    double y0 = rBounds[0];
    double y1 = rBounds[1];
    long double ey0, ey1;
    long double fy0 = eval_quint(a0l, b0l, c0l, d0l, e0l, f0l, y0, &ey0);
    long double fy1 = eval_quint(a0l, b0l, c0l, d0l, e0l, f0l, y1, &ey1);
    int isBoundZero = (fabsl(fy0) <= 4*ey0) || (fabsl(fy1) <= 4*ey1) ? 1 : 0;
    int nr = 2;
    while (nr < rbCount && (fy0 * fy1 > 0 && !isBoundZero)) {
        y0 = y1;
        fy0 = fy1;
        ey0 = ey1;
        y1 = rBounds[nr];
        fy1 = eval_quint(a0l, b0l, c0l, d0l, e0l, f0l, y1, &ey1);
        isBoundZero = isBoundZero || (fabsl(fy1) <= 4*ey1) ? 1 : 0;
        ++nr;
    }
    if (fy0 * fy1 > 0 && !isBoundZero) {
        // CAUTION: It is possible that we will omit a root of
        // even multiplicity due to roundoff.
        return 0;
    }
    // Find the leftmost root of the original polynomial
    double t0d = refine_root(a0l, b0l, c0l, d0l, e0l, f0l, y0, y1);
    // *Polish* the root using a few iterations of Newton-Raphson method
    long double t0 = (long double)t0d;
    long double a0ld = a0l*5;
    long double b0ld = b0l*4;
    long double c0ld = c0l*3;
    long double d0ld = d0l*2;
    long double e0ld = e0l;
    long double fr = ((((a0l * t0 + b0l) * t0 + c0l) * t0 + d0l) * t0 + e0l) * t0 + f0l;
    long double f2r = (((a0ld * t0 + b0ld) * t0 + c0ld) * t0 + d0ld) * t0 + e0ld;
    long double t = t0 - fr/f2r;
    long double fr2 = ((((a0l * t + b0l) * t + c0l) * t + d0l) * t + e0l) * t + f0l;
    if (fabsl(fr2) < fabsl(fr)) {
        t0 = t;
        fr = ((((a0l * t0 + b0l) * t0 + c0l) * t0 + d0l) * t0 + e0l) * t0 + f0l;
        f2r = (((a0ld * t0 + b0ld) * t0 + c0ld) * t0 + d0ld) * t0 + e0ld;
        if (f2r != 0) {
            t = t0 - fr/f2r;
            fr2 = ((((a0l * t + b0l) * t + c0l) * t + d0l) * t + e0l) * t + f0l;
            if (fabsl(fr2) < fabsl(fr)) {
                t0 = t;
                fr = ((((a0l * t0 + b0l) * t0 + c0l) * t0 + d0l) * t0 + e0l) * t0 + f0l;
                f2r = (((a0ld * t0 + b0ld) * t0 + c0ld) * t0 + d0ld) * t0 + e0ld;
                if (f2r != 0) {
                    t = t0 - fr/f2r;
                    fr2 = ((((a0l * t + b0l) * t + c0l) * t + d0l) * t + e0l) * t + f0l;
                    if (fabsl(fr2) < fabsl(fr)) {
                        t0 = t;
                    }
                }
            }
        }
    }
    t0d = (double)t0;
    int isT0Root = t0d >= minB && t0d <= maxB ? 1 : 0;
    if (isT0Root) {
        APPEND_ROOT5(t0d);
    }
    // Use the found root to deflate the quintic into a quartic and solve it
    long double bq = b0l + a0l*t0;
    long double cq = c0l + bq*t0;
    long double dq = d0l + cq*t0;
    long double eq = e0l + dq*t0;
    int qrs = 0;
    if (isT0Root) {
        qrs = poly_quartroots_LD(a0l, bq, cq, dq, eq, minB, maxB, x2, x3, x4, x5, 0);
    } else {
        qrs = poly_quartroots_LD(a0l, bq, cq, dq, eq, minB, maxB, x1, x2, x3, x4, 0);
    }
    rCount = rCount + qrs;
    if (rCount > 0 && rCount <= 5) {
        // Evaluate the original polynomial and try to improve accuracy with a few
        // rounds of Newton-Raphson iterations
        int i;
        long double roots[5] = {*x1, *x2, *x3, *x4, *x5};
        long double zeros[5] = {0, 0, 0, 0, 0};
        char should_polish_any_roots = 0;
        for (i = 0; i < rCount; ++i) {
            long double t = roots[i];
            zeros[i] = ((((a0l * t + b0l) * t + c0l) * t + d0l) * t + e0l) * t + f0l;
            should_polish_any_roots = fabsl(zeros[i]) > EPS_M2_L;
        }
        if (should_polish_any_roots) {
            long double dan = a0l*5;
            long double dbn = b0l*4;
            long double dcn = c0l*3;
            long double ddn = d0l*2;
            long double den = e0l;
            for (i = 0; i < rCount; ++i) {
                long double t = roots[i];
                long double z = zeros[i];
                long double dz = (((dan * t + dbn) * t + dcn) * t + ddn) * t + den;
                if (dz == 0L) { continue; }
                t -= z / dz;
                z = ((((a0l * t + b0l) * t + c0l) * t + d0l) * t + e0l) * t + f0l;
                if (fabsl(z) < fabsl(zeros[i]) && fabsl(z) > EPS_M2_L) {
                    roots[i] = t;
                    // 2nd iteration
                    zeros[i] = z;
                    dz = (((dan * t + dbn) * t + dcn) * t + ddn) * t + den;
                    if (dz == 0L) { continue; }
                    t -= z / dz;
                    z = ((((a0l * t + b0l) * t + c0l) * t + d0l) * t + e0l) * t + f0l;
                    if (fabsl(z) < fabsl(zeros[i]) && fabsl(z) > EPS_M2_L) {
                        roots[i] = t;
                        // 3rd iteration
                        zeros[i] = z;
                        dz = (((dan * t + dbn) * t + dcn) * t + ddn) * t + den;
                        if (dz == 0L) { continue; }
                        t -= z / dz;
                        z = ((((a0l * t + b0l) * t + c0l) * t + d0l) * t + e0l) * t + f0l;
                        if (fabsl(z) < fabsl(zeros[i])) {
                            roots[i] = t;
                        }
                    }
                }
            }
            *x1 = roots[0];
            *x2 = roots[1];
            *x3 = roots[2];
            *x4 = roots[3];
            *x5 = roots[4];
        }
    }
    return rCount;
#undef APPEND_ROOT5
}
