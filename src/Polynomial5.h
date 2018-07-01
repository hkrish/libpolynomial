
#ifndef Polynomial5_h__
#define Polynomial5_h__

/*!
 * File: Polynomial5.h
 *
 * Root finding for upto quintic polynomials.
 *
 * Copyright © 2017 Harikrishnan Gopalakrishnan. All rights reserved.
 *
 */

/**
 * Solve a quadratic equation, using numerically stable methods,
 * given an equation of the form `ax^2 + bx + c = 0`.
 *
 * cf. "Solve Quadratic Equation -Floating-Point Computation Without
 * Extra-Precise Arithmetic" - Kahan W
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
int poly_quadroots(double a, double b, double c,
                   double minB, double maxB, double *x1, double *x2);


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
int poly_cubicroots(double ac , double bc , double cc , double dc,
                    double minB, double maxB, double *x1, double *x2, double *x3);


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
int poly_quartroots(double ac, double bc, double cc, double dc, double ec,
                    double minB, double maxB,
                    double *x1, double *x2, double *x3, double *x4);


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
int poly_quintroots(double a0, double b0, double c0, double d0, double e0, double f0,
                    double minB, double maxB,
                    double *x1, double *x2, double *x3, double *x4, double *x5);


#endif  // Polynomial5_h__
