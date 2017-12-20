
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <mach/mach_time.h>

#include "Polynomial.h"



// ==================================================
// Utilities
//

// Tolerance for comparing the results for equality
#define ACCURACY 1e-10


/**
 * Insertion sort, for sorting roots and root bounds.
 * (number of elements are less that 6)
 */
static void isort(double *arr, int n)
{
    for (int i = 1; i < n; i++) {
        double key = arr[i];
        int j = i-1;
        while (j >= 0 && arr[j] > key) {
            arr[j+1] = arr[j];
            j = j-1;
        }
        arr[j+1] = key;
    }
}


static uint64_t timestamp()
{
    return mach_absolute_time();
}


static uint64_t get_elapsed(const uint64_t t0, const uint64_t t1)
{
    mach_timebase_info_data_t timebase_info = { .denom = 0 };
    mach_timebase_info(&timebase_info);

    return ((t1 - t0) * timebase_info.numer) / timebase_info.denom;
}


static double get_elapsed_seconds(const uint64_t t0, const uint64_t t1)
{
    double nanos = (double) get_elapsed(t0, t1);
    return nanos * 1.0e-9;
}


void print_result(long number, char* name, double t)
{
    long n = number;
    long r = 0;
    char *short_scale_names[] = {"thousand", "million", "billion"};
    int scale_index = 0;
    char str[4096] = {0};

    sprintf(str, "%ld", n);

    while (n >= 1000 && scale_index < 3) {
        long d = n / 1000;
        r = n - (d * 1000);
        n = d;
        if (scale_index > 0 || (scale_index == 0 && n >= 10)) {
            if (r > 0) {
                sprintf(str, "%ld (%-.2lf %s)", number, (float)n + (float)r/1000.0, short_scale_names[scale_index]);
            } else {
                sprintf(str, "%ld (%ld %s)", number, n, short_scale_names[scale_index]);
            }
        }
        ++scale_index;
    }

    printf("Solved %s %s polynomials in %-.3f seconds.\n", str, name, t);
}


int testArrayApprox(double* a, double*b, int n)
{
    if (n == 0) {
        return 1;
    }
    double err = 0, errAvg = 0;
    for (int i=0; i < n; ++i) {
        double diff = fabs(a[i] - b[i]);
        err = fmax(err, diff);
        errAvg += diff;
    }
    /* printf("%d, %25.17lf, %25.17lf\n", n, err, errAvg/(double)n); */
    return err <= ACCURACY && errAvg/(double)n <= ACCURACY ? 1 : 0;
}


// ==================================================
// Tests
//
int testCubic()
{
    // Test data
    int nTests = 13;
    double tests[][4] = {
        {0, 0, 0, 0},
        {-1, 0, 0, -1},
        {1, 3, 3, 1},
        {1, 0, -2, 5},
        {1, -3, 2, 0},
        {1, -3, 2, 2.34e-89},
        {1, -7999999999, -8000000002, 16000000000},
        {16000000000, -8000000002, -7999999999, 1},
        {1, -99999.00001, -99999.00001, 1},
        {0.01, -300, 2990000, -299},
        {-3, 3, -1, 0.1111111111},
        {10000000000, -9999999998, -10000000000, 9999999998},
        {9999999998, -10000000000, -9999999998, 10000000000}
    };
    // Answers
    double roots[][3] = {
        {},
        {-1.0},
        {-1.0},
        {-2.09455148154233},
        {2.0, 1.0, 0.0},
        {2.0, 1.0, -1.17e-89},
        {-2.0, 1.0, 8000000000.0},
        {-0.5, 1.25e-10, 1.0},
        {-1.0, 1e-05, 100000.0},
        {0.000100000001003344},
        {0.333178613737119},
        {-1.0, 1.0, 0.9999999998},
        {-1, 1.0000000002, 1}
    };
    int nRoots[13] = {0, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 3, 3};
    // Run the tests
    int nFailed = 0;
    int nMissed = 0;
    printf("---------------------------------------------\n");
    printf("[Cubic polynomial roots]: Running %d tests \n", nTests);
    for (int i=0; i<nTests; ++i) {
        double xs[3] = {0};
        double *cc = tests[i];
        int rc = poly_cubicroots(cc[0], cc[1], cc[2], cc[3], -DBL_MAX, DBL_MAX,
                                 &xs[0], &xs[1], &xs[2]);
        isort(xs, rc);
        if (rc == nRoots[i]) {
            double *rs = roots[i];
            isort(rs, nRoots[i]);
            if (testArrayApprox(rs, xs, rc) == 0) {
                ++nFailed;
                printf("%d) %d / %d: %24.16lf, %24.16lf, %24.16lf \n", i, nRoots[i], rc, fabs(xs[0]-roots[i][0]), fabs(xs[1]-roots[i][1]), fabs(xs[2]-roots[i][2]));
            }
        } else {
            ++nMissed;
            ++nFailed;
            printf("%d) %d / %d: %24.16lf, %24.16lf, %24.16lf \n", i, nRoots[i], rc, fabs(xs[0]-roots[i][0]), fabs(xs[1]-roots[i][1]), fabs(xs[2]-roots[i][2]));
        }
    }
    if (nMissed==0 && nFailed==0){
        printf("[Cubic polynomial roots]: PASSED.\n");
    } else {
        printf("[Cubic polynomial roots]: Missed roots in %d cases. Failed %d tests in total.\n",
               nMissed, nFailed);
    }
    return nFailed;
}

int testQuartic()
{
    // Test data
    int nTests = 18;
    double tests[][5] = {
        {1.0,-2.63453369092816, 11.47594580628995,-3.01193978276912,0.17801932206896},
        {1.0,-2.63453369092816,2.47594580628995,-0.969679092568218,0.133302415104222},
        {16.0,-16.62910939950919,8.30968497920564,-1.86739743049595,0.12844638447658},
        {1.0, 6, -5, -10, -3},
        {1.0, -2.35220164231774431, 2.00357193581024839, -0.732009944503505094, 0.0963567270541595411},
        {4213952.415, -6257401.32508918, 2907561.77975298, -488990.675178038, 20927.9939160697},
        {57341.571, -161727.801652109, 168247.656267706, -76482.4511413003, 12801.3771995616},
        {1.624401562e-21, -4.54214506896612e-21, 4.72995343184905e-21, -2.17240456536774e-21, 3.70944181855249e-22},
        {1.324862655e-21, -2.25090952994387e-21, 1.2704566968795e-21, -2.55885000154948e-22, 9.1991026101414e-24},
        {3852.061343, -8484.23316822249, 6538.97688624379, -2056.39000615966, 212.050414780893},
        {3340.933089, -4990.31721011812, 3726.99260868093, -1391.74258071197, 194.890453343193},
        {3294.612452, -5742.33193538764, 3753.21565167956, -1090.27543254508, 118.768207298396},
        {1.63064755140304e-16, -2034.51169316325, 1215.54655528735, -322.775442718126, 32.1410930237635},
        {1606.268505, -3167.86593395418, 2156.86431456392, -586.683298630047, 53.9570421132115},
        {3063.145275, -4196.79839872473, 2875.00513659746, -290.139062231677, 5.32413681477471},
        {1836.085081, -3157.51102285065, 2036.23649359839, -378.014005550609, 27.7813065153978},
        {1347.375714, -2372.21108226666, 1566.20719086287, -161.089864668723, 5.63065402793846},
        {463.822363, -787.053152298914, 500.826971130479, -141.64095202131, 1.79148644473404}
    };
    // Answers
    double roots[][4] = {
        {0.08809808973990371, 0.1869183562478567},
        {0.332355001781926, 0.478591097008222},
        {0.117153125076721, 0.308784762589995},
        {-6.54138126514911, -0.618033988749895, -0.458618734850890, 1.61803398874989},
        {0.339953359441404, 0.553720697751012, 0.588050410578735, 0.870477174546594},
        {0.0638539873712612, 0.267956069261243, 0.371231133437524, 0.781883343680119},
        {0.489776436668900, 0.693817355733370, 0.705107127481424, 0.931727590038177},
        {0.5399708576228994, 0.6990489875211234, 0.7435930338679849, 0.8135830710668189},
        {0.0453709962417769, 0.401238808036131, 0.522771212626904, 0.729594738392533},
        {0.2018746976609953, 0.5506293634199524, 0.550629470831883, 0.8993841366181012},
        {0.3734224150244664},
        {0.4357365258470496, 0.4357432884863267},
        {0.1991545128284849},
        {0.202385755872277, 0.316099035619187, 0.669995731256118, 0.783709011002991},
        {0.0237485603932806, 0.0925818978834569},
        {},
        {},
        {0.01325662945878681, 0.8351858078659891}
    };
    int nRoots[18] = {2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 1, 2, 1, 4, 2, 0, 0, 2};
    // Run the tests
    int nFailed = 0;
    int nMissed = 0;
    printf("---------------------------------------------\n");
    printf("[Quartic polynomial roots]: Running %d tests \n", nTests);
    for (int i=0; i<nTests; ++i) {
        double xs[4] = {0};
        double *cc = tests[i];
        int rc = poly_quartroots(cc[0], cc[1], cc[2], cc[3], cc[4], -DBL_MAX, DBL_MAX,
                                 &xs[0] ,&xs[1] ,&xs[2] ,&xs[3]);
        isort(xs, rc);
        if (rc == nRoots[i]) {
            double *rs = roots[i];
            isort(rs, nRoots[i]);
            if (testArrayApprox(rs, xs, rc) == 0) {
                ++nFailed;
                printf("%d) %d / %d: %24.16lf, %24.16lf, %24.16lf, %24.16lf \n", i, nRoots[i], rc, fabs(xs[0]-roots[i][0]), fabs(xs[1]-roots[i][1]), fabs(xs[2]-roots[i][2]), fabs(xs[3]-roots[i][3]));
            }
        } else {
            ++nMissed;
            ++nFailed;
            printf("%d) %d / %d: %24.16lf, %24.16lf, %24.16lf, %24.16lf \n", i, nRoots[i], rc, fabs(xs[0]-roots[i][0]), fabs(xs[1]-roots[i][1]), fabs(xs[2]-roots[i][2]), fabs(xs[3]-roots[i][3]));
        }
    }
    if (nMissed==0 && nFailed==0){
        printf("[Quartic polynomial roots]: PASSED.\n");
    } else {
        printf("[Quartic polynomial roots]: Missed roots in %d cases. Failed %d tests in total.\n",
               nMissed, nFailed);
    }
    return nFailed;
}


int testQuintic()
{
    // Test data
    int nTests = 17;
    double tests[][6] = {
        {1.0, -2.52910939950919052, 2.30968497920563746, -0.896739743049594606, 0.128446384476579654, -0.00275665806066778772},
        {1.0, -2.52910939950919052, 2.30968497920563746, -0.937101799553697019, 0.169404234782936457, -0.0110787163926311995},
        {1.0, -2.63453369092816168, 2.47594580628994667, -0.969679092568217881, 0.133302415104222411, -0.00159027602636815394},
        {1.0, -2.63453369092816168, 2.47594580628994667, -1.01193978276912147, 0.178019322068959658, -0.0112575475542043011},
        // The quintic from the paper "Isolator Polynomials -- Segerberg et. al."
        {1.0, 0, -23, 6, 112, -96},
        {1, -10, 15, 4, -16, 400},
        {1, -3.5, 4.875, -3.375, 1.16015625, -0.158203125},
        {1, 3, -5, -15, 4, 12},
        {1, -1.5, -15.29, 23.895, -11.36, 1.68},
        {1 ,-0.0037841796875, 4.61935997009277e-6, -2.25554686039686e-9, 4.40536496171262e-13, -2.77555756156289e-17},
        {1, -1.921875, 1.02197265625, -0.095947265625, 0.00285911560058594, -2.50339508056641e-5},
        {1.0, -2.52910939950919052, 2.30968497920563746, -0.896739743049594606, 0.128446384476579654, -0.00095665806066778772},
        {1.0, -2.5291093995091905, 2.3096849792056375, -0.89673974304959461, 0.12844638447657965, -0.0069566580606677877},
        {1, 1.470890600490809, 0.6165349007762284, 0.05472411292864265, -0.003358094522770555, -0.003797486716822673},
        {796.692571, -9.27278946051005, 183.852120446035, -936.213438348927, 1630927.19368701, -2068.9112857205},
        {1920000, -4800000, 4400802, -1801203, 307951.125, -13775.0625},
        {22228992000000.0, -55572480000000.0, 57117707539200.0, -30104081308800.0, 7961050818000.0, -815594524200.0}
    };
    // Answers
    double roots[][5] = {
        {0.0258125965473582025, 0.27031189830165379, 0.566533726679526999},
        {0.167473494200687896, 0.225008383579697391, 0.39285214755595044, 0.763169784447613723, 0.980605589725241069},
        {0.0131452033611805842, 0.265990492960650847, 0.573523846559096417, 0.863687532492219689, 0.91818661555501415},
        {0.451580416836080467, 0.793468676205692865, 1.04827269910349},
        {-4.0, -3, 1, 2, 4},
        {-2.2123376379624875742, 3.6437879671398742453, 7.9945186684323459459},
        {0.5, 0.75},
        {-3, -2, -1, 1, 2},
        {-4, 0.3, 0.5, 0.7, 4},
        {0.000122070052903539, 0.000244141111744454, 0.000488280966065109, 0.000976562560843103, 0.0019531249959438}, // *Actual* roots {0.0001220703125, 0.000244140625, 0.00048828125, 0.0009765625, 0.001953125},
        {0.015625, 0.03125, 0.0625, 0.875, 0.9375},
        {0.007871831091751354},
        {0.9521526927712383},
        {0.15215269277155521},
        {0.001268550061539374},
        {0.0669872981077807, 0.3571167434581644, 0.5, 0.6428832565418356, 0.9330127018922193},
        {0.2697317668590528, 0.5, 0.7302682331409472}
    };
    int nRoots[17] = {3, 5, 5, 3, 5, 3, 2, 5, 5, 5, 5, 1, 1, 1, 1, 5, 3};
    // Run the tests
    int nFailed = 0;
    int nMissed = 0;
    printf("---------------------------------------------\n");
    printf("[Quintic polynomial roots]: Running %d tests \n", nTests);
    for (int i=0; i<nTests; ++i) {
        double xs[5] = {0};
        double *cc = tests[i];
        int rc = poly_quintroots(cc[0], cc[1], cc[2], cc[3], cc[4], cc[5], -DBL_MAX, DBL_MAX,
                                 &xs[0], &xs[1], &xs[2], &xs[3], &xs[4]);
        isort(xs, rc);
        if (rc == nRoots[i]) {
            double *rs = roots[i];
            isort(rs, nRoots[i]);
            if (testArrayApprox(rs, xs, rc) == 0) {
                ++nFailed;
                printf("%d) %d / %d: %24.16lf, %24.16lf, %24.16lf, %24.16lf, %24.16lf \n", i, nRoots[i], rc, fabs(xs[0]-roots[i][0]), fabs(xs[1]-roots[i][1]), fabs(xs[2]-roots[i][2]), fabs(xs[3]-roots[i][3]), fabs(xs[4]-roots[i][4]));
            }
        } else {
            ++nMissed;
            ++nFailed;
            printf("%d) %d / %d: %24.16lf, %24.16lf, %24.16lf, %24.16lf, %24.16lf \n", i, nRoots[i], rc, fabs(xs[0]-roots[i][0]), fabs(xs[1]-roots[i][1]), fabs(xs[2]-roots[i][2]), fabs(xs[3]-roots[i][3]), fabs(xs[4]-roots[i][4]));
        }
    }
    if (nMissed==0 && nFailed==0){
        printf("[Quintic polynomial roots]: PASSED.\n");
    } else {
        printf("[Quintic polynomial roots]: Missed roots in %d cases. Failed %d tests in total.\n",
               nMissed, nFailed);
    }
    return nFailed;
}


// ==================================================
// Performance tests
//

long testPerfQuadratic(int N)
{
    int count = 10000;
    double *polys = calloc(count, 3*sizeof(double));
    for(int i=0; i<count; ++i) {
        double a = (double)rand() / (double)RAND_MAX;
        double b = (double)rand() / (double)RAND_MAX;
        double A = (double)rand() / (double)rand();
        double B = (double)rand() / (double)rand();
        double ca = A * B;
        double cb = -(A*b + a*B);
        double cc = a*b;
        polys[i*3] = ca;
        polys[i*3+1] = cb;
        polys[i*3+2] = cc;
    }
    long cr = 0;
    double r1 = 0;
    double r2 = 0, a, b, c;
    for(int i=0; i<N; ++i) {
        a = polys[(i%count)*3];
        b = polys[(i%count)*3+1];
        c = polys[(i%count)*3+2];
        cr += (long)poly_quadroots(a, b, c, 0, 1, &r1, &r2);
    }
    free(polys);
    return cr;
}


long testPerfCubic(int N)
{
    int count = 10000;
    double *polys = calloc(count, 4*sizeof(double));
    for(int i=0; i<count; ++i) {
        double a = (double)rand() / (double)RAND_MAX;
        double b = (double)rand() / (double)RAND_MAX;
        double c = (double)rand() / (double)RAND_MAX;
        double A = (double)rand();
        double cb = A * (-a-b-c);
        double cc = A * (a*b + a*c + b*c);
        double cd = -A * a * b * c;
        polys[i*4] = A;
        polys[i*4+1] = cb;
        polys[i*4+2] = cc;
        polys[i*4+3] = cd;
    }
    long cr = 0;
    double r1 = 0, r2 = 0, r3 = 0, a, b, c, d;
    for(int i=0; i<N; ++i) {
        a = polys[(i%count)*4];
        b = polys[(i%count)*4+1];
        c = polys[(i%count)*4+2];
        d = polys[(i%count)*4+3];
        cr += (long)poly_cubicroots(a, b, c, d, 0, 1, &r1, &r2, &r3);
    }
    free(polys);
    return cr;
}


long testPerfQuartic(int N)
{
    int count = 10000;
    double *polys = calloc(count, 5*sizeof(double));
    for(int i=0; i<count; ++i) {
        double a = (double)rand() / (double)RAND_MAX;
        double b = (double)rand() / (double)RAND_MAX;
        double c = (double)rand() / (double)RAND_MAX;
        double d = (double)rand() / (double)RAND_MAX;
        double A = (double)rand();
        double c1 = -A*(a + b + c + d);
        double c2 = A*(a*b + a*c + b*c + a*d + b*d + c*d);
        double c3 = -A*(a*b*c + a*b*d + a*c*d + b*c*d);
        double c4 = A*a*b*c*d;
        polys[i*5] = A;
        polys[i*5+1] = c1;
        polys[i*5+2] = c2;
        polys[i*5+3] = c3;
        polys[i*5+4] = c4;
    }
    long cr = 0;
    double r1 = 0, r2 = 0, r3 = 0, r4 = 0, a, b, c, d, e;
    for(int i=0; i<N; ++i) {
        a = polys[(i%count)*5];
        b = polys[(i%count)*5+1];
        c = polys[(i%count)*5+2];
        d = polys[(i%count)*5+3];
        e = polys[(i%count)*5+4];
        cr += (long)poly_quartroots(a, b, c, d, e, 0, 1, &r1, &r2, &r3, &r4);
    }
    free(polys);
    return cr;
}


long testPerfQuintic(int N)
{
    int count = 10000;
    double *polys = calloc(count, 6*sizeof(double));
    for(int i=0; i<count; ++i) {
        double a = (double)rand() / (double)RAND_MAX;
        double b = (double)rand() / (double)RAND_MAX;
        double c = (double)rand() / (double)RAND_MAX;
        double d = (double)rand() / (double)RAND_MAX;
        double e = (double)rand() / (double)RAND_MAX;
        double A = (double)rand();
        double ab = a*b;
        double cd = c*d;
        double fc = -A*ab*cd*e;
        double ec = A*(ab*cd + ab*c*e + ab*d*e + a*cd*e + b*cd*e);
        double dc = -A*(ab*c + ab*d + a*cd + b*cd + ab*e + a*c*e + b*c*e + a*d*e + b*d*e + cd*e);
        double cc = A*(ab + a*c + b*c + a*d + b*d + cd + a*e + b*e + c*e + d*e);
        double bc = -A*(a + b + c + d + e);
        polys[i*6] = A;
        polys[i*6+1] = bc;
        polys[i*6+2] = cc;
        polys[i*6+3] = dc;
        polys[i*6+4] = ec;
        polys[i*6+5] = fc;
    }
    long cr = 0;
    double r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0, a, b, c, d, e, f;
    for(int i=0; i<N; ++i) {
        a = polys[(i%count)*6];
        b = polys[(i%count)*6+1];
        c = polys[(i%count)*6+2];
        d = polys[(i%count)*6+3];
        e = polys[(i%count)*6+4];
        f = polys[(i%count)*6+5];
        cr += (long)poly_quintroots(a, b, c, d, e, f, 0, 1, &r1, &r2, &r3, &r4, &r5);
    }
    free(polys);
    return cr;
}




// ==================================================
// Main
//
int main()
{
    // Tests
    printf("\n--------------------------------------------------\n");
    printf("Tests:\n");
    testCubic();
    testQuartic();
    testQuintic();

    // Performance tests
    long ret;
    int N;
    printf("\n--------------------------------------------------\n");
    printf("Performance tests:\n");
    srand(time(NULL));

    // Quadratic
    uint64_t t0 = timestamp();
    N = 1000000;
    ret = testPerfQuadratic(N);
    uint64_t t1 = timestamp();
    double elapsed = get_elapsed_seconds(t0, t1);
    printf("%ld\r", ret);
    print_result(N, "Quadratic", elapsed);

    // Cubic
    t0 = timestamp();
    N = 100000;
    ret = testPerfCubic(N);
    t1 = timestamp();
    elapsed = get_elapsed_seconds(t0, t1);
    printf("%ld\r", ret);
    print_result(N, "Cubic", elapsed);

    // Quartic
    t0 = timestamp();
    N = 100000;
    ret = testPerfQuartic(N);
    t1 = timestamp();
    elapsed = get_elapsed_seconds(t0, t1);
    printf("%ld\r", ret);
    print_result(N, "Quartic", elapsed);

    // Quintic
    t0 = timestamp();
    N = 100000;
    ret = testPerfQuintic(N);
    t1 = timestamp();
    elapsed = get_elapsed_seconds(t0, t1);
    printf("%ld\r", ret);
    print_result(N, "Quintic", elapsed);

    exit(0);
}
