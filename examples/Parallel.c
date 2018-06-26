
#include <math.h>
#include <dispatch/dispatch.h>
#include "Polynomial5.h"


typedef struct poly_roots {
    double **polys;
    double **roots;
} poly_roots;


void poly_par_quadroots(double **input, int n_poly, double **output)
{
    size_t n_threads = 8;
    size_t job_size = 1;
    size_t n_poly_z = (size_t)n_poly;
    if ((size_t)n_poly > n_threads) {
        job_size = (size_t)ceil((double)n_poly/(double)n_threads);
    }
    dispatch_queue_t c_queue = dispatch_queue_create("com.polynomial5.cQueue",
                                                     DISPATCH_QUEUE_CONCURRENT);
    dispatch_group_t group = dispatch_group_create();
    poly_roots pr = {input, output};
    dispatch_set_context(c_queue, &pr);
    size_t begin=0;
    size_t end=job_size;
    size_t i;
    for (i = 0; i < n_threads; ++i) {
        dispatch_group_async(group, c_queue, ^() {
                poly_roots *pr = (poly_roots *)dispatch_get_context(c_queue);
                double **polysl = pr->polys;
                double **rootsl = pr->roots;
                double *p;
                double *r;
                size_t j;
                for (j = begin; j < end && j < n_poly_z; ++j) {
                    p = polysl[j];
                    r = rootsl[j];
                    r[0] = (double) poly_quadroots(p[2], p[3], p[4], p[0], p[1], &r[1], &r[2]);
                }
            });
        if (end >= n_poly_z) {
            break;
        }
        begin = end;
        end += job_size;

    }
    dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
    dispatch_release(group);
    dispatch_release (c_queue);
}

/* void poly_par_cubicroots(double **input, int nPoly, double **output); */

/* void poly_par_quartroots(double **input, int nPoly, double **output); */

/* void poly_par_quintroots(double **input, int nPoly, double **output); */
