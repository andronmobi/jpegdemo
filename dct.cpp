/*
 * @file dct.cpp
 *
 * @author Andrei Mandychev
 * @date Mars 2015
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define CONSTA

#ifdef CONSTA
static const float A[] = {
    51, 52, 51, 50, 50, 52, 50, 52,
    51, 52, 51, 51, 50, 52, 52, 51,
    50, 50, 51, 52, 52, 51, 51, 51,
    51, 50, 50, 50, 52, 50, 50, 51,
    51, 50, 50, 51, 50, 50, 51, 50,
    50, 51, 52, 52, 51, 50, 50, 50,
    51, 52, 51, 50, 52, 50, 52, 50,
    50, 51, 52, 52, 50, 51, 52, 51,
};
#endif

/*
static const float A[] = {
    -122,   49,     66,     41,     41,     43,     40,     38,
    -121,   49,     31,     45,     35,     50,     41,     24,
    -122,   40,     45,     105,    31,     -66,    18,     87,
    -94,    52,     42,     47,     -122,   -122,   8,      51
    -119,   -23,    53,     51,     45,     70,     61,     42,
    -64,    -122,   -25,    -26,    33,      15,    6,      12,
    -76,    -80,    -64,    -122,   53,     64,     38,     -122,
    -78,    -74,    -84,    -122,   57,     43,     41,     -53,
};*/

static gsl_matrix_float* matrixU;

void print_matrix(const gsl_matrix_float* m) {
    for (unsigned int i = 0; i < m->size1; i++) {
        for (unsigned int j = 0; j < m->size2; j++) {
            printf("%.3f\t", gsl_matrix_float_get(m, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void create_dct() {
    matrixU = gsl_matrix_float_alloc(8, 8);
    // fill elements of the matrix U
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            float elem = (i == 0) ? (sqrt(2) / 2) :
               cosf(M_PI * i * (2 * j + 1 )  / 16);
            gsl_matrix_float_set(matrixU, i, j, elem);
        }
    }
    // 1/2 * U
    gsl_matrix_float_scale(matrixU, 0.5);
}

void destroy_dct() {
    gsl_matrix_float_free(matrixU);
}

void compute(const gsl_matrix_float* im, gsl_matrix_float* om) {
    // C = U * A, where A is input matrix
    gsl_matrix_float* matrixC = gsl_matrix_float_alloc(8, 8);
    gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1., matrixU, im, 0., matrixC);
    print_matrix(matrixC);
    // B = C * U', where B is output matrix
    gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1., matrixC, matrixU, 0., om);
    gsl_matrix_float_free(matrixC);
}

void inv_compute(const gsl_matrix_float* im, gsl_matrix_float* om) {
    // C = U' * B, where B is input matrix
    gsl_matrix_float* matrixC = gsl_matrix_float_alloc(8, 8);
    gsl_blas_sgemm(CblasTrans, CblasNoTrans, 1., matrixU, im, 0., matrixC);
    print_matrix(matrixC);
    // A = C * U, where B is output matrix
    gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1., matrixC, matrixU, 0., om);
    gsl_matrix_float_free(matrixC);
}

int main (int argc, char *argv[]) {

    printf("Discrerte Cosine Transformation\n");

    create_dct();

    // Create input matrix A
#ifdef CONSTA
    _gsl_matrix_float_const_view viewA = gsl_matrix_float_const_view_array(A, 8, 8);
    gsl_matrix_float* matrixA = &viewA.matrix;
#else
    gsl_matrix_float* matrixA = gsl_matrix_float_alloc(8, 8);
    gsl_matrix_float_set_all(matrixA, 100);
#endif
    printf("matrix A\n");
    print_matrix(matrixA);

    // Create output matrix B
    gsl_matrix_float* matrixB = gsl_matrix_float_alloc(8, 8);

    // B = U * A * U'
    compute(matrixA, matrixB);
    printf("matrix B = U * A * U'\n");
    print_matrix(matrixB);

    // Create output matrix A, restored
    gsl_matrix_float* matrixAR = gsl_matrix_float_alloc(8, 8);
    inv_compute(matrixB, matrixAR);
    printf("matrix A = U' * B * U\n");
    print_matrix(matrixAR);

    destroy_dct();

#ifndef CONSTA
    gsl_matrix_float_free(matrixA);
#endif
    gsl_matrix_float_free(matrixAR);
    gsl_matrix_float_free(matrixB);

    return 0;
}

