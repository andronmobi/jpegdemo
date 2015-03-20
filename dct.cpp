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
#if 0
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
#else
static const float A[] = {
    -122,   49,     66,     41,     41,     43,     40,     38,
    -121,   49,     31,     45,     35,     50,     41,     24,
    -122,   40,     45,     105,    31,     -66,    18,     87,
    -94,    52,     42,     47,     -122,   -122,   8,      51,
    -119,   -23,    53,     51,     45,     70,     61,     42,
    -64,    -122,   -25,    -26,    33,      15,    6,      12,
    -76,    -80,    -64,    -122,   53,     64,     38,     -122,
    -78,    -74,    -84,    -122,   57,     43,     41,     -53,
};
#endif
#endif

static const float quantizCoefs[] = {
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99,
};

static gsl_matrix_float* matrixU;

void print_matrix_float(const gsl_matrix_float* m) {
    for (unsigned int i = 0; i < m->size1; i++) {
        for (unsigned int j = 0; j < m->size2; j++) {
            printf("%.3f\t", gsl_matrix_float_get(m, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix_int(const gsl_matrix_int* m) {
    for (unsigned int i = 0; i < m->size1; i++) {
        for (unsigned int j = 0; j < m->size2; j++) {
            printf("%d\t", gsl_matrix_int_get(m, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void convert_matrix_float_to_int(const gsl_matrix_float* im, gsl_matrix_int* om) {
    for (unsigned int i = 0; i < im->size1; i++) {
        for (unsigned int j = 0; j < im->size2; j++) {
            int elem = lroundf(gsl_matrix_float_get(im, i, j));
            gsl_matrix_int_set(om, i, j, elem);
        }
    }
}

void copy_matrix_int_to_float(const gsl_matrix_int* im, gsl_matrix_float* om) {
    for (unsigned int i = 0; i < im->size1; i++) {
        for (unsigned int j = 0; j < im->size2; j++) {
            float elem = (float) gsl_matrix_int_get(im, i, j);
            gsl_matrix_float_set(om, i, j, elem);
        }
    }
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
    print_matrix_float(matrixC);
    // B = C * U', where B is output matrix
    gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1., matrixC, matrixU, 0., om);
    gsl_matrix_float_free(matrixC);
}

void inv_compute(const gsl_matrix_float* im, gsl_matrix_float* om) {
    // C = U' * B, where B is input matrix
    gsl_matrix_float* matrixC = gsl_matrix_float_alloc(8, 8);
    gsl_blas_sgemm(CblasTrans, CblasNoTrans, 1., matrixU, im, 0., matrixC);
    print_matrix_float(matrixC);
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
    print_matrix_float(matrixA);

    // Create output matrix B
    gsl_matrix_float* matrixB = gsl_matrix_float_alloc(8, 8);

    // B = U * A * U'
    compute(matrixA, matrixB);
    printf("matrix B = U * A * U'\n");
    print_matrix_float(matrixB);

    // Quantization
    _gsl_matrix_float_const_view viewCoefs = gsl_matrix_float_const_view_array(quantizCoefs, 8, 8);
    gsl_matrix_float* matrixCoefs = &viewCoefs.matrix;
    gsl_matrix_float_div_elements(matrixB, matrixCoefs);
    //print_matrix_float(matrixB);
    gsl_matrix_int* matrixQ = gsl_matrix_int_alloc(8, 8);
    convert_matrix_float_to_int(matrixB, matrixQ);
    printf("quantized matrix B\n");
    print_matrix_int(matrixQ);

    // Dequantization
    copy_matrix_int_to_float(matrixQ, matrixB);
    gsl_matrix_float_mul_elements(matrixB, matrixCoefs);
    printf("dequantized matrix B\n");
    print_matrix_float(matrixB);

    // Create output matrix A, restored
    gsl_matrix_float* matrixAR = gsl_matrix_float_alloc(8, 8);
    inv_compute(matrixB, matrixAR);
    printf("matrix A = U' * B * U\n");
    print_matrix_float(matrixAR);


    destroy_dct();

#ifndef CONSTA
    gsl_matrix_float_free(matrixA);
#endif
    gsl_matrix_float_free(matrixAR);
    gsl_matrix_float_free(matrixB);
    gsl_matrix_int_free(matrixQ);

    return 0;
}

