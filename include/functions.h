#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "math.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"

typedef struct {
  double* csr_data;
  int* col_indices;
  int* row_ptrs;
  int num_rows;
  int num_cols;
  int num_non_zeroes;
} csr_matrix;

void read_csr_matrix(const char* filename, csr_matrix* A);
void spmv_csr(csr_matrix* A, csr_matrix* A_trans, double* x, double* y);

bool fuzzy_equals(double a, double b, double epsilon);
double comp_norm(double* x, int n);
double comp_residual(csr_matrix* A, csr_matrix* A_trans, double* b, double* x);

void row_permute_preconditioner(csr_matrix* A, double* diagonal);
void solver_iter_jacobi(csr_matrix* A, csr_matrix* A_trans, double* b,
                        double* x, const int max_iter, const double thresh,
                        const bool permute);
void solver_iter_SOR(csr_matrix* A, csr_matrix* A_trans, double* b, double* x,
                     const int max_iter, const double thresh,
                     const double omega, const bool permute);

void csr_pretty_print(const csr_matrix* A);
void csr_free(csr_matrix* A);
char csr_triangularity(const csr_matrix* A);
void csr_transpose(csr_matrix* A);
void csr_row_swap(csr_matrix* A, const int row1, const int row2);
bool csr_is_strictly_diagonally_dominant(const csr_matrix* A,
                                         const csr_matrix* A_trans);

#endif
