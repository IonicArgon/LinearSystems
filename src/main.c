#include "functions.h"
#include "sys/time.h"

int main(int argc, char const* argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <matrix-market-filename>\n", argv[0]);
    return 1;
  }

  char method;
  int max_iter;
  bool permute;
  double thresh;
  double omega;

  printf("Select a method:\n");
  printf("1. Jacobi\n");
  printf("2. Successive over-relaxation\n");
  printf("Enter choice: ");
  scanf("%c", &method);

  printf("Enter maximum number of iterations: ");
  scanf("%d", &max_iter);

  printf("Enter threshold: ");
  scanf("%lf", &thresh);

  printf(
      "Do you want to permute the matrix if it is not diagonally dominant? "
      "(1/0): ");
  scanf("%d", &permute);

  if (method == '2') {
    printf("Enter omega: ");
    scanf("%lf", &omega);
  }

  if (method != '1' && method != '2') {
    fprintf(stderr, "Error: Invalid method\n");
    return 1;
  }
  if (max_iter <= 0) {
    fprintf(stderr, "Error: Invalid maximum number of iterations\n");
    return 1;
  }
  if (thresh <= 0.0) {
    fprintf(stderr, "Error: Invalid threshold\n");
    return 1;
  }
  if (method == '2' && (omega <= 0.0 || omega >= 2.0)) {
    fprintf(stderr, "Error: Invalid omega\n");
    return 1;
  }
  if (permute != 0 && permute != 1) {
    fprintf(stderr, "Error: Invalid permute\n");
    return 1;
  }

  struct timeval start, end;
  gettimeofday(&start, NULL);

  const char* filename = argv[1];
  csr_matrix* A = (csr_matrix*)malloc(sizeof(csr_matrix));
  csr_matrix* A_trans = NULL;
  read_csr_matrix(filename, A);
  csr_pretty_print(A);

  char triangularity = csr_triangularity(A);
  if (triangularity == 'N') {
    fprintf(stderr, "Error: Matrix is not triangular\n");
    return 1;
  }

  if (triangularity == 'L') {
    A_trans = (csr_matrix*)malloc(sizeof(csr_matrix));
    A_trans->num_rows = A->num_rows;
    A_trans->num_cols = A->num_cols;
    A_trans->num_non_zeroes = A->num_non_zeroes;
    A_trans->csr_data = (double*)malloc(A->num_non_zeroes * sizeof(double));
    A_trans->col_indices = (int*)malloc(A->num_non_zeroes * sizeof(int));
    A_trans->row_ptrs = (int*)malloc((A->num_rows + 1) * sizeof(int));
    memcpy(A_trans->csr_data, A->csr_data, A->num_non_zeroes * sizeof(double));
    memcpy(A_trans->col_indices, A->col_indices,
           A->num_non_zeroes * sizeof(int));
    memcpy(A_trans->row_ptrs, A->row_ptrs, (A->num_rows + 1) * sizeof(int));
    csr_transpose(A_trans);
  }

  double* b = (double*)malloc(A->num_rows * sizeof(double));
  double* x = (double*)malloc(A->num_rows * sizeof(double));

  for (int i = 0; i < A->num_rows; ++i) {
    b[i] = 1.0;
    x[i] = 0.0;
  }

  if (method == '1')
    solver_iter_jacobi(A, A_trans, b, x, max_iter, thresh, permute);
  else
    solver_iter_SOR(A, A_trans, b, x, max_iter, thresh, omega, permute);

  double residual = comp_residual(A, A_trans, b, x);
  printf("Residual = %lf\n", residual);

  csr_free(A);
  if (triangularity == 'L') {
    csr_free(A_trans);
    free(A_trans);
  }
  free(b);
  free(x);
  free(A);

  gettimeofday(&end, NULL);

  double time_taken =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
  printf("Time taken: %lf seconds\n", time_taken);

  return 0;
}
