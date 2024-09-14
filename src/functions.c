#include "functions.h"

typedef struct {
  int row;
  int col;
  double val;
} _Elem;

int _comp(const void* a, const void* b) {
  return ((const _Elem*)a)->row - ((const _Elem*)b)->row;
}

void read_csr_matrix(const char* filename, csr_matrix* A) {
  FILE* file = fopen(filename, "r");

  if (file == NULL) {
    fprintf(stderr, "Error: Unable to open file %s\n", filename);
    return;
  }

  char line[256];
  do {
    char* _ = fgets(line, sizeof(line), file);
  } while (line[0] == '%');

  sscanf(line, "%d %d %d", &A->num_rows, &A->num_cols, &A->num_non_zeroes);

  A->csr_data = (double*)malloc(A->num_non_zeroes * sizeof(double));
  A->col_indices = (int*)malloc(A->num_non_zeroes * sizeof(int));
  A->row_ptrs = (int*)malloc((A->num_rows + 1) * sizeof(int));

  if (A->csr_data == NULL || A->col_indices == NULL || A->row_ptrs == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for matrix\n");
    return;
  }

  _Elem* elems = (_Elem*)malloc(A->num_non_zeroes * sizeof(_Elem));

  if (elems == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for matrix\n");
    return;
  }

  int index = 0;
  while (fgets(line, sizeof(line), file)) {
    int row, col;
    double val;
    sscanf(line, "%d %d %lf", &row, &col, &val);
    elems[index].row = row - 1;
    elems[index].col = col - 1;
    elems[index].val = val;
    index++;
  }

  qsort(elems, A->num_non_zeroes, sizeof(_Elem), _comp);

  for (int i = 0; i < A->num_non_zeroes; ++i) {
    A->csr_data[i] = elems[i].val;
    A->col_indices[i] = elems[i].col;
  }

  int row = 0;
  for (int i = 0; i < A->num_non_zeroes; ++i) {
    if (elems[i].row == row) {
      A->row_ptrs[row + 1]++;
    } else {
      row++;
      A->row_ptrs[row + 1]++;
    }
  }

  for (int i = 1; i <= A->num_rows; ++i) {
    A->row_ptrs[i] += A->row_ptrs[i - 1];
  }

  free(elems);

  fclose(file);
}

void spmv_csr(csr_matrix* A, csr_matrix* A_trans, double* x, double* y) {
  if (A_trans == NULL) {
    for (int i = 0; i < A->num_rows; ++i) {
      for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
        y[i] += A->csr_data[j] * x[A->col_indices[j]];
      }
    }
  } else {
    for (int i = 0; i < A->num_rows; ++i) {
      for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
        if (A->col_indices[j] != i)
          y[i] += A->csr_data[j] * x[A->col_indices[j]];
      }
    }

    for (int i = 0; i < A_trans->num_rows; ++i) {
      for (int j = A_trans->row_ptrs[i]; j < A_trans->row_ptrs[i + 1]; ++j) {
        y[i] += A_trans->csr_data[j] * x[A_trans->col_indices[j]];
      }
    }
  }
}

bool fuzzy_equals(double a, double b, double epsilon) {
  return fabs(a - b) <= epsilon * fmax(fabs(a), fabs(b));
}

double comp_norm(double* x, int n) {
  double norm = 0.0;
  for (int i = 0; i < n; ++i) {
    norm += x[i] * x[i];
  }
  return sqrt(norm);
}

double comp_residual(csr_matrix* A, csr_matrix* A_trans, double* b, double* x) {
  double* A_time_x = (double*)calloc(A->num_rows, sizeof(double));
  if (A_time_x == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for residual\n");
    return -1.0;
  }

  spmv_csr(A, A_trans, x, A_time_x);

  double* residual = (double*)malloc(A->num_rows * sizeof(double));
  if (residual == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for residual\n");
    return -1.0;
  }

  for (int i = 0; i < A->num_rows; ++i) {
    residual[i] = b[i] - A_time_x[i];
  }

  double norm = comp_norm(residual, A->num_rows);

  free(A_time_x);
  free(residual);

  return norm;
}

void row_permute_preconditioner(csr_matrix* A, double* diagonal) {
  int swapped_rows[A->num_rows];
  for (int i = 0; i < A->num_rows; ++i) {
    swapped_rows[i] = 0;
  }

  for (int i = 0; i < A->num_rows; ++i) {
    if (diagonal[i] == 0.0) {
      printf("[COND] Zero diagonal element in row %d\n", i);
      printf("[COND] Looking for row to swap with...\n");
      int row_to_swap = -1;
      for (int j = 0; j < A->num_rows; ++j) {
        for (int k = A->row_ptrs[i]; k < A->row_ptrs[i + 1]; ++k) {
          if (swapped_rows[j]) continue;

          if (A->col_indices[k] == i &&
              fuzzy_equals(A->csr_data[k], 0.0, 1e-6) == false) {
            row_to_swap = j;
            break;
          }
        }
        if (row_to_swap != -1) {
          printf("[COND] Found row %d to swap with\n", row_to_swap);
          break;
        }
      }

      if (row_to_swap != -1 && !swapped_rows[row_to_swap] && !swapped_rows[i]) {
        if (i > row_to_swap) {
          int temp = i;
          i = row_to_swap;
          row_to_swap = temp;
        }
        swapped_rows[i] = 1;

        csr_row_swap(A, i, row_to_swap);

        for (int j = 0; j < A->num_rows; j++) {
          for (int k = A->row_ptrs[j]; k < A->row_ptrs[j + 1]; k++) {
            if (A->col_indices[k] == j) diagonal[j] = A->csr_data[k];
          }
        }
      } else {
        printf("[COND] Unable to find row to swap with\n");
        printf("[COND] The system is singular\n");
        return;
      }
    }
  }
}

void solver_iter_jacobi(csr_matrix* A, csr_matrix* A_trans, double* b,
                        double* x, const int max_iter, const double thresh,
                        const bool permute) {
  double* diagonal = (double*)malloc(A->num_rows * sizeof(double));
  if (diagonal == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for diagonal\n");
    return;
  }

  for (int i = 0; i < A->num_rows; ++i) {
    for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
      if (A->col_indices[j] == i) {
        diagonal[i] = A->csr_data[j];
        break;
      }
    }
  }

  if (permute) {
    printf("[WARN] Checking for zero diagonal elements...\n");
    printf("[WARN] The original matrix might not be preserved...\n");
    row_permute_preconditioner(A, diagonal);
  }

  for (int i = 0; i < A->num_rows; ++i) {
    if (diagonal[i] == 0.0) {
      fprintf(stderr, "Error: Zero diagonal element in row %d\n", i);
      return;
    }
  }

  for (int i = 0; i < A->num_rows; ++i) x[i] = 0.0;

  double* x_prev = (double*)malloc(A->num_rows * sizeof(double));
  if (x_prev == NULL) {
    fprintf(stderr,
            "Error: Unable to allocate memory for previous solutions\n");
    return;
  }

  if (A_trans == NULL) {
    for (int iter = 0; iter < max_iter; ++iter) {
      for (int i = 0; i < A->num_cols; ++i) x_prev[i] = x[i];

      for (int i = 0; i < A->num_rows; ++i) {
        double sum = 0.0;
        for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
          if (A->col_indices[j] != i)
            sum += A->csr_data[j] * x_prev[A->col_indices[j]];
        }
        x[i] = (b[i] - sum) / diagonal[i];
      }

      double residual = comp_residual(A, A_trans, b, x);
      if (residual < thresh) {
        printf("[INFO] Converged in %d iterations\n", iter + 1);
        break;
      }

      printf("\r[INFO] Iteration %d: Residual = %lf\n", iter + 1, residual);
      fflush(stdout);
    }
  } else {
    for (int iter = 0; iter < max_iter; ++iter) {
      for (int i = 0; i < A->num_cols; ++i) x_prev[i] = x[i];

      for (int i = 0; i < A->num_rows; ++i) {
        double sum = 0.0;
        for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
          if (A->col_indices[j] != i)
            sum += A->csr_data[j] * x_prev[A->col_indices[j]];
        }

        for (int j = A_trans->row_ptrs[i]; j < A_trans->row_ptrs[i + 1]; ++j) {
          if (A_trans->col_indices[j] != i)
            sum += A_trans->csr_data[j] * x_prev[A_trans->col_indices[j]];
        }

        x[i] = (b[i] - sum) / diagonal[i];
      }

      double residual = comp_residual(A, A_trans, b, x);
      if (residual < thresh) {
        printf("[INFO] Converged in %d iterations\n", iter + 1);
        break;
      }

      printf("\r[INFO] Iteration %d: Residual = %lf\n", iter + 1, residual);
      fflush(stdout);
    }
  }

  free(x_prev);
  free(diagonal);
}

void solver_iter_SOR(csr_matrix* A, csr_matrix* A_trans, double* b, double* x,
                     const int max_iter, const double thresh,
                     const double omega, const bool permute) {
  double* diagonal = (double*)malloc(A->num_rows * sizeof(double));
  if (diagonal == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for diagonal\n");
    return;
  }

  for (int i = 0; i < A->num_rows; ++i) {
    for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
      if (A->col_indices[j] == i) {
        diagonal[i] = A->csr_data[j];
        break;
      }
    }
  }

  if (permute) {
    printf("[WARN] Checking for zero diagonal elements...\n");
    printf("[WARN] The original matrix might not be preserved...\n");
    row_permute_preconditioner(A, diagonal);
  }

  for (int i = 0; i < A->num_rows; ++i) {
    if (diagonal[i] == 0.0) {
      fprintf(stderr, "Error: Zero diagonal element in row %d\n", i);
      return;
    }
  }

  for (int i = 0; i < A->num_rows; ++i) x[i] = 0.0;

  double* x_prev = (double*)malloc(A->num_rows * sizeof(double));
  if (x_prev == NULL) {
    fprintf(stderr,
            "Error: Unable to allocate memory for previous solutions\n");
    return;
  }

  if (A_trans == NULL) {
    for (int iter = 0; iter < max_iter; ++iter) {
      for (int i = 0; i < A->num_cols; ++i) x_prev[i] = x[i];

      for (int i = 0; i < A->num_rows; ++i) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
          if (A->col_indices[j] < i)
            sum1 += A->csr_data[j] * x[A->col_indices[j]];
          else if (A->col_indices[j] > i)
            sum2 += A->csr_data[j] * x_prev[A->col_indices[j]];
        }
        x[i] = (1.0 - omega) * x_prev[i] +
               (omega / diagonal[i]) * (b[i] - sum1 - sum2);
      }

      double residual = comp_residual(A, A_trans, b, x);
      if (residual < thresh) {
        printf("[INFO] Converged in %d iterations\n", iter + 1);
        break;
      }

      printf("\r[INFO] Iteration %d: Residual = %lf\n", iter + 1, residual);
      fflush(stdout);
    }
  } else {
    for (int iter = 0; iter < max_iter; ++iter) {
      for (int i = 0; i < A->num_cols; ++i) x_prev[i] = x[i];

      for (int i = 0; i < A->num_rows; ++i) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
          if (A->col_indices[j] < i)
            sum1 += A->csr_data[j] * x[A->col_indices[j]];
          else if (A->col_indices[j] > i)
            sum2 += A->csr_data[j] * x_prev[A->col_indices[j]];
        }

        for (int j = A_trans->row_ptrs[i]; j < A_trans->row_ptrs[i + 1]; ++j) {
          if (A_trans->col_indices[j] < i)
            sum1 += A_trans->csr_data[j] * x[A_trans->col_indices[j]];
          else if (A_trans->col_indices[j] > i)
            sum2 += A_trans->csr_data[j] * x_prev[A_trans->col_indices[j]];
        }

        x[i] = (1.0 - omega) * x_prev[i] +
               (omega / diagonal[i]) * (b[i] - sum1 - sum2);
      }

      double residual = comp_residual(A, A_trans, b, x);
      if (residual < thresh) {
        printf("[INFO] Converged in %d iterations\n", iter + 1);
        break;
      }

      printf("\r[INFO] Iteration %d: Residual = %lf\n", iter + 1, residual);
      fflush(stdout);
    }
  }

  free(x_prev);
  free(diagonal);
}

void csr_pretty_print(const csr_matrix* A) {
  printf("- Number of rows: %d\n", A->num_rows);
  printf("- Number of columns: %d\n", A->num_cols);
  printf("- Number of non-zeroes: %d\n", A->num_non_zeroes);

  int max_elem_width = 0;
  for (int i = 0; i < A->num_non_zeroes; ++i) {
    int width = snprintf(NULL, 0, "%lf", A->csr_data[i]);
    if (width > max_elem_width) max_elem_width = width;
  }

  int max_row_width = snprintf(NULL, 0, "%d", A->num_rows);
  int left_padding = max_row_width + 3;

  printf("\n");
  printf("%*s", left_padding, "");
  for (int i = 0; i < A->num_cols; ++i) printf("%*d", max_elem_width + 1, i);
  printf("\n");

  printf("%*s", left_padding, "+-");
  for (int i = 0; i < A->num_cols; ++i)
    for (int j = 0; j < max_elem_width + 1; ++j) printf("-");
  printf("\n");

  int index = 0;
  for (int i = 0; i < A->num_rows; ++i) {
    printf("%*d |", max_row_width, i);
    for (int j = 0; j < A->num_cols; ++j) {
      if (index < A->row_ptrs[i + 1] && A->col_indices[index] == j) {
        printf(" %*lf", max_elem_width, A->csr_data[index]);
        index++;
      } else {
        printf(" %*s", max_elem_width, "");
      }
    }
    printf("\n");
  }
}

void csr_free(csr_matrix* A) {
  free(A->csr_data);
  free(A->col_indices);
  free(A->row_ptrs);
}

char csr_triangularity(const csr_matrix* A) {
  bool upper = false;
  bool lower = false;
  char triangularity = 'N';

  for (int i = 0; i < A->num_rows; ++i) {
    for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
      if (A->col_indices[j] < i) {
        upper = false;
        break;
      } else
        upper = true;
    }
    if (upper == false) break;
  }

  for (int i = 0; i < A->num_rows; ++i) {
    for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
      if (A->col_indices[j] > i) {
        lower = false;
        break;
      } else
        lower = true;
    }
    if (lower == false) break;
  }

  if (upper == true && lower == false)
    triangularity = 'U';
  else if (upper == false && lower == true)
    triangularity = 'L';
  else if (upper == true && lower == true)
    triangularity = 'D';

  return triangularity;
}

void csr_transpose(csr_matrix* A) {
  _Elem* elems = (_Elem*)malloc(A->num_non_zeroes * sizeof(_Elem));
  if (elems == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for transpose\n");
    return;
  }

  int index = 0;
  for (int i = 0; i < A->num_rows; ++i) {
    for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
      elems[index].row = A->col_indices[j];
      elems[index].col = i;
      elems[index].val = A->csr_data[j];
      index++;
    }
  }

  qsort(elems, A->num_non_zeroes, sizeof(_Elem), _comp);

  for (int i = 0; i < A->num_non_zeroes; ++i) {
    A->col_indices[i] = elems[i].col;
    A->csr_data[i] = elems[i].val;
  }

  for (int i = 0; i < A->num_non_zeroes; ++i) {
    A->row_ptrs[elems[i].row + 1]++;
  }

  for (int i = 1; i < A->num_rows; ++i) {
    A->row_ptrs[i + 1] += A->row_ptrs[i];
  }

  free(elems);
}

void csr_row_swap(csr_matrix* A, const int row1, const int row2) {
  double* temp_csr_data = (double*)malloc(A->num_non_zeroes * sizeof(double));
  int* temp_col_indices = (int*)malloc(A->num_non_zeroes * sizeof(int));
  int* temp_row_ptrs = (int*)malloc((A->num_rows + 1) * sizeof(int));

  if (temp_csr_data == NULL || temp_col_indices == NULL ||
      temp_row_ptrs == NULL) {
    fprintf(stderr, "Error: Unable to allocate memory for row swap\n");
    return;
  }

  int temp_csr_index = 0;
  int temp_row_index = 0;
  while (temp_row_index < row1) {
    for (int i = A->row_ptrs[temp_row_index];
         i < A->row_ptrs[temp_row_index + 1]; ++i) {
      temp_csr_data[temp_csr_index] = A->csr_data[i];
      temp_col_indices[temp_csr_index] = A->col_indices[i];
      temp_csr_index++;
    }
    temp_row_ptrs[temp_row_index + 1] = temp_csr_index;
    temp_row_index++;
  }

  for (int i = A->row_ptrs[row2]; i < A->row_ptrs[row2 + 1]; ++i) {
    temp_csr_data[temp_csr_index] = A->csr_data[i];
    temp_col_indices[temp_csr_index] = A->col_indices[i];
    temp_csr_index++;
  }
  temp_row_ptrs[row1 + 1] = temp_csr_index;
  temp_row_index++;

  while (temp_row_index < row2) {
    for (int i = A->row_ptrs[temp_row_index];
         i < A->row_ptrs[temp_row_index + 1]; ++i) {
      temp_csr_data[temp_csr_index] = A->csr_data[i];
      temp_col_indices[temp_csr_index] = A->col_indices[i];
      temp_csr_index++;
    }
    temp_row_ptrs[temp_row_index + 1] = temp_csr_index;
    temp_row_index++;
  }

  for (int i = A->row_ptrs[row1]; i < A->row_ptrs[row1 + 1]; ++i) {
    temp_csr_data[temp_csr_index] = A->csr_data[i];
    temp_col_indices[temp_csr_index] = A->col_indices[i];
    temp_csr_index++;
  }
  temp_row_ptrs[row2 + 1] = temp_csr_index;
  temp_row_index++;

  while (temp_row_index < A->num_rows) {
    for (int i = A->row_ptrs[temp_row_index];
         i < A->row_ptrs[temp_row_index + 1]; ++i) {
      temp_csr_data[temp_csr_index] = A->csr_data[i];
      temp_col_indices[temp_csr_index] = A->col_indices[i];
      temp_csr_index++;
    }
    temp_row_ptrs[temp_row_index + 1] = temp_csr_index;
    temp_row_index++;
  }

  for (int i = 0; i < A->num_non_zeroes; ++i) {
    A->csr_data[i] = temp_csr_data[i];
    A->col_indices[i] = temp_col_indices[i];
  }

  for (int i = 0; i < A->num_rows + 1; ++i) {
    A->row_ptrs[i] = temp_row_ptrs[i];
  }

  free(temp_csr_data);
  free(temp_col_indices);
  free(temp_row_ptrs);
}

bool csr_is_strictly_diagonally_dominant(const csr_matrix* A,
                                         const csr_matrix* A_trans) {
  if (A_trans == NULL) {
    for (int i = 0; i < A->num_rows; ++i) {
      double sum = 0.0;
      for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
        if (A->col_indices[j] != i) sum += fabs(A->csr_data[j]);
      }
      if (fabs(A->csr_data[A->row_ptrs[i]] - sum) <= sum) return false;
    }
  } else {
    double sum = 0.0;
    for (int i = 0; i < A->num_rows; ++i) {
      for (int j = A->row_ptrs[i]; j < A->row_ptrs[i + 1]; ++j) {
        if (A->col_indices[j] != i) sum += fabs(A->csr_data[j]);
      }
      for (int j = A_trans->row_ptrs[i]; j < A_trans->row_ptrs[i + 1]; ++j) {
        if (A_trans->col_indices[j] != i) sum += fabs(A_trans->csr_data[j]);
      }

      if (fabs(A->csr_data[A->row_ptrs[i]] - sum) <= sum) return false;
    }
  }

  return true;
}
