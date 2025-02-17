#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void print_csr(int n, int *index, int *col_id, double *values) {
    printf("index\n");
    int ind = n * n;

    for (int i = 0; i < ind + 1; ++i) {
        printf("%d ", index[i]);
    }

    printf("\n\n");

    printf("(values col_id)\n");
    for (int i = 0; i < index[ind]; ++i) {
        printf("(%.1lf %d) ", values[i], col_id[i]);
    }

    printf("\n\n");

    printf("A =\n");
    int dim = n * n;
    int count = 0;

    for (int i = 0; i < ind; ++i) {
        for (int j = 0; j < dim; ++j) {
            if (col_id[count] == j) {
                printf("%5.1lf", values[count]);
                count++;
            } else
                printf("%5.1lf", 0.0);

            if ((j + 1) % n == 0 && j != 0 && j != dim - 1)
                printf("%c", ' ');
        }
        printf("\n");

        if ((i + 1) % n == 0 && i != 0 && i != dim - 1) {
            printf("\n");
        }
    }
}

void print_vec(int n, double *a) {
    for (int i = 0; i < n; ++i) {
        printf("%lf ", a[i]);
    }
    printf("\n");
}

void init_matrix_r(int n, double a[n], char m) {
    // Random value per entry
    if (m == 'r' || m == 'R') {
        for (int i = 0; i < n; i++) {
            a[i] = (double)RAND_MAX / (double)rand();
        }
    } else {  // Zeroing up the array
        if (m == 'z' || m == 'Z') {
            for (int i = 0; i < n; i++) {
                a[i] = 0.0;
            }
        } else {  // Same value per entry
            if (m == 'c' || m == 'C') {
                double c = (double)RAND_MAX / (double)rand();
                for (int i = 0; i < n; i++) {
                    a[i] = c;
                }
            }
        }
    }
}

void dgemv(int n, int *index, int *col_id, double *values, double *x, double *result) {
    int ind = n * n;

    for (int i = 0; i < ind; ++i) {
        double sum = 0.0;
        for (int j = index[i]; j < index[i + 1]; ++j) {
            sum += values[j] * x[col_id[j]];
        }

        result[i] = sum;
    }
}

double dotprod(int n, double *a, double *b) {
    int ind = n * n;
    double result = 0;
    for (int i = 0; i < ind; ++i) {
        result += a[i] * b[i];
    }

    return result;
}

void axpy(int n, double *a, double x, double y, double *result) {
    for (int i = 0; i < n; ++i) {
        result[i] = a[i] * x + y;
    }
}

void build_matrix(int n, int *index, int *col_id, double *values) {
    int count = 0;
    index[0] = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > 0) {
                values[count] = -1.0;
                col_id[count] = (i - 1) * n + j;
                count++;
            }

            if (j > 0) {
                values[count] = -1.0;
                col_id[count] = i * n + (j - 1);
                count++;
            }

            values[count] = 4.0;
            col_id[count] = i * n + j;
            count++;

            if (j < n - 1) {
                values[count] = -1.0;
                col_id[count] = i * n + (j + 1);
                count++;
            }

            if (i < n - 1) {
                values[count] = -1.0;
                col_id[count] = (i + 1) * n + j;
                count++;
            }

            index[i * n + j + 1] = count;
        }
    }

    index[n * n] = count;
}

double residual_error(int n, int *index, int *col_id, double *values, double b[n], double x[n], double x_kpp[n]) {
    double res_norm = 0.0;
    double b_norm = 0.0;
    double ax = 0.0;
    int ind = n * n;

    for (int i = 0; i < ind; ++i) {
        b_norm += b[i] * b[i];
        ax = 0.0;
        for (int j = index[i]; j < index[i + 1]; ++j) {
            ax += values[j] * x_kpp[col_id[j]];
        }
        res_norm += (b[i] - ax) * (b[i] - ax);
        x[i] = x_kpp[i];
    }

    return sqrt(res_norm) / sqrt(b_norm);
}

void jacobi(int n, int *index, int *col_id, double *values, double b[n], double x[n], double threshold, int max_iter) {
    int ind = n * n;
    double error = threshold + 1.0;
    int iter = 0;

    double *x_kpp = calloc(n * n, sizeof(double));

    while (error > threshold && iter < max_iter) {
        // #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < ind; ++i) {
            x_kpp[i] = b[i];
            int pos = 0;

            for (int j = index[i]; j < index[i + 1]; ++j) {
                int id = col_id[j];
                if (id != i) {
                    x_kpp[i] -= values[j] * x[id];
                } else {
                    pos = j;
                }
            }
            x_kpp[i] /= values[pos];
        }

        error = residual_error(n, index, col_id, values, b, x, x_kpp);
        iter++;
    }

    printf("\nnb iteration: %d/%d\nerror: %lf (%lf)\n", iter, max_iter, error, threshold);

    free(x_kpp);
}

void gs(int n, int *index, int *col_id, double *values, double b[n], double x[n], double threshold, int max_iter) {
    int ind = n * n;
    double error = threshold + 1.0;
    int iter = 0;

    double *x_kpp = calloc(n * n, sizeof(double));

    while (error > threshold && iter < max_iter) {
        for (int i = 0; i < ind; ++i) {
            x_kpp[i] = b[i];
            int pos = 0;

            for (int j = index[i]; j < index[i + 1]; ++j) {
                int id = col_id[j];
                if (id != i) {
                    if (id < i) {
                        x_kpp[i] -= values[j] * x_kpp[id];
                    } else {
                        x_kpp[i] -= values[j] * x[id];
                    }
                } else {
                    pos = j;
                }
            }
            x_kpp[i] /= values[pos];
        }

        error = residual_error(n, index, col_id, values, b, x, x_kpp);
        iter++;
    }

    printf("\nnb iteration: %d/%d\nerror: %lf (%lf)\n", iter, max_iter, error, threshold);
    free(x_kpp);
}

void conj_grad(int n, int *index, int *col_id, double *values, double *b, double *x, double threshold, int max_iter) {
    int ind = n * n;

    double *r = calloc(ind, sizeof(double));
    double *r_kpp = calloc(ind, sizeof(double));
    double *p = calloc(ind, sizeof(double));
    double *ax = calloc(ind, sizeof(double));

    int iter = 0;

    if (!r || !r_kpp || !p || !ax) {
        perror("Error allocation");
        exit(EXIT_FAILURE);
    }

    dgemv(n, index, col_id, values, x, r);

    for (int i = 0; i < ind; ++i) {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    while (iter < max_iter) {
        dgemv(n, index, col_id, values, p, ax);

        double dot_r = dotprod(n, r, r);
        double alpha = dot_r / dotprod(n, ax, p);

        for (int i = 0; i < ind; ++i) {
            x[i] += alpha * p[i];
            r[i] = r[i] - alpha * ax[i];
        }

        double dot_rpp = dotprod(n, r, r);

        if (dot_rpp < threshold * threshold) {
            break;
        }

        double beta = dot_rpp / dot_r;

        for (int i = 0; i < ind; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        iter++;
    }

    double error = residual_error(n, index, col_id, values, b, x, x);

    printf("\nnb iteration: %d/%d\nerror: %lf (%lf)\n", iter, max_iter, error, threshold);

    free(r);
    free(r_kpp);
    free(p);
    free(ax);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        perror("expected 3 argument (dim max_iter threshold)");
        exit(EXIT_FAILURE);
    }

    int n = strtol(argv[1], NULL, 10);
    int max_iter = strtol(argv[2], NULL, 10);
    double threshold = strtod(argv[3], NULL);

    srand(getpid());

    int size = n * n;
    int total = 5 * (n * n) - 4 * n;
    int count = 0;

    int *index = aligned_alloc(32, sizeof(int) * size + 1);
    int *col_id = aligned_alloc(32, sizeof(int) * total);
    double *values = aligned_alloc(32, sizeof(double) * total);

    double *x = aligned_alloc(32, sizeof(double) * size);
    double *b = aligned_alloc(32, sizeof(double) * size);

    double *ax = aligned_alloc(32, sizeof(double) * size);

    if (!index || !col_id || !values || !x || !b || !ax) {
        perror("Error allocation CSR");
        exit(EXIT_FAILURE);
    }

    init_matrix_r(size, b, 'r');

    printf("Building matrix %dx%d ...\n", size, size);
    build_matrix(n, index, col_id, values);

    //------------------------ Jacobi ----------------------------------------

    init_matrix_r(size, x, 'z');

    printf("\nComputing Jacobi ...");
    jacobi(n, index, col_id, values, b, x, threshold, max_iter);

    printf("Checking ... ");
    dgemv(n, index, col_id, values, x, ax);

    for (int i = 0; i < size; ++i) {
        if (fabs(ax[i] - b[i]) < threshold * max_iter) {
            count++;
        }
    }

    if (count == size) {
        printf("passed\n");
    } else {
        printf("failed\n");
        printf("\ndiff:\n");
        int error = 0;
        for (int i = 0; i < size; ++i) {
            if (fabs(ax[i] - b[i]) >= threshold * max_iter) {
                error++;

                if (error == 6) {
                    printf("...\n");
                    break;
                }
                
                printf("%f %f\n", ax[i], b[i]);
            }
        }
    }

    //------------------------ Gauss-Seidel ------------------------------------------

    init_matrix_r(size, x, 'z');

    printf("\nComputing Gauss-Seidel ...");
    gs(n, index, col_id, values, b, x, threshold, max_iter);

    printf("Checking ... ");
    dgemv(n, index, col_id, values, x, ax);

    count = 0;
    for (int i = 0; i < size; ++i) {
        if (fabs(ax[i] - b[i]) < threshold * max_iter) {
            count++;
        }
    }

    if (count == size) {
        printf("passed\n");
    } else {
        printf("failed\n");
        printf("\ndiff:\n");
        int error = 0;
        for (int i = 0; i < size; ++i) {
            if (fabs(ax[i] - b[i]) >= threshold * max_iter) {
                error++;

                if (error == 6) {
                    printf("...\n");
                    break;
                }
                
                printf("%f %f\n", ax[i], b[i]);
            }
        }
    }

    //----------------------- Conjugate Gradient -----------------------------------------

    init_matrix_r(size, x, 'z');

    printf("\nComputing Conjugate Gradient ...");
    conj_grad(n, index, col_id, values, b, x, threshold, max_iter);

    printf("Checking ... ");
    dgemv(n, index, col_id, values, x, ax);

    count = 0;
    for (int i = 0; i < size; ++i) {
        if (fabs(ax[i] - b[i]) < threshold * max_iter) {
            count++;
        }
    }

    if (count == size) {
        printf("passed\n");
    } else {
        printf("failed\n");
        printf("\ndiff:\n");
        int error = 0;
        for (int i = 0; i < size; ++i) {
            if (fabs(ax[i] - b[i]) >= threshold * max_iter) {
                error++;

                if (error == 6) {
                    printf("...\n");
                    break;
                }
                
                printf("%f %f\n", ax[i], b[i]);
            }
        }
    }

    //-------------------------------------------------------------------------------------

    free(values);
    free(col_id);
    free(index);
    free(x);
    free(b);
    free(ax);

    return 0;
}
