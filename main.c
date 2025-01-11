#include <stdio.h>
#include <stdlib.h>

void print_csr(size_t n, size_t ind, size_t total, int *index, int *col_id, int *values)	{
	printf("index\n");
	for (size_t i = 0; i < ind+1; ++i)	{
		printf("%d ", index[i]);
	}

	printf("\n\n");

	printf("(values col_id)\n");
	for (size_t i = 0; i < total; ++i)	{
		printf("(%d %d) ", values[i], col_id[i]);
	}

	printf("\n\n");
	
	printf("sparse matrix\n");	
	size_t dim = n*n;
	int count = 0;

	for (size_t i = 0; i < ind; ++i)	{
		for (size_t j = 0; j < dim; ++j)	{
			if (col_id[count] == j)	{
				printf("%3d", values[count]);
				count++;
			}
			else
				printf("%3d", 0);

			if ((j+1)%n == 0 && j != 0 && j != dim - 1)
				printf("%3c", '|');
		}
		printf("\n");
		
		if ((i+1)%n == 0 && i != 0 && i != dim - 1)	{
			for (size_t k = 0; k < dim + dim / (n + 1); ++k)	{
				printf("%3c", '-');
			}
			printf("\n");
		}
	}
}

void build_matrix(size_t n, int *index, int *col_id, int *values)	{
	int count = 0;
    index[0] = 0;
    
	for (size_t i = 0; i < n; ++i)	{
        for (size_t j = 0; j < n; ++j)	{
        	if (i > 0) { 
            	values[count] = -1;
				col_id[count] = (i-1) * n + j; 
                count++;
            }

			if (j > 0)	{
				values[count] = -1;
                col_id[count] = i * n + (j-1);
				count++;
			}

			values[count] = 4;
			col_id[count] = i * n + j;
            count++;

			if (j < n-1)	{
				values[count] = -1;
				col_id[count] = i * n + (j+1);
                count++;
			}

			if (i < n-1)	{
				values[count] = -1;
				col_id[count] = (i+1) * n + j;
                count++;
            }

			index[i * n + j + 1] = count;
        }
    }

	index[n*n] = count+1;
}

void jacobi(size_t n, int *index, int *col_id, int *values, int b[n], int x[n], int threshold, size_t max_iter)	{
	//
}

void gs(size_t n, int *index, int *col_id, int *values, int b[n], int x[n], int threshold, size_t max_iter)	{
	//
}

int main(int argc, char *argv[])	{
	if (argc != 2)	{
		perror("expected 1 argument (dim)");
		exit(EXIT_FAILURE);
	}

	size_t n = strtol(argv[1], NULL, 10);

	size_t size = n*n;
	size_t total = 5 * (n*n) - 4 * n;
	
	int *index = aligned_alloc(32, sizeof(int) * size + 1);
	int *col_id = aligned_alloc(32, sizeof(int) * total);
	int *values = aligned_alloc(32, sizeof(int) * total);
	
	if (!index || !col_id || !values)	{
		perror("Error allocation CSR");
		free(index);
		free(col_id);
		free(values);
		exit(EXIT_FAILURE);
	}

	build_matrix(n, index, col_id, values);
	
	print_csr(n, size, total, index, col_id, values);	
	
	free(values);
	free(col_id);
	free(index);

    return 0;
}

