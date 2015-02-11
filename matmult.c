//Import Statements
# include <omp.h>
#include<stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

// Function Declarations

void matrix_multiply_row_major(int n1, int n2, int n3, double *a, double *b, double *c);

void matrix_multiply_col_major(int n1, int n2, int n3, double *a, double *b, double *c);

void matrix_multiply_block(int n1, int n2, int n3, double *a, double *b, double *c, int typeBlock);

void matrix_fill_random(int n1, int n2, double *a);

void matrix_print_col_major_toFL(int n1, int n2, double *a);

void matrix_transpose(int n1, int n2, double *a);

void get_walltime(double *wct);


// Initialize variables

FILE *f1; 

FILE *f2;

FILE *f3;

int n1, n2, n3, print_to_file, numEL;

double timed, timed1, timed2;

double *a, *b, *c, *d;

double start, endtime,start1, endtime1, start2, endtime2; 

// Define Constants

# define type 2 // Type of matrix multiplication - 1 = "naive", 2 = Transpose a and multiply, all others = use 2D blocking

# define iterations 100 // Number of iterations to run the multiplication

# define blockSize 4 // Block size
# define numThreads 8 // THe number of threads

/*
 * Main function that takes in either the three matrices dimensions or the 
 * three matrix dimensions, the filenames of a and b and the filename to write c.
 * It then calls the matrix multiplication function based on type and computes c = a*b.
 * It times the required time to run each multiplication and averages it over the number of iterations.
 * n1 - int row dimension of matrix a
 * n2 - int column dimension of matrix a and row dimension of matrix b
 * n3 - int column dimension of matrix b
 * a - an n1 by n2 matrix
 * b - an n2 by n3 matrix
 * c - an n1 by n3  matrix
 * one - 0 if matrix_fill_random correctly populates matrix a, 1 otherwise
 * two - 0 if matrix_fill_random correctly populates matrix b, 1 otherwise
 * three - 0 if matrix_multiply executes correctly, 1 otherwise
 * four - output of matrix_print
 * five - output of matrix_print
 * six - output of matrix_print
 */
 
int main(int argc, const char* argv[]) {
    omp_set_num_threads(numThreads); //Define Number of threads
    // Check input
    if(argc != 4 && argc != 7) {
        printf ("Error - Incorrect Input");
        exit (EXIT_FAILURE);
    }
    
    // Get n1, n2, and n3 from command line
    sscanf (argv[1], "%d", &n1);
    sscanf (argv[2], "%d", &n2);
    sscanf (argv[3], "%d", &n3);
    
    // Create arrays using malloc
    a = malloc(sizeof(double) * (n1 * n2));
    if (a == NULL) {
        printf("malloc error");
        exit(1);
    }
    b = malloc(sizeof(double) * (n2 * n3));
    if (b == NULL) {
        printf("malloc error");
        exit(1);
    }
    c = malloc(sizeof(double) * (n1 * n3));
    if (c == NULL) {
        printf("malloc error");
        exit(1);
    }
    
    // Populate the matrices
    if (argc == 4) {
        matrix_fill_random(n1, n2, a);
        matrix_fill_random(n2, n3, b);
    } else {
        print_to_file = 1;
        
        // Open Files
        f1 = fopen(argv[4], "r");
        f2 = fopen(argv[5], "r");
        f3 = fopen(argv[6], "w");
        
        //Populate the matrix while reading from file
        if (f1 != NULL && f2 != NULL && f3 != NULL) {
            int i = 0;
            double num1, num2, num3;
            while (fscanf(f1, "%lf", &num1) == 1) {
                a[i] = num1;
                i++;
            }
            int j = 0;
            while (fscanf(f2, "%lf", &num2) != EOF) {
                b[j] = num2;
                j++;
            }
        } else {
            printf("Error - failure to open file, %d %d %d\n",f1 != NULL , f2 != NULL , f3 != NULL);
            exit (EXIT_FAILURE);
        }
    }
    
    // Conduct Multiplication 
    
    for(int count = 0; count < iterations; count ++) {
        // Start Timer
         get_walltime(&start);
        
        // Check type of multiplication
        if (type == 1) { 
        
            // Multiply a and b and put it into c using "naive" way.
            matrix_multiply_col_major(n1, n2, n3,a, b, c);
            
            get_walltime(&endtime);
            
            // If output file name given, print to file.
            if ( print_to_file == 1) {
                fprintf(f3,"Matrix c =\n");
                matrix_print_col_major_toFL(n1, n3, c);    
            }
        } else if (type == 2) {
        
            // Multiply a and b and put it into c using transpose a method.
            get_walltime(&start1); // Time transposition
            matrix_transpose(n1, n2, a); // Transpose a
            get_walltime(&endtime1); 
            get_walltime(&start2); // Time multiplication
            matrix_multiply_row_major(n1, n2, n3,a, b, c);
            get_walltime(&endtime2);
            get_walltime(&endtime);
            // If output file name given, print to file.
            if ( print_to_file == 1) {
                fprintf(f3,"Matrix c =\n");
                matrix_print_col_major_toFL(n1, n3, c);    
            }
        } else if (type == 3) {
            
            //// Multiply a and b and put it into c using 2D blocking method.
            matrix_multiply_block(n1, n2, n3, a, b, c, 0);
            get_walltime(&endtime);
            
            // If output file name given, print to file.
            if ( print_to_file == 1) {
                fprintf(f3,"Matrix c =\n");
                matrix_print_col_major_toFL(n1, n3, c);    
            }
        } else {
            matrix_transpose(n1, n2, a); // Transpose a 
            
            //// Multiply a and b and put it into c using 2D blocking method.
            matrix_multiply_block(n1, n2, n3, a, b, c, 1);
            get_walltime(&endtime);
            
            // If output file name given, print to file.
            if ( print_to_file == 1) {
                fprintf(f3,"Matrix c =\n");
                matrix_print_col_major_toFL(n1, n3, c);    
            }
        }
        
        
        // Sum times
        timed += (endtime - start);
        timed1 += (endtime1 - start1);
        timed2 += (endtime2 - start2);
    }
    
    // Free malloc variables
    free(a);
    free(b);
    free(c);
    // Print average times
    printf("********************************************************************************\n");
    printf("To multiply the given two matrices, it takes this computer \n%.12f seconds per matrix multiplication.\n", timed / ( iterations));
    
    // If used transpose print time of transpose and multiplication seperately.
    if (type ==2) {
        printf("Time of transposition: %.12f.\n", (timed1/ ( iterations)));
        printf("To multiply two %d by %d matrices, it takes this computer \n%.12f seconds per matrix multiplication.\n", n1, n3, (timed2/ (iterations)));
    }
}

/*
 * Function that takes in matrix dimensions of two matrices and three matrix pointers.
 * It matrix multiplies the first two matrices and puts the result in the third matrix.
 * It calculates matrix multiplication in a "naive" way and assumes column based storage.
 * n1 - int row dimension of matrix a
 * n2 - int column dimension of matrix a and row dimension of matrix b
 * n3 - int column dimension of matrix b
 * a - double matrix pointer
 * b - double matrix pointer
 * c - double matrix pointer
 * i - row counter of matrix a
 * j - column counter of matrix b
 * k - column counter of a/ row counter of b
 */
 
void matrix_multiply_row_major(int n1, int n2, int n3, double *a, double *b, double *c){
    if(n1 <= 0 || n2 <= 0 || n3 <= 0 || a == NULL || b == NULL || c == NULL) {
                printf ("Error - Can't Multiply Matrices");
                exit (EXIT_FAILURE);
    } else {
        #pragma omp parallel 
        {
        int id,nthrds;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        for(int i = id; i < n1; i+=nthrds) {
            for(int j = 0; j < n3; j++) {
                for(int k = 0; k < n2; k++) {
                    c[i + (j * n1)] += (a[k + (n2 * i)] * b[(j * n2) + k]);
                }
            }
        }
        }
    }
}

/*
 * Function that takes in matrix dimensions of two matrices and three matrix pointers.
 * It matrix multiplies the first two matrices and puts the result in the third matrix.
 * It calculates matrix multiplication by transposing matrix a to achieve better preformance
 * and assumes column based storage and square matrices.
 * n1 - int row dimension of matrix a
 * n2 - int column dimension of matrix a and row dimension of matrix b
 * n3 - int column dimension of matrix b
 * a - double matrix pointer
 * b - double matrix pointer
 * c - double matrix pointer
 * i - row counter of matrix a
 * j - column counter of matrix b
 * k - column counter of a/ row counter of b
 */
 
void matrix_multiply_col_major(int n1, int n2, int n3, double *a, double *b, double *c){
    if(n1 <= 0 || n2 <= 0 || n3 <= 0 || a == NULL || b == NULL || c == NULL) {
        printf ("Error - Can't Multiply Matrices");
                exit (EXIT_FAILURE);
    } else {
    #pragma omp parallel 
    {
            int id,nthrds;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        for(int i = id; i < n1; i+=nthrds) {
            for(int j = 0; j < n3; j++) {
                c[(j * n1) + i] =0;
                for(int k = 0; k < n2; k++) {
                    c[(j * n1) + i] += (a[(k * n1) + i] * b[(j * n2) + k]);
                }

            }
        }}
    }
}

/*
 * Function that takes in matrix dimensions of two matrices and three matrix pointers.
 * It matrix multiplies the first two matrices and puts the result in the third matrix.
 * It calculates matrix multiplication by 2 dimensional blocking and assumes column based storage
 * and square matrices.
 * n1 - int row dimension of matrix a
 * n2 - int column dimension of matrix a and row dimension of matrix b
 * n3 - int column dimension of matrix b
 * a - double matrix pointer
 * b - double matrix pointer
 * c - double matrix pointer
 * i - row counter of matrix a
 * j - column counter of matrix b
 * k - column counter of a/ row counter of b
 * numEL - number of elements in each row of the blockL
 * blockL - number of blocks in a row
 * remainder - number of remaining rows after blocking 
 * index_a - index of first entry in a block of matrix a
 * index_b - index of first entry in a block of matrix b 
 */
 
 void matrix_multiply_block(int n1, int n2, int n3, double *a, double *b, double *c, int typeBlock){
    // Check inputs (make sure correct
    if(n1 <= 0 || n2 <= 0 || n3 <= 0 || a == NULL || b == NULL || c == NULL) {
        printf ("Error - Can't Multiply Matrices");
                exit (EXIT_FAILURE);
    } else {
        int eq = blockSize * ((n2 / blockSize) + ((n2 % blockSize) ? 1 : 0));
        int en = blockSize * ((n3 / blockSize) + ((n3 % blockSize) ? 1 : 0));
#pragma omp parallel 
                    {
        int id,nthrds;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
# pragma omp for        
for (int kk = 0; kk < eq; kk += (blockSize))
        {
            for (int jj = 0; jj < en; jj += blockSize)
            {
                    
                for (int i = 0; i < n1; i++)
                {
                    for (int j = jj; (j < jj + blockSize) && (j < n3); j++)
                    {
                        double sum = c[i + j*n1];
                        for (int k = kk; (k < kk + blockSize) && (k < n2); k++)
                        {
                            if (typeBlock == 1) {
                                sum += a[i*n2 + k]*b[j*n2 + k]; //for matrix A stored in row major order
                            } else {
                                sum += a[i + k*n1]*b[j*n2 + k];// for matrix A stored in column major order
                            }
                        }

                        c[i + j*n1] = sum;
                    }
                }
            }
        }}
    }
}

/*
 * Function that takes in matrix dimensions and a matrix pointer and randomly fills the matrix
 * with numbers between 10 and -10. It returns 0 if the function successfully runs
 * and 1 otherwise.
 * n1 - int row dimanesion
 * n2 - int column dimension
 * a - double matrix pointer
 */
 
void matrix_fill_random(int n1, int n2, double *a) {
    if(n1 <= 0 || n2 <= 0 || a == NULL) {
        printf ("Error - FAilure to Populate");
                exit (EXIT_FAILURE);
    } else {
        for (int i = 0; i < n1 * n2; i++) {
            a[i] = ((rand() / (float)RAND_MAX) * 20) -10;
        }
    }
}

/*
 * Function that takes in matrix dimensions and matrix pointer and prints it
 * in two dimensional form to file. 
 * n1 - int row dimension
 * n2 - int column dimension
 * a - double matrix pointer
 * i - row counter
 * j - column counter
 */
 
void matrix_print_col_major_toFL(int n1, int n2, double *a){
    if(n1 <= 0 || n2 <= 0 || a == NULL) {
        printf ("Error - Can't Print Matrices");
        exit (EXIT_FAILURE);
    } else {
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            fprintf(f3, " %f ", a[i + (j * n1)]);
        }
        fprintf(f3, "\n");
    }
    fprintf(f3,"\n");
    }
}

/*
 * Function that takes in matrix dimensions and matrix pointer and transposes it. 
 * n1 - int row dimension
 * n2 - int column dimension
 * a - double matrix pointer
 * i - row counter
 * j - column counter
 * atrans - temporary storage for transposed values
 */
 
void matrix_transpose(int n1, int n2, double *a){
    if(n1 <= 0 || n2 <= 0 || a == NULL) {
        printf ("Error - Can't Transpose Matrice");
        exit (EXIT_FAILURE);
    } else {
    double* atrans = malloc(sizeof(double) * (n1 * n2));
    if (a == NULL) {
        printf("malloc error");
        exit(1);
    }
    #pragma omp parallel 
    {
    int id,nthrds;
    id = omp_get_thread_num();
    nthrds = omp_get_num_threads();
    for (int i = id; i < n1; i+=nthrds) {
        for (int j = 0; j < n2; j++) {
            atrans[j + (i * n2)] = a[i + (j * n1)];
        }
    }}
    #pragma omp parallel 
    {
    int id,nthrds;
    id = omp_get_thread_num();
    nthrds = omp_get_num_threads();
    for(int i = id; i < n1 * n2; i+=nthrds) {
        a[i] = atrans[i];
    }}
    free(atrans);
	}
}

/*
 * Function that gets walltime. 
 * wct - input storage for time of day

 */
 
 void get_walltime(double *wct) {
  struct timeval tp;
  gettimeofday(&tp,NULL);
  *wct = (double)(tp.tv_sec+tp.tv_usec/1000000.0);
}