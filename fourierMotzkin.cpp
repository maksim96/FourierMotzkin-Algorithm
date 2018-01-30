#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <limits>

using std::vector;

void release_memory(int       m,
                    double ** A,
                    double *  b,
                    double *  c)
{
   if (A) {
      int i;
      for (i = 0; i < m; i++) {
         if (A[i]) {
            free(A[i]);
         } else {
            /* If A[i] is NULL, then the next entries won't have been
               allocated and their pointer values will be undefined. */
            break;
         }
      }
      free (A);
   }

   if (b) {
      free(b);
   }
   if (c) {
      free(c);
   }
}

/* The function reads an LP instance from filename. The file
   format is expected to be exactly as in the problem specification.
   On return, *m and *n will be the number of rows and columns,
   respectively; A, b and c will be the the matrix, right-hand side
   vector and objective function vector, respectively. Note that
   the order of the data in the input file is c, then b, then A.
   Memory will be allocated for A, b and c; it is the caller's
   responsibility to release the memory later.
*/
int read_LP(const char * filename,
            int *        m,
            int *        n,
            double ***   A,
            double **    b,
            double **    c)
{
   FILE *fp;

   if (!(fp = fopen(filename, "r")))
   {
      fprintf(stderr, "Could not open \"%s\".\n", filename);
      return EXIT_FAILURE;
   } else {
      int i;
      int j;

      fscanf(fp, "%d %d\n", m, n);

      /* Memory allocation. */
      if ((*A = (double**) malloc(*m * sizeof(double*))) &&
          (*b = (double*)  malloc(*m * sizeof(double)))  &&
          (*c = (double*)  malloc(*n * sizeof(double)))) {
         for (i = 0; i < *m; i++) {
            if (!((*A)[i] = (double*) malloc(*n * sizeof(double)))) {
               fprintf(stderr, "Memory allocation failure.\n");
               release_memory(i - 1, *A, *b, *c);
               fclose(fp);
               return EXIT_FAILURE;
            }
         }
      } else {
         fprintf(stderr, "Memory allocation failure.\n");
         release_memory(0, *A, *b, *c);
         fclose(fp);
         return EXIT_FAILURE;
      }

      /* Copying the values into A, b and c. */
      for (j = 0; j < *n; j++) {
         fscanf(fp, "%lf", *c + j);
      }
      for (i = 0; i < *m; i++) {
         fscanf(fp, "%lf", *b + i);
      }
      for (i = 0; i < *m; i++) {
         for (j = 0; j < *n; j++) {
         /* fscanf(fp, "%lf", (*A)[i]+j); */
            fscanf(fp, "%lf", *(*A+i)+j);
         }
      }

      fclose(fp);
      return EXIT_SUCCESS;
   }
}


void test_it(const char * lp_file)
{
   int       m = 0;
   int       n = 0;
   double ** A;
   double *  b;
   double *  c;

   if (read_LP(lp_file, &m, &n, &A, &b, &c) != EXIT_SUCCESS) {

      fprintf(stderr, "Error parsing \"%s\".\n", lp_file);

   } else {

      int i;
      int j;

      printf("A has %d row%s and %d column%s.\n",
             m, (m == 1 ? "" : "s"),
             n, (n == 1 ? "" : "s"));
      printf("c = [");
      for (j = 0; j < n; j++) {
         printf("%.1lf%s", c[j], (j == n-1 ? "]\n" : ", "));
      }
      printf("transpose(b) = [");
      for (i = 0; i < m; i++) {
         printf("%.1lf%s", b[i], (i == m-1 ? "]\n" : ", "));
      }
      for (i = 0; i < m; i++) {
         printf("A[%d] = [", i);
         for (j = 0; j < n; j++) {
            printf("%5.1lf%s", A[i][j], (j == n-1 ? "]\n" : ", "));
         }
      }

      release_memory(m, A, b, c);
   }
}

void fourierMotzkin(const char * lp_file) {
   int       m = 0;
   int       n = 0;
   double ** A;
   double *  b;
   double *  c;

   if (read_LP(lp_file, &m, &n, &A, &b, &c) != EXIT_SUCCESS) {

      fprintf(stderr, "Error parsing \"%s\".\n", lp_file);

   } else {

      int i;
      int j;

      printf("A has %d row%s and %d column%s.\n",
             m, (m == 1 ? "" : "s"),
             n, (n == 1 ? "" : "s"));
      printf("c = [");
      for (j = 0; j < n; j++) {
         printf("%.1lf%s", c[j], (j == n-1 ? "]\n" : ", "));
      }
      printf("transpose(b) = [");
      for (i = 0; i < m; i++) {
         printf("%.1lf%s", b[i], (i == m-1 ? "]\n" : ", "));
      }
      for (i = 0; i < m; i++) {
         printf("A[%d] = [", i);
         for (j = 0; j < n; j++) {
            printf("%5.1lf%s", A[i][j], (j == n-1 ? "]\n" : ", "));
         }
      }

      //the row count of the input matrix A
      int firstM = m;      

      //the combination of the rows of the last run
      double **yTemp = NULL;

      double *d;

      //the row count of the last (empty) matrix 
      int lastM = 0;


      double ***allD = new double**[n+1];
      double **alld = new double*[n+1];
      double *allm = new double[n+1];

      allD[0] = A;
      alld[0] = b;
      allm[0] = m;

      for (int z = 0; z <= n-1; z++) {

        std::vector<int> N;
        N.reserve(m);
        std::vector<int> Z;
        Z.reserve(m);
        std::vector<int> P;
        P.reserve(m);
        
        for (i = 0; i < m; i++) {
                if (A[i][z]<0) {
                    N.push_back(i);
                } else if (A[i][z] == 0){
                    Z.push_back(i);
                } else {
                    P.push_back(i);
                }
        }

        double** D = new double*[N.size()*P.size()+Z.size()];
        for (i = 0;i<N.size()*P.size()+Z.size();i++) {
            D[i] = new double[n];
        }

        double** y = new double*[N.size()*P.size()+Z.size()];
        for (i = 0;i<N.size()*P.size()+Z.size();i++) {
            y[i] = new double[firstM];
        }

        d = new double[N.size()*P.size()+Z.size()];
        
        if (N.size() > 0 && P.size() > 0) {
            //all pairs of negative and positive rows
            for (i = 0; i < N.size(); i++)  {
                for (j = 0; j < P.size(); j++) {
                        int s = N[i];
                        int t = P[j];
                        for (int k = 0; k < n; k++) {
                            D[i*P.size()+j][k] = A[t][z]*A[s][k] - A[s][z]*A[t][k];
                            if (yTemp == NULL) {
                                y[i*P.size()+j][s] = A[t][z];
                                y[i*P.size()+j][t] = -A[s][z];
                            } else {
                                for (int p = 0; p < firstM; p++) {
                                    y[i*P.size()+j][p] = A[t][z]*yTemp[s][p]-A[s][z]*yTemp[t][p];
                                   // printf("(%d %d): %f = %f*%f - %f*%f\n", s, t, y[i*P.size()+j][k], A[t][z],yTemp[s][k],A[s][z],yTemp[t][k]);
                                }
                                
                            }
                            
                        }
                        d[i*P.size()+j] = A[t][z]*b[s] - A[s][z]*b[t];
                }
            }
        } 

        //just copy the rows with 0 coefficient
        if (Z.size() > 0) {
            int startingPosition = N.size()*P.size();
           // printf("startPos %d\n", startingPosition);
            for (i = 0; i < Z.size(); i++)  {
                
                for (int k = 0; k < n; k++) {
                    D[startingPosition+i][k] = A[Z[i]][k];
                }
                if (yTemp == NULL) {
                    y[startingPosition+i][Z[i]] = 1;   
                } else {
                    for (int k = 0; k < firstM; k++) {
                        y[startingPosition+i][k] = yTemp[Z[i]][k];   
                    }
                }
               
                d[startingPosition+i] = b[Z[i]];
                
            }
        }
        for (i = 0; i < N.size()*P.size()+Z.size(); i++) {
            //check if number is very small. Since it is almost surely an rounding error.
            for (int k = 0; k < n; k++) {
                if (fabs(D[i][k]) < pow(10,-50)) {
                    D[i][k] = 0;
                }
            }
            if (fabs(d[i]) < pow(10,-50)) {
                d[i] = 0;
            }
          /* printf("D[%d] = [", i);
            for (j = 0; j < n; j++) {
                printf("%5.1lf%s", D[i][j], (j == n-1 ? "] " : ", "));
            }
            printf(" | %5.1lf || %d\n",d[i], d[i] < 0);*/
        }
      //  printf("\n");

        A = D;
        b = d;
        yTemp = y;
    /*  printf("================y=================\n");
        for (i = 0; i < N.size()*P.size()+Z.size(); i++) {
            printf("y[%d] = [", i);
            for (j = 0; j < firstM; j++) {
                printf("%5.1lf %s", y[i][j], (j == firstM-1 ? "]\n" : ", "));
            }
        }*/

        allD[z+1] = D;
        alld[z+1] = d;

        m = (int) N.size()*P.size() + Z.size();
      //  printf("m = %d\n", m);
        lastM = m;
        allm[z+1] = m;

        int foundNonZero = 0;


        //"pruning"
        for (i = 0; i < N.size()*P.size()+Z.size(); i++) {
            for (int k = 0; k < n; k++) {
                if (D[i][k] != 0) {
                    foundNonZero = 1;
                    break;
                }
            }   
            if (foundNonZero) {
                break;
            }

        }
        if (!foundNonZero ) {
            break;
        }

    }

   
   

   for (int i = 0; i < lastM; i++) {
       if (d[i] < 0) {

            printf("empty!\n");
            printf("y = (");
                for (int j = 0; j < firstM; j++) {
                    printf("%f %s", yTemp[i][j], (j == firstM-1 ? ")\n" : ", "));
                   
                }
                return;
            }
           
   }
   
   
   double* solutionX = new double[n];

   for (int z = n - 1; z >= 0; z--) {//over all matrices
    double xmax = -std::numeric_limits<double>::infinity();
    double xmin = std::numeric_limits<double>::infinity();
        for (int i = 0; i < allm[z]; i++) { //over all rows
           /* if (z < n-1) {
                printf("%f = %f  ",  alld[z][i],alld[z][i]);
                for (int k = n - 1; k > z; k--) { 
                    printf("-%f * %f ",  allD[z][i][k],solutionX[k]);
                }
                printf("\n");
            }*/
            for (int k = n - 1; k > z; k--) { //substitute and subtract all already calculated x_i from the right hand side
                alld[z][i] -= allD[z][i][k]*solutionX[k];
               // printf("%f = -%f * %f\n",  alld[z][i], allD[z][i][k],solutionX[k]);
            }
          // printf("%f <= %f \n", allD[z][i][z], alld[z][i]);
            if (allD[z][i][z] != 0) {
                    alld[z][i] /=  allD[z][i][z];
            }
            //find feasible interval for x_i
            if (allD[z][i][z] < 0 && alld[z][i] > xmax) {
                xmax = alld[z][i];
            } else if (allD[z][i][z] > 0 && alld[z][i] < xmin) {
                xmin = alld[z][i];
            }

            if  (xmax > -std::numeric_limits<double>::infinity()) {
                solutionX[z] = xmax;
            } else if (xmin < std::numeric_limits<double>::infinity()) {
                    solutionX[z] = xmin;
            }  else {
                solutionX[z] = 0;
            }
            

            
         //  printf("\t%s %s %f \n", (allD[z][i][z] != 0 ? "x" : "0"), ( allD[z][i][z] >= 0 ? "<=" : ">="), alld[z][i]);
        }
        //printf("solution x_n:  %f\n", solutionX[z]);
        
    }
    printf("solvable!\n");
   printf("solution is (");
   for (int z = 0; z < n; z++) {
       printf("%f %s", solutionX[z], (z <= n-2 ? ", ": ")"));
   }
   printf("\n");


/*
   printf("\n");

   for (int i = 0; i < lastM; i++) {
    //  printf("%f <= %f \n", allD[n-2][i][n], alld[n-2][i]);
    for (j = 0; j < n; j++) {
       printf("%f ", allD[1][i][j]);
   }
   printf("\n");
  }*/


   release_memory(m, A, b, c);
  }
}

int main(int argc, const char * argv[])
{
   if (argc < 2) {
      fprintf(stderr, "Usage:  %s  <lp file>\n", argv[0]);
      return EXIT_FAILURE;
   }
   fourierMotzkin(argv[1]);
   return EXIT_SUCCESS;
}