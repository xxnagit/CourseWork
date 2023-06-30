#include<stdio.h>  
#include <stdlib.h>  
#include <time.h>
/********************************************************/  
int main() {  
    static int SIZE = 10;


    int i, j, k,temp;  
    int n = SIZE; 
    int m = n * n;
   
    double a[m];
    double b[m];
    double c[m];

    printf("m = %d\n", m);
  
    for(k=0; k<m; k++) {
        a[k] = 1.23; 
    }

    for(k=0; k<m; k++) {
        b[k] = 4.56; 
    }

     for(k=0; k<m; k++) {
        c[k] = 0; 
    }
    for (k=0; k<m; k++) {
        printf("a[%d]= %f\n", k, a[k]);    
    }

    for (k=0; k<m; k++) {
        printf("b[%d]= %f\n", k, b[k]);   
    }

    for (k=0; k<m; k++) {
    printf("c[%d]= %f\n", k, c[k]);   
    }
    clock_t start, end;
    double cpu_time_used;
    
    start = clock();
    for (k = 0; k < n; k++){ //kij order
        /* For each column j of B */
        for (i = 0; i < n; i++) 
        {
            /* Compute C(i,j) */
            double aik = a[i+k*n];
            
            for(j = 0; j < n; j++)
                c[i+j*n] += aik * b[k+j*n];
        }
    }
    end = clock();

    for (k=0; k<m; k++) {
        printf("c[%d]= %f\n", k, c[k]);
        
    }

   
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("It took this long to calculate matrix mutiplication in a naive way: %f\n",cpu_time_used); 
 
  

    return 0;  
}  
    


    
  
 