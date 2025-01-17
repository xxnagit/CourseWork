#include<stdio.h>  
#include <stdlib.h>  
#include <time.h>
#define min(a, b) (((a) < (b)) ? (a) : (b))
/********************************************************/  
int main() {  
    static int SIZE = 600;
    static int BLOCK_SIZE = 5*5;
    int en = BLOCK_SIZE*(SIZE/BLOCK_SIZE);


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
    /*
    for (k=0; k<m; k++) {
        printf("a[%d]= %f\n", k, a[k]);    
    }

    for (k=0; k<m; k++) {
        printf("b[%d]= %f\n", k, b[k]);   
    }

    for (k=0; k<m; k++) {
    printf("c[%d]= %f\n", k, c[k]);   
    }
    */
    clock_t start, end;
    double cpu_time_used;
    
    start = clock();

    int ii, kk, jj;
    for (i = 0; i<n; i += BLOCK_SIZE){
        for (j = 0; j<n; j += BLOCK_SIZE){
            for (k = 0; k<n; k += BLOCK_SIZE){
                int i_min = min(i+BLOCK_SIZE, n);
                int k_min = min(k+BLOCK_SIZE, n);
                int j_min = min(j+BLOCK_SIZE, n);
                for (ii = i; ii<i_min; ii++){
                    for (jj = j; jj<j_min; jj++){
                        for (kk = k; kk<k_min; kk++){
                            c[ii+jj*n] += a[ii+kk*n] * b[kk+jj*n];
            
                        }    
            
                    }    
        
                }
            }
        }
    }



    end = clock();
/*
    for (k=0; k<m; k++) {
        printf("c[%d]= %f\n", k, c[k]);
        
    }
*/
   
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("It took this long to calculate matrix mutiplication in a Block way: %f\n",cpu_time_used); 
   
  

    return 0;  
}  
    


    
  
 