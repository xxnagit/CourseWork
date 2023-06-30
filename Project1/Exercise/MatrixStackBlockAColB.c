#include<stdio.h>  
#include <stdlib.h>  
#include <time.h>
#define min(a, b) (((a) < (b)) ? (a) : (b))
/********************************************************/  
int main() {  
    static int SIZE = 512;
    static int BLOCK_SIZE = 256;


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
    for (k = 0; k < n; k += BLOCK_SIZE){
        for (i = 0; i < n; i += BLOCK_SIZE){
            for (j = 0; j < n; j++){             //jik order
                int i_min = min(i+BLOCK_SIZE, n);
                int k_min = min(k+BLOCK_SIZE, n);
                for (int ii = i; ii < i_min; ii++){ 
                // compute block a[ii][kk] (a[i][k]-a[i_min][k_min]) * col block b[kk][j](b[k][j]-b[k_min][j])
                    for (int kk = k; kk < k_min; kk++) {
                    
                        
                        
                        c[ii+j*n] += a[ii+kk*n] * b[kk+j*n];
                        
                    }//c[ii][j] += a[ii][kk] * b[kk][j]
                }
            }

        }
    }
    end = clock();

    for (k=0; k<m; k++) {
        printf("c[%d]= %f\n", k, c[k]);
        
    }

   
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("It took this long to calculate matrix mutiplication in a BlockA*ColB way: %f\n",cpu_time_used); 
 
  

    return 0;  
}  
    


    
  
 