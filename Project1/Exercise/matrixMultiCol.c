#include<stdio.h>  
#include <stdlib.h>  
#include <time.h>
/********************************************************/  
int main() {  
//  double **str1;//创建二维指针来指向数组  
//  double **str2;  
//  double **str3;
    double *a;
    double *b;
    double *c;
    double temp;
    int i, j, k;  
    int n; 
    int m;
    
    printf("Please input the size n of the n*n matrix, n = ");
    scanf("%d", &n);//str1行数  
    m = n * n;
    printf("m = %d\n", m);
    
    
   

    a = (double *)malloc(sizeof(double) * m);//a[] to store the data of A[][] is column order  
    for(k=0; k<m; k++) {
        a[k] = 1.23; 
    }

    b = (double *)malloc(sizeof(double) * m);//b[] to store the data of B[][] is column order  
    for(k=0; k<m; k++) {
        b[k] = 4.56; 
    }
    c = (double *)malloc(sizeof(double) * m);//b[] to store the data of B[][] is column order  
    for(k=0; k<m; k++) {
        c[k] = 0; 
    }
    for (k=0; k<m; k++) {
        printf("a[%d]= %f\n", k, a[k]);
        
    }
    for (k=0; k<m; k++) {
        printf("b[%d]= %f\n", k, b[k]);
        
    }
    

    clock_t start, end;
    double cpu_time_used;
    
    start = clock();
    for (i = 0; i < n/2; i+= 2){
        /* For each column j of B */
        for (j = 0; j < n/2; j+= 2) 
        {
            /* Compute C(i,j) */
            double cij = c[i+j*n];
            double ciij = c[i+1+j*n];
            double cijj = c[i+(j+1)*n];
            double ciijj = c[i+1+(j+1)*n];
            for(k = 0; k < n/2; k+= 2){
                cij += a[i+k*n] * b[k+j*n];
                ciij += a[i+1+k*n] * b[k+j*n];
                cijj += a[i+k*n] * b[k+(j+1)*n];
                ciijj += a[i+1+(k+1)*n] * b[k+1+(j+1)*n];
            }   
            c[i+j*n] = cij;
            c[i+1+j*n] = ciij;
            c[i+(j+1)*n] = cijj;
            c[i+1+(j+1)*n] = ciijj;
        }
    }

    for (k=0; k<m; k++) {
        printf("c[%d]= %f\n", k, c[k]);
        
    }
    end = clock();
   
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("It took this long to calculate matrix mutiplication in a naive way: %f\n",cpu_time_used); 
 
  
    //(3)二维数组内存释放  
 
  
    free(a);
    free(b);
    free(c);
   


    return 0;  
}  
    
 