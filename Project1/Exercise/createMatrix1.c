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
    int i, j, k,temp;  
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
    
/*    //(1)二维数组动态内存申请  
    str1 = (double **)malloc(sizeof(double*) * n);//行数  
    for (i=0; i<n; i++){  
        str1[i] = (double *)malloc(sizeof(double) * n);//列数  
    }  
  
    str2 = (double **)malloc(sizeof(double*) * n);//行数  
    for (i=0; i<n; i++){  
        str2[i] = (double *)malloc(sizeof(double) * n);//列数  
    }

    str3 = (double **)malloc(sizeof(double*) * n);//行数  
    for (i=0; i<n; i++){  
        str3[i] = (double *)malloc(sizeof(double) * n);//列数  
    }
      

    //(2)二维数组数据输入  
    for (j=0; j<n; j++){  
        for (i=0; i<n; i++){  
            str1[i][j] = a[i + j*n];  
        }  
    }  
  
    for (j=0; j<n; j++){  
        for (i=0; i<n; i++){  
            str2[i][j] = b[i + j*n];  
        }  
    }    
    
    for (j=0; j<n; j++){  
        for (i=0; i<n; i++){  
            str3[i][j] = 0;  
        }  
    }
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
             printf("str1[%d][%d]= %f\n", i,j, str1[i][j]);
             printf("str2[%d][%d]= %f\n", i,j, str2[i][j]);
             printf("str3[%d][%d]= %f\n", i,j, str3[i][j]);
        }
    }

 */  
    clock_t start, end;
    double cpu_time_used;
    
    start = clock();
    for (i = 0; i < n; i++)
    /* For each column j of B */
    for (j = 0; j < n; j++) 
    {
      /* Compute C(i,j) */
      double cij = c[i+j*n];
      for(k = 0; k < n; k++)
           cij += a[i+k*n] * b[k+j*n];
      c[i+j*n] = cij;
    }

    for (k=0; k<m; k++) {
        printf("c[%d]= %f\n", k, c[k]);
        
    }
    end = clock();
   
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("it took this long to call printf: %f\n",cpu_time_used); 
 
  
    //(3)二维数组内存释放  
 
  
    free(a);
    free(b);
    free(c);
   

/*    for (i=0; i<n; i++){//行数  
        free(str1[i]);//先释放一维指针   
    }  
    free(str1);//最后释放我二维指针  

    for (i=0; i<n; i++){//行数  
        free(str2[i]);//先释放一维指针    
    }  
    free(str2);//最后释放我二维指针 

    for (i=0; i<n; i++){//行数  
        free(str3[i]);//先释放一维指针    
    }  
    free(str3);//最后释放我二维指针  
*/ 
    return 0;  
}  
    
 /*   #include <time.h>
     
     clock_t start, end;
     double cpu_time_used;
     
     start = clock();
     ... // Do the work. 
     end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

*/




    
  
 