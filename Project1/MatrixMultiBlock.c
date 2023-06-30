#include<stdio.h>  
#include <stdlib.h>  
#include <time.h>

#define n 100
#define BlockSize  10
int main()
{
	int i, j, k;

	double a[n][n],b[n][n],c[n][n];
	for(i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			a[i][j] = 1.23;
		}
         
    }

    for(i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			b[i][j] = 4.56;
		}
         
    }

    for(i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			c[i][j] = 0.00;
		}
         
    clock_t start, end;
    double cpu_time_used;
    
    start = clock();
	
	for( i1=0;i1<(n/BlockSize);++i1)
	{
		for(j1=0;j1<(n/BlockSize);++j1)
		{
			for(k1=0;k1<(n/BlockSize);++k1)
			{
				for(i=i1=0;i<min(i1+BlockSize-1);++i)
				{
					for(j=j1=0;j<min(j1+BlockSize-1);++j)
					{
						for(k=k1;k<min(k1+BlockSize-1);++k)
						{               
							c[i][j] = c[i][j] + a[i][k] * b[k][j]
						}
					}
				}
			}
		}
	end = clock();

	for(i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			printf("c[%d][%d]= %f\n", i, j, c[i][j]);
		}
         
   
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("It took this long to calculate matrix mutiplication in a Block way: %f\n",cpu_time_used); 
 
           }
 return 0;
}