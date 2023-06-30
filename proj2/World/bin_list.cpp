#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common-v2.h"

using std::vector;

#define CUTOFF 0.01    //Value copied from common.cpp
#define DENSITY 0.0005

double spaceDim; //parameters for bins
int binNo;
double binDim = CUTOFF * 2;   //cutoff*2 works the best comparing with cutoff*1 and cutoff*3

void binsInit(vector<bin_type>& bins, particle_t* particles, int n) {
    spaceDim = sqrt(n * DENSITY);
    binNo = int(ceil(spaceDim / binDim)); 
    int binSpace = binNo * binNo;
    bins.resize(binSpace);

    for (int i = 0; i < n; i++) { //put particles into bins
        int x = int(particles[i].x / binDim);
        int y = int(particles[i].y / binDim);
        bins.[x*binNo + y] = NULL;
    }
    for (int i = 0; i < n; i++) { //put particles into bins
        int x = int(particles[i].x / binDim);
        int y = int(particles[i].y / binDim);
        int binsIndex = x*binNo + y;
        linkedlist_t binsList;
        int count;
        if (bins[binsIndex] == NULL) {
            binsList.head = &particles[i];
            binsList.tail = NULL;
            binsList.nodeNo = 1;
            bins[binsIndex] = binsList.head;
        }
        else {
            node_t* nodeTemp;
            nodeTemp = (node_t*)malloc(sizeof(node_t));
            nodeTemp = &particles[i];
            binsList.tail = nodeTemp;
            binsList.tail->next = NULL;
            binsList.nodeNo++;
        //bins[x*binNo + y].push_back(particles[i]); //push a few particles(2-3) into one bin with index[x*binNo+y]
        }
        count++;
        printf("bins[%d].No = %d\n",binsIndex, binsList.nodeNo);
        printf("Total:%d", count);
    }
}
int main( int argc, char **argv ) 
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    vector<bin_type> binsWithParcles;
    
    bin_type *binTemp;

    set_size( n );
    init_particles( n, particles );
    
    binsInit(binsWithParcles, particles, n);

}