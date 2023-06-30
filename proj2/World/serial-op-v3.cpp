#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <forward_list>
#include <iterator>
#include "common-v3.h"

using std::vector;
using std::forward_list;

#define CUTOFF 0.01    //Value copied from common.cpp
#define DENSITY 0.0005

double spaceDim; //parameters for bins
int binNo;
double binDim = CUTOFF * 2;   //cutoff*2 works the best comparing with cutoff*1 and cutoff*3
//forward_list<particle_t>::iterator it;

void binsInit(vector<bin_type>& bins, particle_t* particles, int n) {
    spaceDim = sqrt(n * DENSITY);
    binNo = int(ceil(spaceDim / binDim)); 
    int binSpace = binNo * binNo;
    bins.resize(binSpace);

    for (int i = 0; i < n; i++) { //put particles into bins
        int x = int(particles[i].x / binDim);
        int y = int(particles[i].y / binDim);
        bins[x*binNo + y].push_front(particles[i]); //push a few particles(2-3) into one bin with index[x*binNo+y]
    }
/* checking for binning initiations with particles
    for (int i = 0; i < binNo; i++) {
        for (int j = 0; j < binNo; j++) {  
            int count = 0;
            for (it = bins[i*binNo + j].begin(); it!= bins[i*binNo + j].end(); it++) {
                count++;
            }
            printf("bins[%d] contains %d particles.\n",i*binNo + j, count);
        }
    } 
*/      
}

//
//  benchmarking program
//
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
    
    bin_type temp;

    set_size( n );
    init_particles( n, particles );
    
    binsInit(binsWithParcles, particles, n);
 // checking for binning initiations with particles
    /*int count = 0;
        for (int i = 0; i < binNo; i++) {
            for (int j = 0; j < binNo; j++) {  
                
                for (it = binsWithParcles[i*binNo + j].begin(); it!= binsWithParcles[i*binNo + j].end(); it++) {
                    count++;
                }
                
            }
        } 
        printf("binsWithParcles initially contains %d particles.\n", count);
    */
    //
    //  simulate a number of time steps
    //
    
    double simulation_time = read_timer( );
    
    for (int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
     
        for (int i = 0; i < binNo; i++) {
            for (int j = 0; j < binNo; j++) {
                 //go through all bins
                
                    forward_list<particle_t>::iterator itTarget; 
                    forward_list<particle_t>::iterator itNeighbor;
             
                    for (itTarget = binsWithParcles[i*binNo+j].begin(); itTarget!= binsWithParcles[i*binNo+j].end(); itTarget++) {
                         //for each particle inside the bin, of which calculate the force, initiate its ax=ay=0
                        (*itTarget).ax = (*itTarget).ay = 0;
                        for (int xNeighbor = -1; xNeighbor <= 1; xNeighbor++) {   //go through 8 adjacent bins and itself: Xneighbor -1~1, yNeighbor -1~1
                            for (int yNeighbor = -1; yNeighbor <= 1; yNeighbor++) {
                                if (i + xNeighbor >= 0 && i + xNeighbor < binNo && j + yNeighbor >= 0 && j + yNeighbor < binNo) {
                                  
                                    int indexNeighbor = (i+xNeighbor) * binNo + j + yNeighbor;
                                    
                                    for (itNeighbor = binsWithParcles[indexNeighbor].begin(); itNeighbor!= binsWithParcles[indexNeighbor].end(); itNeighbor++) {//for each particle inside the adjacent bins
                                        apply_force(*itTarget, *itNeighbor, &dmin, &davg, &navg);
                                    }
                                }
                            }
                        }
                    }

            }
        }
        /*int count_2 = 0;
        for (int i = 0; i < binNo; i++) {
            for (int j = 0; j < binNo; j++) {  
             
                for (it = binsWithParcles[i*binNo + j].begin(); it!= binsWithParcles[i*binNo + j].end(); it++) {
                    count_2++;
                }
                
            }
        } 
        if (count_2 != 1000)
            printf("binsWithParcles count_2 contains %d.\n", count_2);
 */
        for (int i = 0; i < binNo; i++) {
            for(int j = 0; j < binNo; j++) {
                forward_list<particle_t>::iterator itCurr= binsWithParcles[i*binNo + j].begin(); 
                int countCurr = 0;
                //int countPush = 0;
                while ( itCurr!= binsWithParcles[i*binNo + j].end() ) {
                    move(*itCurr);
                    int xIndex = int((*itCurr).x / binDim);  //See if the particle is still inside the previous bin?
                    int yIndex = int((*itCurr).y / binDim);
                    //printf("xIndex: %d, i : %d\n", xIndex,i );
                    if (xIndex == i && yIndex == j) { // YES, go on to check the next particle
                        itCurr++;
                        countCurr++;
                    }
                    else {
                        temp.push_front(*itCurr);
                        //countPush++;
                        itCurr = binsWithParcles[i*binNo + j].erase_after(itCurr--);
                        
                    }
                }
                //binsWithParcles[i*binNo + j].resize(countCurr); //even not resize is working fine
               /*
                if (countPush != 0)
                    printf("There are %d elements pushed in temp\n", countPush );
                */
                /*
                bin_type& vecCurr = binsWithParcles[i * binNo + j];
                int vecTail = vecCurr.size();
                int k = 0;
                
                while(k < vecTail) { // move particles inside the Current bin
                    
                    move( vecCurr[k] ); 
                    int xIndex = int(vecCurr[k].x / binDim);  //See if the particle is still inside the previous bin?
                    int yIndex = int(vecCurr[k].y / binDim);
                    if (xIndex == i && yIndex == j)  // YES, go on to check the next particle
                        k++;
                    else {
                        temp.push_back(vecCurr[k]);  // NO, push particle that belongs to a diff bin into temp vector
                        vecTail--;
                        vecCurr.remove(k); //Delete it from the current bin
                        k++;
                    }
                }
                vecCurr.resize(vecTail); //resize the current bin
                */
            }
        }
/*
        int count_3 = 0;
        for (int i = 0; i < binNo; i++) {
            for (int j = 0; j < binNo; j++) {  
             
                for (it = binsWithParcles[i*binNo + j].begin(); it!= binsWithParcles[i*binNo + j].end(); it++) {
                    count_3++;
                }
                
            }
        } 
        
        printf("binsWithParcles count_3 contains %d.\n", count_3);
 */
        forward_list<particle_t>::iterator itTemp; 
        particle_t pTemp;
        //int count = 0;
        for (itTemp = temp.begin(); itTemp!= temp.end(); itTemp++ ) { // Reallocate the particles into a diff bin
            pTemp = *itTemp;
            int xIndex = int((*itTemp).x / binDim);
            int yIndex = int((*itTemp).y / binDim);
            binsWithParcles[xIndex * binNo + yIndex].push_front(*itTemp);
            //count++;
        }
        temp.clear(); //free the items in temp list

        /*printf("%d elements changed places\n", count);
        ///checking for binning initiations with particles
        int count_over = 0;
        for (int i = 0; i < binNo; i++) {
            for (int j = 0; j < binNo; j++) {  
                
                for (it = binsWithParcles[i*binNo + j].begin(); it!= binsWithParcles[i*binNo + j].end(); it++) {
                    count_over++;
                }
                
            }
        } 
        printf("binsWithParcles afterall %d contains %d.\n", step, count_over);
        printf("-------------------------\n");
        */
        if( find_option( argc, argv, "-no" ) == -1 )
        {
            if (navg) 
            {
                absavg +=  davg/navg;
                nabsavg++;
            }
              
            if (dmin < absmin) 
                absmin = dmin;

            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
        if (nabsavg) absavg /= nabsavg;
        // 
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum ) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    
    if( fsave )
        fclose( fsave );
    
    return 0;
    
}