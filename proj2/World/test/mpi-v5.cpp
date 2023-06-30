#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <cmath>
#include <signal.h>
#include <unistd.h>
#include "common-v5.h"

using std::vector;
using std::map;


#define CUTOFF 0.01    //Value copied from common.cpp
#define DENSITY 0.0005

double binDim, spaceDim;
int binNo;

inline void binsInit(vector<bin_t>& bins, particle_t* particles, int n) {
    spaceDim = sqrt(n * DENSITY);
    binDim = CUTOFF;  
    binNo = int(ceil(spaceDim / binDim)); 
    int binSpace = binNo * binNo;
    bins.resize(binSpace);

    for (int i = 0; i < n; i++) {
        int x = int(particles[i].x / binDim);
        int y = int(particles[i].y / binDim);
        bins[x*binNo + y].push_back(particles[i]);
    }
}

inline void applyForceToBin(vector<bin_t>& bins, int i, int j, double& dmin, double& davg, int& navg) {
    bin_t& vecTarget = bins[i * binNo + j];
    for (int k = 0; k < vecTarget.size(); k++) {
        vecTarget[k].ax = vecTarget[k].ay = 0;
        for (int xNeighbor = -1; xNeighbor <= 1; xNeighbor++) {
            for (int yNeighbor = -1; yNeighbor <= 1; yNeighbor++) {
                if (i + xNeighbor >= 0 && i + xNeighbor < binNo && j + yNeighbor >= 0 && j + yNeighbor < binNo) {
                    bin_t& vecAdjacent = bins[(i+xNeighbor) * binNo + j + yNeighbor];
                    for (int l = 0; l < vecAdjacent.size(); l++) {
                        apply_force( vecTarget[k], vecAdjacent[l], &dmin, &davg, &navg);
                    }    
                }
            }
        }
    }
}

void binningParticle(particle_t& particle, vector<bin_t>& bins) {//push particles into bins 
    int x = particle.x / binDim;
    int y = particle.y / binDim;
    bins[x*binNo + y].push_back(particle);
}


//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //signal(SIGSEGV, sigsegv);

    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //
    //  process command line parameters
    //
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
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


    particle_t *particles = new particle_t[n];
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 ) 
        init_particles( n, particles );
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
    //Each process has all particles

   

    vector<bin_t> binsWithParcles;
    binsInit(binsWithParcles, particles, n);

    delete[] particles;
    particles = NULL;

    int binsProc = binNo / n_proc; // bins per processor

    // Though each process has all particles, only access particles within private range[ binsStartPr, binsEndPr] for each process.

    int binsStartPr = binsProc * rank;
    int binsEndPr = binsProc * (rank + 1);

    if (rank == n_proc - 1)
        binsEndPr = binNo;

    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        // compute local forces
        for (int i = binsStartPr; i < binsEndPr; ++i) {
            for (int j = 0; j < binNo; ++j) {
                applyForceToBin(binsWithParcles, i, j, dmin, davg, navg);
            }
        }

        if (find_option( argc, argv, "-no" ) == -1) {
            
            MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
            
          if (rank == 0){
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) 
                absmin = rdmin;
          }
        }

        bin_t lMove; //seperate move depending on the range[binsStartPr, binsEndPr] into localMove and remoteMove
        bin_t rMove;

        for (int i = binsStartPr; i < binsEndPr; ++i) {
            for (int j = 0; j < binNo; ++j) {
                bin_t& vecCurr = binsWithParcles[i * binNo + j];
                int vecTail = vecCurr.size();
                int k = 0;
                while(k < vecTail) { // move particles inside the Current bin
                    move( vecCurr[k] ); 
                    int xIndex = int(vecCurr[k].x / binDim);  //See if the particle is still inside the previous bin?
                    int yIndex = int(vecCurr[k].y / binDim);
                    if (binsStartPr <= xIndex && xIndex < binsEndPr) { //check if the move is still within the local process charging area
                        if (xIndex == i && yIndex == j)  // YES, go on to check the next particle
                            k++;
                        else{
                            lMove.push_back(vecCurr[k]);
                            vecTail--;
                            vecCurr[k] = vecCurr[vecTail];
                            k++;
                        }
                    } else {
                        rMove.push_back(vecCurr[k]);  // NO, push particle that belongs to a remote vector that would gather later
                        vecTail--;
                        vecCurr[k] = vecCurr[vecTail]; //Delete it from the current bin
                        k++;
                    }
                }
                vecCurr.resize(vecTail); //resize the current bin
            }
        }

        for (int i = 0; i < lMove.size(); ++i) { //rebin the particles inside local_move for all processes
            binningParticle(lMove[i], binsWithParcles);
        }

        if (rank != 0) { //prepare for all remote move of all processes
            for (int i = binsStartPr - 1, j = 0; j < binNo; ++j) {
                bin_t& binTemp = binsWithParcles[i * binNo + j];
                binTemp.clear();
            }
            for (int i = binsStartPr, j = 0; j < binNo; ++j) {
                bin_t& binTemp = binsWithParcles[i * binNo + j];
                rMove.insert(rMove.end(), binTemp.begin(), binTemp.end());
                binTemp.clear();
            }
        }

        if (rank != n_proc - 1) {
            for (int i = binsEndPr, j = 0; j < binNo; ++j) {
                bin_t& binTemp = binsWithParcles[i * binNo + j];
                binTemp.clear();
            }
            for (int i = binsEndPr - 1, j = 0; j < binNo; ++j) {
                bin_t& binTemp = binsWithParcles[i * binNo + j];
                rMove.insert(rMove.end(), binTemp.begin(), binTemp.end());
                binTemp.clear();
            }
        }

        bin_t inMove; //move coming from other process
        int send_count = rMove.size();
        int recv_counts[n_proc];
        
        MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //gather the entire no. of elements of remote move through all processes


        // Got recv_counts 

        int displs[n_proc]; //displs[i] represents the data in the recv comes from i process
        int total = 0;

        if (rank == 0) {
            displs[0] = 0; 
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recv_counts[i-1];
            }
            total = recv_counts[n_proc-1] + displs[n_proc-1];
            
            inMove.resize(total);
        }

        // Got total num 

       
        MPI_Gatherv(rMove.data(), send_count, PARTICLE, inMove.data(), recv_counts, displs, PARTICLE, 0, MPI_COMM_WORLD);
        //the root process got the remote move from other processes
       
        vector<bin_t> scParticles; //scatter particles , each process has a vector of scParticles
        scParticles.resize(n_proc);

        if (rank == 0) {
            for (int i = 0; i < inMove.size(); ++i) {
                int x = int(inMove[i].x / binDim);
                int id = min(x / binsProc, n_proc-1);
                scParticles[id].push_back(inMove[i]); //put the incoming move to the corresponding section of the bin

                int row = x % binsProc;
                if (row == 0 && id != 0)
                    scParticles[id - 1].push_back(inMove[i]);
                if (row == binsProc-1 && id != n_proc-1)
                    scParticles[id + 1].push_back(inMove[i]);
            }
            for (int i = 0; i < n_proc; ++i) {
                recv_counts[i] = scParticles[i].size();
            }
            displs[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recv_counts[i-1];
            }
        }

        send_count = 0;
        
        MPI_Scatter(recv_counts, 1, MPI_INT, &send_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        bin_t outMove;
        outMove.resize(send_count);

        bin_t scParticlesVec;
        for (int i = 0; i < scParticles.size(); ++i) {
            scParticlesVec.insert(scParticlesVec.end(), scParticles[i].begin(), scParticles[i].end());
        }

        MPI_Scatterv(scParticlesVec.data(), recv_counts, displs, PARTICLE, outMove.data(), send_count, PARTICLE, 0, MPI_COMM_WORLD);

        //scatter outgoing move to each process
    
        for (int i = 0; i < send_count; ++i) { //rebin outgoing move
            particle_t &p = outMove[i];
            binningParticle(p, binsWithParcles);
        }

   
    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
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
      if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
