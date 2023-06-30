#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <cmath>
#include <signal.h>
#include <unistd.h>
#include <list>
#include "common.h"

using std::vector;
using std::map;
using std::list;

#define CUTOFF 0.01    //Value copied from common.cpp
#define DENSITY 0.0005

double binDim, spaceDim;
int binNo;

inline void binsInit(vector<bin_type>& bins, particle_t* particles, int n) {
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

inline void applyForceToBin(vector<bin_type>& bins, int i, int j, double& dmin, double& davg, int& navg) {
    list<particle_t>::iterator itTarget; 
    list<particle_t>::iterator itNeighbor;

    for (itTarget = bins[i*binNo+j].begin(); itTarget!= bins[i*binNo+j].end(); itTarget++) {
         //for each particle inside the bin, of which calculate the force, initiate its ax=ay=0
        (*itTarget).ax = (*itTarget).ay = 0;
        for (int xNeighbor = -1; xNeighbor <= 1; xNeighbor++) {   //go through 8 adjacent bins and itself: Xneighbor -1~1, yNeighbor -1~1
            for (int yNeighbor = -1; yNeighbor <= 1; yNeighbor++) {
                if (i + xNeighbor >= 0 && i + xNeighbor < binNo && j + yNeighbor >= 0 && j + yNeighbor < binNo) {
                  
                    int indexNeighbor = (i+xNeighbor) * binNo + j + yNeighbor;
                    
                    for (itNeighbor = bins[indexNeighbor].begin(); itNeighbor!= bins[indexNeighbor].end(); itNeighbor++) {//for each particle inside the adjacent bins
                        apply_force(*itTarget, *itNeighbor, &dmin, &davg, &navg);
                    }
                }
            }
        }
    }
}

void binningParticle(particle_t& particle, vector<bin_type>& bins) {//push particles into bins 
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

   

    vector<bin_type> binsWithParcles;
    binsInit(binsWithParcles, particles, n);

    delete[] particles;
    particles = NULL;

    int binsProc = binNo / n_proc; // bins per processor

    // Though each process has all particles, only access particles within private range[ binsStartPr, binsEndPr].
    //Each process owns several rows of the entire space

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

        bin_type lMove; //seperate move depending on the range[binsStartPr, binsEndPr] into localMove and remoteMove
        bin_type rMove;

        for (int i = binsStartPr; i < binsEndPr; ++i) {
            for (int j = 0; j < binNo; ++j) {
                list<particle_t>::iterator itCurr= binsWithParcles[i*binNo + j].begin();
                int vecTail = binsWithParcles[i*binNo + j].size();
                int k = 0;
                while(k < vecTail) { // move particles inside the Current bin
                    move( *itCurr ); 
                    int xIndex = int((*itCurr).x / binDim);  //See if the particle is still inside the previous bin?
                    int yIndex = int((*itCurr).y / binDim);
                    if (binsStartPr <= xIndex && xIndex < binsEndPr) { //check if the move is still within the local process charging area
                        if (xIndex == i && yIndex == j) {  // YES, go on to check the next particle
                            k++;
                            itCurr++;
                        } else{
                            lMove.push_back(*itCurr);
                            vecTail--;
                            itCurr = binsWithParcles[i*binNo + j].erase(itCurr);
                            k++;
                        }
                    } else {
                        rMove.push_back(*itCurr);  // NO, push particle that belongs to a remote vector that would gather later
                        vecTail--;
                        itCurr = binsWithParcles[i*binNo + j].erase(itCurr);//Delete it from the current bin
                        k++;
                    }
                }
            }
        }

        list<particle_t>::iterator itLMove;
        vector<particle_t> rMoveVec, inMoveVec;
        for (itLMove = lMove.begin(); itLMove != lMove.end(); itLMove ++) {
            binningParticle(*itLMove, binsWithParcles);
        }
        if (rank != 0) { //prepare for all remote move of all processes
            for (int i = binsStartPr - 1, j = 0; j < binNo; ++j) {
                bin_type& binTemp = binsWithParcles[i * binNo + j];
                binTemp.clear();
            }
            for (int i = binsStartPr, j = 0; j < binNo; ++j) {
                bin_type& binTemp = binsWithParcles[i * binNo + j];
                list<particle_t>::iterator itBinTemp;
                for (itBinTemp = binsWithParcles[i * binNo + j].begin(); itBinTemp != binsWithParcles[i * binNo + j].end(); itBinTemp++) {
                    rMove.push_back(*itBinTemp);
                    rMoveVec.push_back(*itBinTemp);
                }
                binTemp.clear();
            }
        }

        if (rank != n_proc - 1) {
            for (int i = binsEndPr, j = 0; j < binNo; ++j) {
                bin_type& binTemp = binsWithParcles[i * binNo + j];
                binTemp.clear();
            }
            for (int i = binsEndPr - 1, j = 0; j < binNo; ++j) {
                bin_type& binTemp = binsWithParcles[i * binNo + j];
                list<particle_t>::iterator itBinTemp;
                for (itBinTemp = binsWithParcles[i * binNo + j].begin(); itBinTemp != binsWithParcles[i * binNo + j].end(); itBinTemp++) {
                    rMove.push_back(*itBinTemp);
                    rMoveVec.push_back(*itBinTemp);
                }
                binTemp.clear();
            }
        }

        bin_type inMove; //move coming from other process
        int send_count = rMoveVec.size();
        int recv_counts[n_proc];
        
        MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD); //gather the entire no. of elements of remote move through all processes

        // Got recv_counts 

        int displs[n_proc]; //displs[i] represents the data in the recv comes from i process
        int total = 0;

        if (rank == 0) {
            displs[0] = 0; 
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recv_counts[i-1];
            }
            total = recv_counts[n_proc-1] + displs[n_proc-1];
            
            inMoveVec.resize(total);
        }

        // Got total number

       
        MPI_Gatherv(rMoveVec.data(), send_count, PARTICLE, inMoveVec.data(), recv_counts, displs, PARTICLE, 0, MPI_COMM_WORLD);
        
       
        vector< vector<particle_t> > scParticles; //scatter particles , each process has a vector of scParticles
        scParticles.resize(n_proc);

        if (rank == 0) {
            for (int i = 0; i < inMoveVec.size(); ++i) {
                int x = int(inMoveVec[i].x / binDim);
                int id = min(x / binsProc, n_proc-1);
                scParticles[id].push_back(inMoveVec[i]); //put the incoming move to the corresponding process

                int row = x % binsProc;
                if (row == 0 && id != 0)
                    scParticles[id - 1].push_back(inMoveVec[i]);
                if (row == binsProc-1 && id != n_proc-1)
                    scParticles[id + 1].push_back(inMoveVec[i]);
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
        
        vector<particle_t> outMoveVec;
        outMoveVec.resize(send_count);

        vector<particle_t> scParticlesVec;
        for (int i = 0; i < scParticles.size(); ++i) {
            scParticlesVec.insert(scParticlesVec.end(), scParticles[i].begin(), scParticles[i].end());
        }

        MPI_Scatterv(scParticlesVec.data(), recv_counts, displs, PARTICLE, outMoveVec.data(), send_count, PARTICLE, 0, MPI_COMM_WORLD);

        //scatter outgoing move to each process
    
        for (int i = 0; i < send_count; ++i) { //rebin outgoing move
            particle_t &p = outMoveVec[i];
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
