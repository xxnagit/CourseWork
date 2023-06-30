#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <signal.h>
#include <unistd.h>
#include "common-v2.h"
#include <list>
#include <iterator>

using std::list;
using std::vector;
using std::map;
using std::set;

#define CUTOFF 0.01    //Value copied from common.cpp
#define DENSITY 0.0005

double binDim, spaceDim;
int binNo;

inline void bins_init(vector<bin_type>& bins, particle_t* particles, int n) {
    spaceDim = sqrt(n * DENSITY);
    binDim = CUTOFF * 2;  
    binNo = int(ceil(spaceDim / binDim)); // Should be around sqrt(N/2)
    int binSpace = binNo * binNo;
    bins.resize(binSpace);

    for (int i = 0; i < n; i++) {
        int x = int(particles[i].x / binDim);
        int y = int(particles[i].y / binDim);
        bins[x*binNo + y].push_back(particles[i]);
    }
}

inline void bin_forces(vector<bin_type>& bins, int i, int j, double& dmin, double& davg, int& navg) {
   //go through all bins

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

void bin_particle(particle_t& particle, vector<bin_type>& bins) {
    int x = particle.x / binDim;
    int y = particle.y / binDim;

    bins[x*binNo + y].push_back(particle);
}

/*
inline void get_neighbors(int i, int j, vector<int>& neighbors) {
//push neighbour bins index into a neighbour vector
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0)
                continue;
            if (i + dx >= 0 && i + dx < binNo && j + dy >= 0 && j + dy < binNo) {
                int index = (i + dx) * binNo + j + dy;
                neighbors.push_back(index);
            }
        }
    }
}
*/

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
    
    MPI_Datatype PARTICLE; //Initiate PARTICLE type in MPI
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    //
    //  initialize and distribute the particles into processes
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );

    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

    vector<bin_type> binsWithParcles;
    bins_init(binsWithParcles, particles, n);

    delete[] particles;
    particles = NULL;

    int bins_processor = binNo / n_proc;

    // Though every process has all particles, just access particles within the range [my_bins_start, my_bins_end].

    int my_bins_start = bins_processor * rank;
    int my_bins_end = bins_processor * (rank + 1);

    if (rank == n_proc - 1)
        my_bins_end = binNo;
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        
        // compute local forces
        for (int i = my_bins_start; i < my_bins_end; ++i) {
            for (int j = 0; j < binNo; ++j) {
                bin_forces(binsWithParcles, i, j, dmin, davg, navg);
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

        // move, but not rebin
        bin_type localMove; //move that happens in the local process
        bin_type remoteMove; // moves that happen in the remote processes
/*
        for (int i = my_bins_start; i < my_bins_end; ++i) {
            for (int j = 0; j < binNo; ++j) {
                bin_type& bin = bins[i * binNo + j];
                int tail = bin.size(), k = 0;
                for (; k < tail; ) {
                    move(bin[k]);
                    int x = int(bin[k].x / binDim);
                    int y = int(bin[k].y / binDim);
                    if (my_bins_start <= x && x < my_bins_end) {
                        if (x == i && y == j)
                            ++k;
                        else {
                            local_move.push_back(bin[k]);
                            bin[k] = bin[--tail];
                        }
                    } else {
                        //int who = x / bins_processor;
                        remote_move.push_back(bin[k]);
                        bin[k] = bin[--tail];
                    }
                }
                bin.resize(k);
            }
        }
*/
        for (int i = my_bins_start; i < my_bins_end; ++i) {
            for(int j = 0; j < binNo; j++) {
                list<particle_t>::iterator itCurr= binsWithParcles[i*binNo + j].begin(); 
                int countCurr = 0;
                while ( itCurr!= binsWithParcles[i*binNo + j].end() ) {
                    move(*itCurr);
                    int xIndex = int((*itCurr).x / binDim);  //See if the particle is still inside the previous bin?
                    int yIndex = int((*itCurr).y / binDim);
                    if (my_bins_start <= xIndex && xIndex < my_bins_end) { // check if there is move in local, go on to check the next particle
                        if (xIndex == i && yIndex == j) {
                            itCurr++;
                            countCurr++;
                        }
                        else {
                            localMove.push_back(*itCurr);
                            itCurr = binsWithParcles[i*binNo + j].erase(itCurr);
                        }
                    }
                    else { // not in the range, must be the remote move
                        remoteMove.push_back(*itCurr);
                        itCurr = binsWithParcles[i*binNo + j].erase(itCurr);
                    }
                }
                //binsWithParcles[i*binNo + j].resize(countCurr); //resize the current bin, without resize still works
            }
        }

        list<particle_t>::iterator itLocal;
        for (itLocal = localMove.begin(); itLocal != localMove.end(); itLocal ++) {
            bin_particle(*itLocal, binsWithParcles);
        }

        if (rank != 0) {
            for (int i = my_bins_start - 1, j = 0; j < binNo; ++j) {
                binsWithParcles[i * binNo + j].clear();
               
            }
            for (int i = my_bins_start, j = 0; j < binNo; ++j) {
                //bin_type& bin = bins[i * binNo + j];
                remoteMove.insert(remoteMove.end(), binsWithParcles[i * binNo + j].begin(), binsWithParcles[i * binNo + j].end());
                binsWithParcles[i * binNo + j].clear();
            }
        }

        if (rank != n_proc - 1) {
            for (int i = my_bins_end, j = 0; j < binNo; ++j) {
                binsWithParcles[i * binNo + j].clear();
               
            }
            for (int i = my_bins_end - 1, j = 0; j < binNo; ++j) {
                //bin_type& bin = bins[i * binNo + j];
                remoteMove.insert(remoteMove.end(), binsWithParcles[i * binNo + j].begin(), binsWithParcles[i * binNo + j].end());
                binsWithParcles[i * binNo + j].clear();
            }
        }

        // int len_ = remote_move.size();
        // int total_ = 0;
        // MPI_Reduce(&len_, &total_, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        vector<particle_t> incomingMove_vector;
        
        int send_count = remoteMove.size();
        int recv_counts[n_proc];

        // printf("worker: %d. MPI_Gather.\n", rank);
        MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // now root gets recv_counts

        int disc[n_proc];
        int total_num = 0;

        if (rank == 0) {
            disc[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                disc[i] = disc[i-1] + recv_counts[i-1];
            }
            total_num = recv_counts[n_proc-1] + disc[n_proc-1];
            // printf("worker: %d, 1. %d / %d.\n", rank, total_, total_num);
            // assert(total_ == total_num);
            incomingMove_vector.resize(total_num);
        }

        // now root knows total_num.

        //printf("worker: %d. MPI_Gatherv.\n", rank);
        vector<particle_t> remoteMove_vector;
        remoteMove_vector.resize(remoteMove.size());
        list<particle_t>::iterator itRemoteMove = remoteMove.begin();;
        for (int i = 0; i < remoteMove.size(); i++) {
            remoteMove_vector[i] = (*itRemoteMove);
            if (itRemoteMove != remoteMove.end())
                itRemoteMove++;
        }

        MPI_Gatherv(remoteMove_vector.data(), send_count, PARTICLE, 
            incomingMove_vector.data(), recv_counts, disc, PARTICLE, 
            0, MPI_COMM_WORLD);

        //printf("worker: %d. Classify.\n", rank);

        vector<bin_type> scatter_particles; //every process has a list of scatter_particles
        list<particle_t>::iterator itScatterP;
        scatter_particles.resize(n_proc);
        bin_type incomingMove;
        list<particle_t>::iterator itInMove;
        for (int i = 0; i < incomingMove_vector.size(); i++) { //put into another container
           incomingMove.push_back(incomingMove_vector[i]);
        }


        if (rank == 0) {
            for (itInMove = incomingMove.begin(); itInMove!= incomingMove.end(); itInMove++) {
            //for (int i = 0; i < incoming_move.size(); ++i) {
                int x = int((*itInMove).x / binDim);

                assert((*itInMove).x >= 0 && (*itInMove).y >= 0 && (*itInMove).x <= spaceDim && (*itInMove).y <= spaceDim);

                int which = min(x /bins_processor, n_proc-1);
                scatter_particles[which].push_back((*itInMove));

                int row = x % bins_processor;
                if (row == 0 && which != 0)
                    scatter_particles[which - 1].push_back((*itInMove));
                if (row == bins_processor-1 && which != n_proc-1)
                    scatter_particles[which + 1].push_back((*itInMove));
            }
            for (int i = 0; i < n_proc; ++i) {
                recv_counts[i] = scatter_particles[i].size();
            }
            disc[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                disc[i] = disc[i-1] + recv_counts[i-1];
            }
            // printf("worker: %d, 2. %d / %d.\n", rank, total_, displs[n_proc-1] + recv_counts[n_proc-1]);
            // assert(total_ == displs[n_proc-1] + recv_counts[n_proc-1]);
        }

        // printf("worker: %d. MPI_Scatter.\n", rank);
        send_count = 0;
        MPI_Scatter(recv_counts, 1, MPI_INT, &send_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        bin_type outgoingMove;
        list<particle_t>::iterator itOutMove;
        outgoingMove.resize(send_count);
        vector<particle_t> outgoingMove_vector;
        outgoingMove_vector.resize(send_count);

        bin_type scatter_particles_list;
        list<particle_t>::iterator itScatterPL = scatter_particles_list.begin();
        for (int i = 0; i < scatter_particles.size(); i++) {
            for (itScatterP = scatter_particles[i].begin(); itScatterP!= scatter_particles[i].end(); itScatterP++) {
                scatter_particles_list.push_back(*itScatterP);
            }
        }

        // printf("worker: %d. MPI_Scatterv.\n", rank);
        vector<particle_t> scatter_particles_vector;
        scatter_particles_vector.resize(scatter_particles_list.size());

        for (int i = 0; i < scatter_particles_list.size(); i++) { //put the list of scatter_particle into a vector container to get ready for MPI_Scatterv
            scatter_particles_vector[i] = (*itScatterPL);
            if (itScatterPL != scatter_particles_list.end())
                itScatterPL++; 
        }

        MPI_Scatterv(scatter_particles_vector.data(), recv_counts, disc, PARTICLE, 
            outgoingMove_vector.data(), send_count, PARTICLE, 0, MPI_COMM_WORLD);
        for (int i = 0; i < outgoingMove_vector.size(); i++) {
            outgoingMove.push_back(outgoingMove_vector[i]);
        }
        
        for (itOutMove = outgoingMove.begin(); itOutMove != outgoingMove.end(); itOutMove++) {
            particle_t p = (*itOutMove);
            assert(p.x >= 0 && p.y >= 0 && p.x <= spaceDim && p.y <= spaceDim);
            bin_particle(p, binsWithParcles);
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
