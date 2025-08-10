#include <fstream>
#include <ios>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h> //memset
#include <math.h>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "utils.h"

// row-major order (from some slides on MPI) // k is blockX
#define ind(i,j,k) (j)*(k+2)+i



int main(int argv, char** argc)
{
    const int dimensions = 2;
  

    int rank;
    int procs = atoi(argc[3]);
    MPI_Init(&argv, &argc);
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Comm_rank(com, &rank);
    MPI_Comm_size(com, &procs);
    int iterations = atoi(argc[2]);
    int size = atoi(argc[1]);
    if(rank == 0){
        if(argv < 4)
        {
            std::cout <<"Please run as: main.out <size_m*n> <n_iters> <n_procs>\n";
            MPI_Finalize();
            exit(1);
        }
	int arguments[3] = {size, iterations, procs};
	MPI_Bcast(arguments, 3, MPI_INT, 0, com);
	
    }
    //define dimensions 
    int procDims[dimensions];
    memset(procDims, 0, dimensions*sizeof(int));
    MPI_Dims_create(procs, 2, procDims);
    int prox = procDims[0];
    int proy = procDims[1];
    

    //cartesian topology (removes all the boiler plate code and uses default x,y dimensions denoted above, )
    //think of this as thinking of procs/ranks not in a single line but in a x,y pattern. 
    int periods[dimensions];
    memset(periods, 0, dimensions*sizeof(int));
    MPI_Comm cartesionCom;
   //MPI_Cart_create(MPI_Comm old_comm, int ndims, const int *dims, const int *periods, int reorder, MPI_Comm *comm_cart)
    MPI_Cart_create(com, dimensions,  procDims, periods, 0, &cartesionCom);
    //now define actual XY that a proc can use:
    int myCoords[2];
    MPI_Cart_coords(cartesionCom, rank, dimensions, myCoords);
    //my cords in this topology
    int rankX = myCoords[0];
    int rankY = myCoords[1];

    // lets split up problem here:
    int BlockX = size / prox;
    int BlockY = size / proy;
    int offsetX = size%prox;
    int offsetY = size%proy;
    /*if(rank == procs-1)
    {
        BlockX+=offsetX;
        BlockY+=offsetY;
    }*/
    //https://www.mpich.org/static/docs/v3.2/www3/MPI_Cart_shift.html
    int north, south, east, west;
    MPI_Cart_shift(cartesionCom, 0, 1, &west, &east);
    MPI_Cart_shift(cartesionCom, 1, 1, &north, &south);
    /* NS talk to each other and EW talk to each other (lets rank X,Y know who its nearest neigbor is)
    N
    |
    |
    |
    S

    E--------W
    */
  
    int sz = (size+2)*(size+2)*sizeof(double);//(BlockX+2)*(BlockY+2)*sizeof(double);
    int sza = (BlockX+2)*(BlockY+2)*sizeof(double);
    int szb = (BlockX+2)*(BlockY+2);
    double *ArrayOld = (double*)malloc(sz);
     double *Arrayinit = (double*)malloc(sz);
       //define ghosting
    if(rank == 0)
    {
   
    //define edges E/W
    for (int i=0; i<size+2; i++) {
        for (int k =0; k < size+2; k++) {
            if(i == 0 || i == size+1)
                Arrayinit[ind(k, i, size)] = 5.0; //overflow/underflow at 1 for high iterations
            else
                Arrayinit[ind(k, i, size)] = 1.0; //0 results in error (no idea) //possibly cpu or something else.
        }

        
    }

    
    }
     MPI_Scatter(Arrayinit, szb, MPI_DOUBLE, ArrayOld, szb, MPI_DOUBLE, 0, cartesionCom);
  
    double *ArrayNew = (double*)malloc(sz);
    double *extra  = (double*)malloc(sz);
    //finally now the maths part
    double ttime = -MPI_Wtime();
    std::cout << "chunk ranks: {rankX, rankY, GR:global Rank} {" << rankX << "," << rankY << ", GR:" << rank << "} Size: BX : BY " << BlockX << " : " << BlockY << std::endl;
   for(int it = 0; it<iterations; ++it)
   {
       
     // exchange data with neighbors
    double *sbuf =  (double*)malloc(sza); // send buffer (west, east, north, south)
    double *rbuf =  (double*)malloc(sza); // receive buffer (w, e, n, s)
    for(int i=0; i<BlockY; ++i)
        sbuf[i] = ArrayOld[ind(1,i, BlockX)];  //packup West
    if(it ==0)
    /*for(int i=0; i<BlockY; ++i)
        std::cout << sbuf[i] << " ";*/
    for(int i=0; i<BlockY; ++i) 
        sbuf[BlockY+i] = ArrayOld[ind(BlockX,i+1, BlockX)]; // pack east
    for(int i=0; i<BlockX; ++i) 
        sbuf[2*BlockY+i] = ArrayOld[ind(i+1,1, BlockX)]; // pack north
    for(int i=0; i<BlockX; ++i) 
        sbuf[2*BlockY+BlockX+i] = ArrayOld[ind(i+1,BlockY, BlockX)]; // pack south
    int counts[4] = {BlockY, BlockY, BlockX, BlockX};
    int displs[4] = {0, BlockY, 2*BlockY, 2*BlockY+BlockX};
    MPI_Request req;
    MPI_Status status;
   /* if(it ==0)
    std::cout << std::endl;*/
    MPI_Ineighbor_alltoallv(sbuf, counts, displs, MPI_DOUBLE, rbuf, counts, displs, MPI_DOUBLE, cartesionCom, &req);
    MPI_Wait(&req, &status);
    //DEF: MPI_Ineighbor_alltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)

    
    for(int i=0; i<BlockY; ++i) 
        ArrayOld [ind(0,i+1, BlockX)] = rbuf[i];
    for(int i=0; i<BlockY; ++i) 
        ArrayOld [ind(BlockX+1,i+1, BlockX)] = rbuf[BlockY+i]; 
    for(int i=0; i<BlockX; ++i) 
        ArrayOld [ind(i+1,0, BlockX)] = rbuf[2*BlockY+i]; 
    for(int i=0; i<BlockX; ++i) 
        ArrayOld [ind(i+1,BlockY+1, BlockX)] = rbuf[2*BlockY+BlockX+i];

   //update 
    for(int i = 1; i < BlockX+1; i++)
    {
      for(int k = 1; k < BlockY+1; k++)
      {
        //ArrayOld[ind(i,k, BlockX)] + 
        double heat = (ArrayOld[ind(i,k-1, BlockX)] + ArrayOld[ind(i-1,k, BlockX)] + ArrayOld[ind(i+1,k, BlockX)] + ArrayOld[ind(i,k+1, BlockX)]);
        ArrayNew[ind(i,k, BlockX)] =  heat/4.0;
      
      }

  
    }

       
    
    //set arrayNew back to 0 and swap new with old
    extra=ArrayNew;
    ArrayNew = ArrayOld;
    ArrayOld = extra;
  
      
   }

    MPI_Request req2;
    MPI_Status status2;
    int szs = (BlockX+2)*(BlockY+2);
    // std::cout << "Before Gather: " << rankX << "," << rankY << " RR: " << rank << std::endl;
     double *recvArray =  (double*)malloc( sz);
    int err = MPI_Igather(ArrayOld, szb, MPI_DOUBLE, recvArray, szb, MPI_DOUBLE, 0, cartesionCom, &req2); //because catesian
        MPI_Wait(&req2, &status2);
     
//DEF: MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
   std::cout << " P: " << procs << " EE: " << err << std::endl;

        
    
 
   //skip for high iterations (on cluster)&& iterations == 1
   if(rankX ==0 && rankY == 0 && iterations == 1)
   {
     ttime+=MPI_Wtime();
     std::cout << "took: " << ttime<< std::endl;
 

    std::cout << "FILE \n";
    std::ofstream of;
    of.open("output.csv", std::ios::out);
  
    for(int i = 0; i < size; i++)
    {
        for(int j=0; j < size; j++)
        {
 
            of << std::setprecision(9) << std::fixed <<  recvArray[ind(i,j, size)] << ",";

        }
     
        of << "\n"; 
    }
    of.close();
   //more worried about times on cluster and not so much output
   }else if(rankX ==0 && rankY == 0 && iterations > 1) {
      ttime+=MPI_Wtime();
     std::cout << "took: " << ttime<< std::endl;
       std::ofstream of;
        of.open("time.csv", std::ios::out | std::ios::app);
        of << procs << "," << ttime << "\n";
        of.close();
   }

   MPI_Finalize();

return 0;
}

void printArray(double * arr, int blocXize, int blocYize)
{
    for (int i = 0; i < blocXize+2; i++) {
        for (int k = 0; k < blocYize+2; k++) {
            std::cout << arr[(k)*(blocXize+2)+i] << " ";
        }
        std::cout << std::endl;
    
    }
     std::cout << "\n*******************\n" << std::endl;
}
