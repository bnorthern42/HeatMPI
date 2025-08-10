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
#define send_tag  1
#define recvTag  2

double **allocate2dArray(int size);
double **initArr(double** arr, int size);

using namespace std;
int main(int argv, char *argc[])
{

  const int dimensions = 2;
  

    int rank;
    int procs = atoi(argc[3]);
    MPI_Init(&argv, &argc);
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    int iterations = atoi(argc[2]);
    int size = atoi(argc[1]);
    int BlockX = size/ procs;
    int offset = size % procs;
    if(rank == procs-1)
    {
        BlockX+=offset;
        
    }
    MPI_Status status;
    if(rank == 0){
        if(argv < 4)
        {
            std::cout <<"Please run as: main.out <size_m*n> <n_iters> <n_procs>\n";
            MPI_Finalize();
            exit(1);
        }
           //  double **arr = allocate2dArray(size);
            // arr = initArr( arr, size);
             
           
           
    }

    int sz = (BlockX+2)*(BlockX+2)*sizeof(double);
    double *ArrayOld = (double*)malloc((BlockX+2)*(BlockX+2)*sizeof(double));
     double *Arrayinit = (double*)malloc((size+2)*(size+2)*sizeof(double));
      if(rank == 0)
    {

   
    //define edges E/W
     for(int i=0; i<size+2; ++i)
            Arrayinit[ind(i,0,size+2)] = 1.0; 
        for(int i=0; i<size+2; ++i) 
            Arrayinit[ind(i,size+1,size+2)] = 1.0;
    /*for (int i=0; i<size+2; i++) {
        for (int k =0; k < size+2; k++) {
    
           
            if(k == 0 || k == size+1)
                Arrayinit[ind(i, k, size)] = 1.0; //overflow/underflow at 1 for high iterations
            else
                Arrayinit[ind(i, k, size)] = 0.0; //0 results in error (no idea) //possibly cpu or something else.
        
        }
    

        
    }*/
   // if(rank ==0)
    //    printArray(Arrayinit, size, size);
    }

      std::cout << "before scatter\n";
      MPI_Scatter(Arrayinit,(BlockX+2)*(BlockX+2), MPI_DOUBLE, ArrayOld, (BlockX+2)*(BlockX+2), MPI_DOUBLE, 0, com);
    
      std::cout << "after scatter\n";
      double *ArrayNew = (double*)malloc((BlockX+2)*(BlockX+2)*sizeof(double));
    double *extra  = (double*)malloc((BlockX+2)*(BlockX+2)*sizeof(double));
    double ttime = -MPI_Wtime();
   for(int itr = 0; itr < iterations; itr++)
    {
        double *sbufnorth = (double*)calloc(size+2,sizeof(double)); // send buffers
        double *sbufsouth = (double*)calloc(size+2,sizeof(double)); 
        double *rbufnorth = (double*)calloc(size+2,sizeof(double));  // send buffers
        double *rbufsouth =(double*)calloc(size+2,sizeof(double)); 
    for(int i=0; i<size+2; ++i) 
        sbufnorth[i] = ArrayOld[(i*size)]; // pack loop - last valid region
    for(int i=0; i<size+2; ++i) 
        sbufsouth[i] = ArrayOld[ind(size+2,i+1,size)]; // pack loop


   if(!rank && !itr)
    {
        cout << "si " << size+1;
        cout << "\nnorth" << std::endl;
        for(int i=0; i<size+2; ++i) 
        std::cout << " " << i*(size+2);
    std::cout << std::endl;
         for(int i=0; i<size+2; ++i) 
            cout << sbufnorth[i] << " "; // pack loop
        cout << "\nsouth" << std::endl;
         for(int i=0; i<size+2; ++i) 
            cout << sbufsouth[i] << " ";; // pack loop
        cout << "\n---------" << std::endl;
    }
    



    //send bottom south
        if(rank!=procs-1){
            MPI_Send(sbufsouth, size+2, MPI_DOUBLE, rank+1, send_tag, com);
            
        }
       
//send top north
        if(rank!=0)
            MPI_Send(sbufnorth, size+2, MPI_DOUBLE, rank-1, send_tag, com);
        
       //recieve both now:
        //send bottom south
        if(rank!=procs-1){
        
            MPI_Recv(sbufnorth, size+2, MPI_DOUBLE, rank+1, send_tag, com, &status);
            for(int i=0; i<size+2; ++i)
                ArrayOld[ind(i,size+1,size)] = rbufnorth[i]; // unpack loop - into ghost cells
        }
       
        //send top north
        if(rank!=0)
        {
            MPI_Recv(sbufsouth, size+2, MPI_DOUBLE, rank-1, send_tag, com, &status); 
            for(int i=0; i<size+2; ++i)
                ArrayOld[ind(i,0,size)] = rbufsouth[i]; // unpack loop - into ghost cells
        }


           // if(rank==0 && itr == 0)
              //  printArray(ArrayOld, size, size);
            //update values
        for(int i = 1; i < BlockX+1; i++)
        {
        for(int k = 1; k < BlockX+1; k++)
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
     int szs = (BlockX+2)*(BlockX+2);
    double *recvArray =  (double*)malloc( (BlockX+2)*(BlockX+2)*sizeof(double));
     int err = MPI_Igather(ArrayOld, szs, MPI_DOUBLE, recvArray, szs, MPI_DOUBLE, 0, com, &req2); //because catesian
        MPI_Wait(&req2, &status);

    if(rank==0 && iterations == 1)
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
   
   }else if(rank==0 && iterations > 1){
    ttime+=MPI_Wtime();
    std::cout << "took" << ttime << std::endl;
   }

   MPI_Finalize();

return 0;


}
// and initilize it for problem
double **allocate2dArray(int size)
{
    double *dat = (double*)malloc(size*size*sizeof(double));
    double **arr = (double**)malloc(size*sizeof(double));
    for(int i = 0; i < size; i++)
        arr[i] = &(dat[size+1]);
   

    return arr;
}
double **initArr(double** arr, int size)
{
     for (int i = 0; i < size; i++) {
        for (int k = 0; k < size; k++) {
            if(i==0 || i == size-1)
                arr[i][k] = 1;
            else
                arr[i][k] = 0;
        }
    
    }
    return arr;
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
