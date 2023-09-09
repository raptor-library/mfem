// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Rosenbrock benchmark example

/**
 * xval: design variables
 * fval: objective value
 * dfdx: gradient of objective
 * gx: constraint values
 * dgdx: gradient of constraints
 * xmin: lower bound on design variables
 * xmax: upper bound on design variables
 * norm2: norm of KKT residual
 * normInf: infinity norm of KKT residual
 * xo1: old design variables
 * xo2: old design variables
 * U: upper asymptote
 * L: lower asymptote
 * init: initial asymptote
 * decrease: decrease factor
 * increase: increase factor
 * n: number of design variables
 * m: number of constraints
 * a: Rosenbrock parameter
 * b: Rosenbrock parameter
*/

#include "mtop_MMA.hpp"
#include <fstream>
//#include <mpi.h>

//using namespace mfem;
using namespace mma;

void Rosenbrock(double* xval, double a, double b, double* const fval,
                double* const dfdx, double* const gx, double* const dgdx);
void MultiVar_Rosenbrock(double* xval, double a, double b, int nVar, int nCon, double* const fval,
                double* const dfdx, double* const gx, double* const dgdx);

int main(int argc, char *argv[])
{
   int nVar = 20;
   int nCon = 10;

   // mfem::MPI_Session mpi(argc, argv); // Initialize MPI (MPI Init)
   // int rank, size;
   // MPI_Comm communicator = MPI_COMM_WORLD;
   // MPI_Comm_size(communicator, &size);
   // MPI_Comm_rank(communicator, &rank);

   int nVarGlobal = nVar;
   // nVar = nVar / size;
   // if (nVar % size != 0)
   // {
   //    if (rank == 0)
   //    {
   //       std::cout << "Number of variables must be divisible by number of processors" << std::endl;
   //    }
   //    MPI_Finalize();
   //    exit(1);
   // }
   
   
   //Rosenbrock parameters
   double aR = 1.0;
   double bR = 100.0;
   // Objective parameters
   double epsimin = 0.0000001;
   int sx = 1;

   double* xval = new double[nVar];

   double* const fval = new double[1];
   double* const dfdx = new double[nVarGlobal];
   double* const gx = new double[nCon];
   double* const dgdx = new double[nVarGlobal * nCon];
   double* xmin = new double[nVarGlobal];
   double* xmax = new double[nVarGlobal];
   double* upper = new double[nVarGlobal];
   double* lower = new double[nVarGlobal];
   // Simulation parameters
   int maxiter = 2;
   int restart = maxiter + 1;
   double norm2 = 0.0;
   double normInf = 0.0;
   double kkttol = 0.0;
   double tol = 1.0;
   // -- RESTART --
   // if a restart from a previous run is desired, set iter to the iteration
   // number of the previous run
   int iter = 0;
   std::string name = "Restart.dat";

   // ----------- SET UP PROBLEM ------------
   for (int i = 0; i < nVarGlobal; i++)
   {
      xval[i] = 0.0;
      xmin[i] = -2.0;
      xmax[i] = 2.0;
   }

   // ---------- CALL CONSTRUCTOR -----------
   MMA MMAmain(nVarGlobal, nCon, xval, xmin, xmax, "classic" );
   std::ofstream mma;
   mma.open("mma.dat");
   //if (rank == 0)
   //{
   mma << xval[0] << "\n" << xval[1] << "\n" << xmax[0] << "\n" << xmax[1] << "\n"
      << xmin[0] << "\n" << xmin[1] << std::endl;
   //}
   
   
   //Rosenbrock(xval, aR, bR, fval, dfdx, gx, dgdx);
   MultiVar_Rosenbrock(xval, aR, bR, nVarGlobal, nCon, fval, dfdx, gx, dgdx);


   while (iter < maxiter)
   {

      // ----------------- RUN MMA --------------------
      MMAmain.Update(iter, dfdx, gx, dgdx, xval);
      // ---------- UPDATE DESIGN & CONSTRAINTS -------
      //Rosenbrock(xval, aR, bR, fval, dfdx, gx, dgdx);
      MultiVar_Rosenbrock(xval, aR, bR, nVarGlobal, nCon, fval, dfdx, gx, dgdx);
      // ---------- PRINT FOR VISUALIZATION------------
      lower = MMAmain.getLow();
      upper = MMAmain.getUpp();
      // if (rank == 0)
      // {
         mma << xval[0] << "\n" << xval[1] << "\n" << upper[0] << "\n" << upper[1] <<
         "\n" << lower[0] << "\n" << lower[1] << std::endl;
         // ------------- CHECK CONVERGENCE --------------
         if (MMAmain.getKKT() < tol)
         {
            printf("KKT satisfied: Norm = %f\n", MMAmain.getKKT());
            break;
         }
         // -------------- RESTART if desired ------------
         if (iter % restart == 0)
         {
            MMAmain.outputRestart(xval, iter, name);
         }
      //}
      iter++;
   }
   mma.close();

   return 0;

   //MPI_Finalize();
}


void Rosenbrock(double* xval, double a, double b, double* const fval,
                double* const dfdx,
                double* const gx, double* const dgdx)
{
   // f = b*(y-x^2)^2 + (a-x)^2
   *fval = b * (xval[1] - xval[0]*xval[0])*(xval[1] - xval[0]*xval[0]) +
           (a - xval[0])*(a - xval[0]);

   dfdx[0] = -4.0*b*xval[0]*(xval[1] - xval[0]*xval[0]) - 2.0*(a - xval[0]);
   dfdx[1] = 2.0*b*(xval[1] - xval[0]*xval[0]);

   gx[0] = (xval[0] - 1.0)*(xval[0] - 1.0) + (xval[1] - 1.0)*(xval[1] - 1.0) - 4.0;
   //gx[1] = xval[0] - 2;
   //gx[2] = xval[1] - 1;

   dgdx[0] = 2.0*(xval[0] - 1.0);
   dgdx[1] = 2.0*(xval[1] - 1.0);
   //dgdx[2] = 1.0;
   //dgdx[3] = 0.0;
   //dgdx[4] = 0.0;
   //dgdx[5] = 1.0;
}


void MultiVar_Rosenbrock(double* xval, double a, double b, int nVar, int nCon, double* const fval,
                double* const dfdx, double* const gx, double* const dgdx)
{
   
   printf("\n New run \n\n");
   *fval = 0.0;
   dfdx[0] = -4.0 * b * xval[0] * (xval[1] - xval[0] * xval[0]) - 2.0 * (a - xval[0]);
   *fval = b * (xval[1] - xval[0] * xval[0]) * (xval[1] - xval[0] * xval[0]) + (a - xval[0]) * (a - xval[0]);
   printf("xval = %f\n", xval[0]);
   for (int i = 1; i < (nVar - 1); i++) 
   {
      *fval += b * (xval[i + 1] - xval[i] * xval[i]) * (xval[i + 1] - xval[i] * xval[i]) + (a - xval[i]) * (a - xval[i]);
      dfdx[i] = -4.0 * b * xval[i] * (xval[i + 1] - xval[i] * xval[i]) - 2.0 * (a - xval[i]) + 2.0 * b * (xval[i] - xval[i - 1] * xval[i - 1]);
   }
   dfdx[nVar - 1] = 2.0 * b * (xval[nVar - 1] - xval[nVar - 2] * xval[nVar - 2]);

    // Compute gx and dgdx
   for (int i = 0; i < nCon; i++)
   {
      gx[i] = 0.0;
      for (int j = 0; j < nVar; j++) 
      {
         gx[i] += (xval[j] - 1.0) * (xval[j] - 1.0);
         dgdx[i * nVar + j] = 2.0 * (xval[j] - 1.0);
         printf("dgdx[%d] = %f\n", i * nVar + j, dgdx[i * nVar + j]);
      }
      gx[i] -= 4.0;
      printf("gx[%d] = %f\n", i, gx[i]);

   }
}
