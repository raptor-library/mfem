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

using namespace mfem;
using namespace mma;

void Rosenbrock(double* xval, double a, double b, double* fval, double* dfdx,
                double* gx, double* dgdx);

int main(int argc, char *argv[])
{
   //Rosenbrock parameters
   double aR = 1.0;
   double bR = 100.0;
   // Objective parameters
   int n = 2;
   int m = 2;
   double epsimin = 0.0000001;
   int sx = 1;
   double* xval = new double[n];
   double* x = new double[n];

   double* fval = new double;
   double* dfdx = new double[n];
   double* gx = new double[m];
   double* dgdx = new double[n*m];
   double* xmin = new double[n];
   double* xmax = new double[n];
   // Simulation parameters
   int iter = 0;
   int maxiter = 5;
   int restart = 10;
   double norm2 = 0.0;
   double normInf = 0.0;
   double kktnorm;
   double kkttol = 0.0;
   double* low = new double[n];
   double* upp = new double[n];
   double* xo1 = new double[n];
   double* xo2 = new double[n];
   double* c = new double[n];
   double* d = new double[n];
   double* a = new double[n];
   // Results
   double* xmma = new double[n];
   double* ymma = new double[m];
   double* zmma = new double[m];
   double* lam = new double[m];
   double* xsi = new double[n];
   double* eta = new double[n];
   double* mu = new double[m];
   double zet;
   double* s = new double[m];
   //Initialize
   for (int i = 0; i < n; i++)
   {
      xval[i] = 0.0;
      x[i] = 0.0;
      xo1[i] = 0.0;
      xo2[i] = 0.0;
      xmin[i] = -1.0;
      xmax[i] = 2.0;
      low[i] = xmin[i];
      upp[i] = xmax[i];
   }
   c[0] = 1000.0;
   c[1] = 1000.0;
   d[0] = 1.0;
   d[1] = 1.0;
   double a0 = 0.0;
   a[0] = 0.0;
   a[1] = 0.0;

   if (iter < 0.5)
   {
      Rosenbrock(xval, aR, bR, fval, dfdx, gx, dgdx);
   }
   kktnorm = kkttol + 10;
   while (kktnorm > kkttol && iter < maxiter)
   {
      iter++;
      // Run MMA
      MMA MMAmain(n,m,xval,sx);
      MMAmain.mmasub(n, m, iter, xval, xmin, xmax, xo1, xo2, *fval, dfdx, gx, dgdx,
                     low, upp, a0, a, c, d, xmma, ymma, zmma, lam, xsi, eta, mu, zet, s);
      // Update design variables
      printf("iter = %d, xval = %f, %f, fval = %f\n", iter, xval[0], xval[1], *fval);
      for (int i = 0; i < n; i++)
      {
         xo2[i] = xo1[i];
         xo1[i] = x[i];
         x[i] = xval[i];
      }
      // Compute objective and constraints
      Rosenbrock(xval, aR, bR, fval, dfdx, gx, dgdx);
      // Compute KKT residual
      //MMAmain.KKTresidual(xval, dfdx, gx, dgdx, xmin, xmax, &norm2, &normInf);
      // Check convergence
      if (norm2 < 1e-6 && normInf < 1e-6)
      {
         break;
      }
   }


   return 0;
}


void Rosenbrock(double* xval, double a, double b, double* fval, double* dfdx,
                double* gx, double* dgdx)
{
   // f = b*(y-x^2)^2 + (a-x)^2
   *fval = b * (xval[1] - xval[0]*xval[0])*(xval[1] - xval[0]*xval[0]) +
            (a - xval[0])*(a - xval[0]);

   dfdx[0] = -4.0*b*xval[0]*(xval[1] - xval[0]*xval[0]) - 2.0*(a - xval[0]);
   dfdx[1] = 2.0*b*(xval[1] - xval[0]*xval[0]);

   gx[0] = xval[0] + xval[1] - 1.0;
   gx[1] = 1.0 - xval[0] - xval[1];

   dgdx[0] = 1.0;
   dgdx[1] = -1.0;
}
