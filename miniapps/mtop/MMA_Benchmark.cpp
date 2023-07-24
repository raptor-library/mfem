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
   int nVar = 2;
   int nCon = 3;
   double epsimin = 0.0000001;
   int sx = 1;
   double* xval = new double[nVar];
   double* x = new double[nVar];

   double* fval = new double[1];
   double* dfdx = new double[nVar];
   double* gx = new double[nCon];
   double* dgdx = new double[nVar * nCon];
   double* xmin = new double[nVar];
   double* xmax = new double[nVar];
   // Simulation parameters
   int iter = 0;
   int maxiter = 100;
   int restart = 10;
   double norm2 = 0.0;
   double normInf = 0.0;
   double kktnorm = 0.0;
   double kkttol = 0.0;
   double* low = new double[nVar];
   double* upp = new double[nVar];
   double* xo1 = new double[nVar];
   double* xo2 = new double[nVar];
   double* c = new double[nCon];
   double* d = new double[nCon];
   double* a = new double[nCon];
   // Results
   double* xmma = new double[nVar];
   double* ymma = new double[nCon];
   double* zmma = new double[nCon];
   double* lam = new double[nCon];
   double* xsi = new double[nVar];
   double* eta = new double[nVar];
   double* mu = new double[nCon];
   double zet;
   double* s = new double[nCon];
   //Initialize
   std::ofstream mma;
   mma.open("mma.dat");

   for (int i = 0; i < nVar; i++)
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
   c[2] = 1000.0;
   d[0] = 1.0;
   d[1] = 1.0;
   d[2] = 1.0;
   double a0 = 1.0;
   a[0] = 0.0;
   a[1] = 0.0;
   a[2] = 0.0;

   if (iter < 0.5)
   {
      Rosenbrock(xval, aR, bR, fval, dfdx, gx, dgdx);
   }


   kktnorm = kkttol + 10;
   while (kktnorm > kkttol && iter < maxiter)
   {
      //print results
      for (int i = 0; i < nVar; i++)
      {
         mma << xval[i] << "\n";
      }
      for (int i = 0; i < nVar; i++)
      {
         mma << upp[i] << "\n";
      }
      for (int i = 0; i < nVar; i++)
      {
         mma << low[i] << "\n";
      }

      iter++;
      // Run MMA
      MMA MMAmain(nVar,nCon,xval,sx);
      MMAmain.mmasub(nVar, nCon, iter, xval, xmin, xmax, xo1, xo2, fval, dfdx, gx,
                     dgdx,
                     low, upp, a0, a, c, d, xmma, ymma, zmma, lam, xsi, eta, mu, zet, s);
      // Update design variables
      //printf("iter = %d, xval = %f, %f, fval = %f\n", iter, xmma[0], xmma[1], *fval);
      for (int i = 0; i < nVar; i++)
      {
         xo2[i] = xo1[i];
         xo1[i] = xval[i];
         xval[i] = xmma[i];
      }
      // Compute objective and constraints
      Rosenbrock(xval, aR, bR, fval, dfdx, gx, dgdx);
      // Compute KKT residual
      MMAmain.kktcheck(nCon, nVar, xmma, ymma, *zmma, lam, xsi, eta, mu, zet, s, xmin,
                       xmax, dfdx, gx, dgdx, a0, a, c, d, &kktnorm);
      //printf("kktnorm = %f\n", kktnorm);
   }
   printf("xval = %f, %f\n", xval[0], xval[1]);
   printf("kktnorm = %f\n", kktnorm);

   delete[] xval;
   delete[] x;
   delete[] fval;
   delete[] dfdx;
   delete[] gx;
   delete[] dgdx;
   delete[] xmin;
   delete[] xmax;
   delete[] low;
   delete[] upp;
   delete[] xo1;
   delete[] xo2;
   delete[] c;
   delete[] d;
   delete[] a;
   delete[] xmma;
   delete[] ymma;
   delete[] zmma;
   delete[] lam;
   delete[] xsi;
   delete[] eta;
   delete[] mu;
   delete[] s;

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

   gx[0] = (xval[0] - 1.0)*(xval[0] - 1.0) + (xval[1] - 1.0)*(xval[1] - 1.0) - 4.0;
   gx[1] = xval[0] - 2;
   gx[2] = xval[1] - 1;

   dgdx[0] = 2.0*(xval[0] - 1.0);
   dgdx[1] = 2.0*(xval[1] - 1.0);
   dgdx[2] = 1.0;
   dgdx[3] = 0.0;
   dgdx[4] = 0.0;
   dgdx[5] = 1.0;
}
