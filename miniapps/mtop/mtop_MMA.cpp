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

/**
 * xval: design variables
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
*/

#include "mtop_MMA.hpp"
#include <fstream>
#include <math.h>

#ifndef MFEM_USE_LAPACK
mfem::mfem_error("MMA relies on LAPACK. Please use MFEM_USE_LAPACK");
#endif

extern "C" void dgesv_(int* nLAP, int* nrhs, double* AA, int* lda, int* ipiv,
                       double* bb, int* ldb, int* info);

namespace mma
{

void MMA::Update(int iter, double* const fval, double* const dfdx, double* const gx, double* const dgdx, double* xval)
{
   for (int i = 0; i < nVar; i++)
   {
      //factor[i] = asyincr;
      p0[i] = 0.0;
      q0[i] = 0.0;
      pq0[i] = 0.0;
      xmami[i] = 0.0;
      ux1[i] = 0.0;
      xl1[i] = 0.0;
      alfa[i] = 0.0;
      beta[i] = 0.0;
   }
   for (int i = 0; i < (nCon*nVar); i++)
   {
      p[i] = 0.0;
      q[i] = 0.0;
      pq[i] = 0.0;
      P[i] = 0.0;
      Q[i] = 0.0;
      PQ[i] = 0.0;
   }
   for (int i = 0; i < nCon; i++)
   {
      b[i] = 0.0;
   }
   // Calculation of the asymptotes low and upp
   if (iter < 2.5)
   {
      for (int i = 0; i < nVar; i++)
      {
         low[i] = xval[i] - asyinit * (xmax[i] - xmin[i]);
         upp[i] = xval[i] + asyinit * (xmax[i] - xmin[i]);
      }
   }
   else
   {

      for (int i = 0; i < nVar; i++)
      {
         //Determine sign
         zz = (xval[i] - xo1[i]) * (xo1[i] - xo2[i]);
         if ( zz > 0.0)
         {
            factor[i] = asyincr;
         }
         else if ( zz < 0.0)
         {
            factor[i] = asydecr;
         }
         else
         {
            factor[i] = 1.0;
         }
         

         //Find new asymptote
         low[i] = xval[i] - factor[i] * (xo1[i] - low[i]);
         upp[i] = xval[i] + factor[i] * (upp[i] - xo1[i]);

         lowmin = xval[i] - 10.0 * (xmax[i] - xmin[i]);
         lowmax = xval[i] - 0.01 * (xmax[i] - xmin[i]);
         uppmin = xval[i] + 0.01 * (xmax[i] - xmin[i]);
         uppmax = xval[i] + 10.0 * (xmax[i] - xmin[i]);

         low[i] = std::max(low[i], lowmin);
         low[i] = std::min(low[i], lowmax);
         upp[i] = std::min(upp[i], uppmax);
         upp[i] = std::max(upp[i], uppmin);
      }
   }
   for (int i = 0; i < nVar; i++)
   {
      // Calculation of bounds alfa and beta according to:
      // alfa = max{xmin, low + 0.1(xval-low), xval-0.5(xmax-xmin)}
      // beta = min{xmax, upp - 0.1(upp-xval), xval+0.5(xmax-xmin)}
      
      alfa[i] = std::max(std::max(low[i] + albefa * (xval[i] - low[i]),
                                  xval[i] - move * (xmax[i] - xmin[i])), xmin[i]);
      beta[i] = std::min(std::min(upp[i] - albefa * (upp[i] - xval[i]),
                                  xval[i] + move * (xmax[i] - xmin[i])), xmax[i]);
      xmami[i] = std::max(xmax[i] - xmin[i], xmamieps);

      // Calculations of p0, q0, P, Q, and b
      ux1[i] = upp[i] - xval[i];
      if (std::fabs(ux1[i]) <= machineEpsilon || ux1[i] == 0.0)
      {
         ux1[i] = machineEpsilon;
      }
      xl1[i] = xval[i] - low[i];
      if (std::fabs(xl1[i]) <= machineEpsilon || xl1[i] == 0.0)
      {
         xl1[i] = machineEpsilon;
      }
      p0[i] = std::max(dfdx[i], 0.0);
      q0[i] = std::max(-dfdx[i], 0.0);
      pq0[i] = 0.001 * (p0[i] + q0[i]) + raa0 / xmami[i];
      p0[i] += pq0[i];
      q0[i] += pq0[i];
      p0[i] *= ux1[i] * ux1[i];
      q0[i] *= xl1[i] * xl1[i];
   }
   // P = max(dgdx,0)
   // Q = max(-dgdx,0)
   // PQ = 0.001(P+Q) + raa0/xmami
   // P = P + PQ
   // Q = Q + PQ
   for (int i = 0; i < nCon; i++)
   {
      for (int j = 0; j < nVar; j++)
      {
         p[i * nVar + j] = std::max(dgdx[i * nVar + j], 0.0);
         q[i * nVar + j] = std::max(-1*dgdx[i * nVar + j], 0.0);
         pq[i * nVar + j] = 0.001 * (p[i * nVar + j] + q[i * nVar + j]) + raa0 /
                            xmami[j];
         p[i * nVar + j] += pq[i * nVar + j];
         q[i * nVar + j] += pq[i * nVar + j];
         // P = P * spdiags(ux2,0,n,n)
         // Q = Q * spdiags(xl2,0,n,n)
         P[i * nVar + j] = p[i * nVar + j] * ux1[j] * ux1[j];
         Q[i * nVar + j] = q[i * nVar + j] * xl1[j] * xl1[j];
         // b = P/ux1 + Q/xl1 - gx
         b[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
      }
      b[i] -= gx[i];
   }

   //subsolv();
   mSubProblem->Perform();

   // Update design variables
   for (int i = 0; i < nVar; i++)
   {
      xo2[i] = xo1[i];
      xo1[i] = xval[i];
      xval[i] = x[i];
   }
   //mSubProblem->VariableA();
   
   //mmasub(iter, fval, dfdx, gx, dgdx, xval);
   kktcheck(y, dfdx, gx, dgdx, xval);
}

void MMA::mmasub(int iter, double* const fval, double* const dfdx, double* const gx, double* const dgdx, double* xval)
{
   // for (int i = 0; i < nVar; i++)
   // {
   //    //factor[i] = asyincr;
   //    p0[i] = 0.0;
   //    q0[i] = 0.0;
   //    pq0[i] = 0.0;
   //    xmami[i] = 0.0;
   //    ux1[i] = 0.0;
   //    xl1[i] = 0.0;
   //    alfa[i] = 0.0;
   //    beta[i] = 0.0;
   // }
   // for (int i = 0; i < (nCon*nVar); i++)
   // {
   //    p[i] = 0.0;
   //    q[i] = 0.0;
   //    pq[i] = 0.0;
   //    P[i] = 0.0;
   //    Q[i] = 0.0;
   //    PQ[i] = 0.0;
   // }
   // for (int i = 0; i < nCon; i++)
   // {
   //    b[i] = 0.0;
   // }
   // // Calculation of the asymptotes low and upp
   // if (iter < 2.5)
   // {
   //    for (int i = 0; i < nVar; i++)
   //    {
   //       low[i] = xval[i] - asyinit * (xmax[i] - xmin[i]);
   //       upp[i] = xval[i] + asyinit * (xmax[i] - xmin[i]);
   //    }
   // }
   // else
   // {

   //    for (int i = 0; i < nVar; i++)
   //    {
   //       //Determine sign
   //       zz = (xval[i] - xo1[i]) * (xo1[i] - xo2[i]);
   //       if ( zz > 0.0)
   //       {
   //          factor[i] = asyincr;
   //       }
   //       else if ( zz < 0.0)
   //       {
   //          factor[i] = asydecr;
   //       }
   //       else
   //       {
   //          factor[i] = 1.0;
   //       }
         

   //       //Find new asymptote
   //       low[i] = xval[i] - factor[i] * (xo1[i] - low[i]);
   //       upp[i] = xval[i] + factor[i] * (upp[i] - xo1[i]);

   //       lowmin = xval[i] - 10.0 * (xmax[i] - xmin[i]);
   //       lowmax = xval[i] - 0.01 * (xmax[i] - xmin[i]);
   //       uppmin = xval[i] + 0.01 * (xmax[i] - xmin[i]);
   //       uppmax = xval[i] + 10.0 * (xmax[i] - xmin[i]);

   //       low[i] = std::max(low[i], lowmin);
   //       low[i] = std::min(low[i], lowmax);
   //       upp[i] = std::min(upp[i], uppmax);
   //       upp[i] = std::max(upp[i], uppmin);
   //    }
   // }
   // for (int i = 0; i < nVar; i++)
   // {
   //    // Calculation of bounds alfa and beta according to:
   //    // alfa = max{xmin, low + 0.1(xval-low), xval-0.5(xmax-xmin)}
   //    // beta = min{xmax, upp - 0.1(upp-xval), xval+0.5(xmax-xmin)}
      
   //    alfa[i] = std::max(std::max(low[i] + albefa * (xval[i] - low[i]),
   //                                xval[i] - move * (xmax[i] - xmin[i])), xmin[i]);
   //    beta[i] = std::min(std::min(upp[i] - albefa * (upp[i] - xval[i]),
   //                                xval[i] + move * (xmax[i] - xmin[i])), xmax[i]);
   //    xmami[i] = std::max(xmax[i] - xmin[i], xmamieps);

   //    // Calculations of p0, q0, P, Q, and b
   //    ux1[i] = upp[i] - xval[i];
   //    if (std::fabs(ux1[i]) <= machineEpsilon || ux1[i] == 0.0)
   //    {
   //       ux1[i] = machineEpsilon;
   //    }
   //    xl1[i] = xval[i] - low[i];
   //    if (std::fabs(xl1[i]) <= machineEpsilon || xl1[i] == 0.0)
   //    {
   //       xl1[i] = machineEpsilon;
   //    }
   //    p0[i] = std::max(dfdx[i], 0.0);
   //    q0[i] = std::max(-dfdx[i], 0.0);
   //    pq0[i] = 0.001 * (p0[i] + q0[i]) + raa0 / xmami[i];
   //    p0[i] += pq0[i];
   //    q0[i] += pq0[i];
   //    p0[i] *= ux1[i] * ux1[i];
   //    q0[i] *= xl1[i] * xl1[i];
   // }
   // // P = max(dgdx,0)
   // // Q = max(-dgdx,0)
   // // PQ = 0.001(P+Q) + raa0/xmami
   // // P = P + PQ
   // // Q = Q + PQ
   // for (int i = 0; i < nCon; i++)
   // {
   //    for (int j = 0; j < nVar; j++)
   //    {
   //       p[i * nVar + j] = std::max(dgdx[i * nVar + j], 0.0);
   //       q[i * nVar + j] = std::max(-1*dgdx[i * nVar + j], 0.0);
   //       pq[i * nVar + j] = 0.001 * (p[i * nVar + j] + q[i * nVar + j]) + raa0 /
   //                          xmami[j];
   //       p[i * nVar + j] += pq[i * nVar + j];
   //       q[i * nVar + j] += pq[i * nVar + j];
   //       // P = P * spdiags(ux2,0,n,n)
   //       // Q = Q * spdiags(xl2,0,n,n)
   //       P[i * nVar + j] = p[i * nVar + j] * ux1[j] * ux1[j];
   //       Q[i * nVar + j] = q[i * nVar + j] * xl1[j] * xl1[j];
   //       // b = P/ux1 + Q/xl1 - gx
   //       b[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
   //    }
   //    b[i] -= gx[i];
   // }

   // subsolv();

   // // Update design variables
   // for (int i = 0; i < nVar; i++)
   // {
   //    xo2[i] = xo1[i];
   //    xo1[i] = xval[i];
   //    xval[i] = x[i];
   // }
}

/**
 * Subproblem functionality
*/

void MMA::SubProblemClassic::Perform()
{
   ittt = 0;
   itto = 0;
   epsi = 1;
   itera = 0;
   z = 1;
   zet = 1;

   for (int i = 0; i < nVar; i++)
   {
      epsvecn[i] = epsi;
      x[i] = 0.5 * (alfa[i] + beta[i]);
      xsi[i] = 1.0/(x[i] - alfa[i]);
      xsi[i] = std::max(xsi[i], 1.0);
      eta[i] = 1.0/(beta[i] - x[i]);
      eta[i] = std::max(eta[i], 1.0);
   }
   for (int i = 0; i < nCon; i++)
   {
      epsvecm[i] = epsi;
      y[i] = 1.0;
      lam[i] = 1.0;
      mu[i] = std::max(1.0, 0.5 * c[i]);
      s[i] = 1.0;
   }

   //for (int Ik = 0; true; Ik++)
   while (epsi > epsimin)
   {

      rez = a0 - zet;
      for (int i = 0; i < nVar; i++)
      {
         epsvecn[i] = epsi;
         ux1[i] = upp[i] - x[i];
         if (ux1[i] == 0)
         {
            ux1[i] = machineEpsilon;
         }
         else
         {
            while (std::fabs(ux1[i]) <= machineEpsilon)
            {
               ux1[i] *= 10;
            }
         }
         xl1[i] = x[i] - low[i];
         if (xl1[i] == 0)
         {
            xl1[i] = machineEpsilon;
         }
         else
         {
            while (std::fabs(xl1[i]) <= machineEpsilon)
            {
               xl1[i] *= 10;
            }
         }
         // plam = P' * lam, qlam = Q' * lam
         plam[i] = p0[i];
         qlam[i] = q0[i];
         for (int j = 0; j < nCon; j++)
         {
            plam[i] += P[j * nVar + i] * lam[j];
            qlam[i] += Q[j * nVar + i] * lam[j];
         }
         dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
         rex[i] = dpsidx[i] - xsi[i] + eta[i];
         rez -= a[i] * lam[i];
         rexsi[i] = xsi[i] * (x[i] - alfa[i]) - epsvecn[i];
         if (std::fabs(x[i]-alfa[i]) == 0)
         {
            rexsi[i] = xsi[i] * machineEpsilon - epsvecn[i];
         }
         reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
         if (std::fabs(beta[i] - x[i]) == 0)
         {
            reeta[i] = eta[i] * machineEpsilon - epsvecn[i];;
         }
         residu1[i] = rex[i];
         residu2[nCon + i] = rexsi[i];
         residu2[nCon + nVar + i] = reeta[i];
      }
      for (int i = 0; i < nCon; i++)
      {
         epsvecm[i] = epsi;
         // gvec = P/ux + Q/xl
         for (int j = 0; j < nVar; j++)
         {
            gvec[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
         }
         rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
         relam[i] = gvec[i] - a[i] * z - y[i] + s[i] - b[i];
         remu[i] = mu[i] * y[i] - epsvecm[i];
         res[i] = lam[i] * s[i] - epsvecm[i];
         residu2[i] = relam[i];
         residu1[nVar + i] = rey[i];
         residu2[nCon + 2 * nVar + i] = remu[i];
         residu2[2 * nVar + 2 * nCon + 1 + i] = res[i];
      }
      rezet = zet * z - epsi;
      residu1[nVar + nCon] = rez;
      residu2[2 * nVar + 2 * nCon] = rezet;

      // Concatenate the residuals
      for (int i = 0; i < (nVar + nCon + 1); i++)
      {
         residu[i] = residu1[i];
      }
      for (int i = 0; i < (2 * nVar + 3 * nCon + 1); i++)
      {
         residu[nVar + nCon + 1 + i] = residu2[i];
      }

      //Get vector product and maximum absolute value
      residunorm = 0.0;
      residumax = 0.0;
      for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
      {
         residunorm += residu[i] * residu[i];
         residumax = std::max(residumax, std::abs(residu[i]));
      }
      // Norm of the residual
      residunorm = std::sqrt(residunorm);

      ittt = 0;

      while (residumax > 0.9 * epsi && ittt < 200)
      {
         ittt++;
         for (int i = 0; i < nVar; i++)
         {
            ux1[i] = upp[i] - x[i];
            if (ux1[i] == 0)
            {
               ux1[i] = machineEpsilon;
            }
            else
            {
               while (std::fabs(ux1[i]) <= machineEpsilon)
               {
                  ux1[i] *= 10;
               }
            }
            xl1[i] = x[i] - low[i];
            if (xl1[i] == 0)
            {
               xl1[i] = machineEpsilon;
            }
            else
            {
               while (std::fabs(xl1[i]) <= machineEpsilon)
               {
                  xl1[i] *= 10;
               }
            }
            // plam = P' * lam, qlam = Q' * lam
            plam[i] = p0[i];
            qlam[i] = q0[i];
            for (int j = 0; j < nCon; j++)
            {
               plam[i] += P[j * nVar + i] * lam[j];
               qlam[i] += Q[j * nVar + i] * lam[j];
            }
            dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
            // NaN-Avoidance
            if (std::fabs(x[i] - alfa[i]) < machineEpsilon)
            {
               if (std::fabs(beta[i] - x[i]) < machineEpsilon)
               {
                  delx[i] = dpsidx[i];
                  diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / machineEpsilon + eta[i] / machineEpsilon;
               }
               else
               {
                  delx[i] = dpsidx[i] - epsvecn[i] / machineEpsilon + epsvecn[i] / (beta[i] - x[i]);
                  diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / (x[i] - alfa[i]) + eta[i] / (beta[i] - x[i]);
               }
            }
            else if (std::fabs(beta[i] - x[i]) < machineEpsilon)
            {
               delx[i] = dpsidx[i] - epsvecn[i] / (x[i] - alfa[i]) + epsvecn[i] / machineEpsilon;
               diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / (x[i] - alfa[i]) + eta[i] / machineEpsilon;
            }
            else
            {
               delx[i] = dpsidx[i] - epsvecn[i] / (x[i] - alfa[i]) + epsvecn[i] / (beta[i] - x[i]);
               diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / (x[i] - alfa[i]) + eta[i] / (beta[i] - x[i]);
            }
         }

         delz = a0 - epsi / z;
         for (int i = 0; i < nCon; i++)
         {
            gvec[i] = 0.0;
            // gvec = P/ux + Q/xl
            for (int j = 0; j < nVar; j++)
            {
               gvec[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
               Puxinv[i * nVar + j] = P[i * nVar + j] / (ux1[j] * ux1[j]);
               Qxlinv[i * nVar + j] = Q[i * nVar + j] / (xl1[j] * xl1[j]);

               GG[i * nVar + j] = Puxinv[i * nVar + j] - Qxlinv[i * nVar + j];
            }

            dely[i] = c[i] + d[i] * y[i] - lam[i] - epsvecm[i] / y[i];
            delz -= a[i] * lam[i];
            dellam[i] = gvec[i] - a[i] * z - y[i] - b[i] + epsvecm[i] / lam[i];
            diagy[i] = d[i] + mu[i] / y[i];
            //diagyinv[i] = 1.0 / diagy[i];
            diaglam[i] = s[i] / lam[i];
            diaglamyi[i] = diaglam[i] + 1.0 / diagy[i];
         }

         if (nCon < nVar)
         {
            // blam = dellam + dely./diagy - GG*(delx./diagx);
            // bb = [blam; delz];
            for (int j = 0; j < nCon; j++)
            {
               sum = 0.0;
               for (int i = 0; i < nVar; i++)
               {
                  sum += GG[j * nVar + i] * (delx[i] / diagx[i]);
               }
               blam[j] = dellam[j] + dely[j] / diagy[j] - sum;
               bb1[j] = blam[j];
            }
            bb1[nCon] = delz;            

            // Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
            for (int i = 0; i < nCon; i++)
            {
               // Axx = GG*spdiags(diagxinv,0,n,n);
               for (int k = 0; k < nVar; k++)
               {
                  Axx[i * nVar + k] = GG[k * nCon + i] / diagx[k];
               }
            }
            // Alam = spdiags(diaglamyi,0,m,m) + Axx*GG';
            for (int i = 0; i < nCon; i++)
            {
               for (int j = 0; j < nCon; j++)
               {
                  Alam[i * nCon + j] = 0.0;
                  for (int k = 0; k < nVar; k++)
                  {
                     Alam[i * nCon + j] += Axx[i * nVar + k] * GG[j * nVar + k];
                  }
                  if (i == j)
                  {
                     Alam[i * nCon + j] += diaglamyi[i];
                  }
               }
            }
            // AA1 = [Alam     a
            //       a'    -zet/z];
            for (int i = 0; i < nCon; i++)
            {
               for (int j = 0; j < nCon; j++)
               {
                  AA1[i * nCon + j] = Alam[i * nCon + j];
               }
               AA1[nCon * nCon + i] = a[i];
               AA1[nCon * nCon + nCon + i] = a[i];
            }
            AA1[nCon * nCon + nCon + nCon] = -zet / z;                              

            // ----------------------------------------------------------------------------
            //solut = AA1\bb1 --> solve linear system of equations using LAPACK
            int info;
            int nLAP = nCon + 1;
            int nrhs = 1;
            int lda = nLAP;
            int ldb = nLAP;
            int* ipiv = new int[nLAP];
            dgesv_(&nLAP, &nrhs, AA1, &lda, ipiv, bb1, &ldb, &info);
            delete[] ipiv;
            for (int i = 0; i < (nCon + 1); i++)
            {
               solut[i] = bb1[i];
            }
            // ----------------------------------------------------------------------------
            //dlam = solut(1:nCon);
            for (int i = 0; i < nCon; i++)
            {
               dlam[i] = solut[i];
            }
            dz = solut[nCon];
            //dx = -(GG'*dlam)./diagx - delx./diagx;
            for (int i = 0; i < nVar; i++)
            {
               sum = 0.0;
               for (int j = 0; j < nCon; j++)
               {
                  sum += GG[j * nVar + i] * dlam[j];
               }
               dx[i] = -sum / diagx[i] - delx[i] / diagx[i];           
            }
         }
         else
         {
            sum = 0.0;
            for (int i = 0; i < nCon; i++)
            {
               //diaglamyiinv[i] = 1.0 / diaglamyi[i];
               dellamyi[i] = dellam[i] + dely[i] / diagy[i];
               // azz = zet/z + a'*(a./diaglamyi)
               sum += a[i] * (a[i] / diaglamyi[i]);
            }
            azz = zet / z + sum;
            // Axx = spdiags(diagx,0,nVar,nVar) + GG'*spdiags(diaglamyiinv,0,nCon,nCon)*GG;
            // AA = [Axx      axz
            //       axz'     azz];
            for (int i = 0; i < nVar; i++)
            {
               // Axx =  GG'*spdiags(diaglamyiinv,0,nCon,nCon);
               for (int k = 0; k < nCon; k++)
               {
                  Axx[i * nCon + k] = GG[k * nVar + i] / diaglamyi[k];
               }
               sum = 0.0;
               // axz = -GG'*(a./diaglamyi)
               for (int j = 0; j < nCon; j++)
               {
                  sum -= GG[j * nVar + i] * (a[j] / diaglamyi[j]);
               }
               axz[i] = sum;
            }
            //Assemble matrix AA
            for (int i = 0; i < (nVar + 1); i++)
            {
               for (int j = 0; j < (nVar + 1); j++)
               {
                  // AA = [Axx  .
                  //       .    .]
                  AA[i * (nVar + 1) + j] = 0.0;
                  if (i < nVar && j < nVar)
                  {
                     // Axx = Axx*GG
                     for (int k = 0; k < nCon; k++)
                     {
                        AA[i * (nVar + 1) + j] += Axx[i * nCon + k] * GG[k * nVar + j];
                     }
                     // Axx = Axx + spdiags(diagx,0,nVar,nVar)
                     if (i == j)
                     {
                        AA[i * (nVar + 1) + j] += diagx[j];
                     }
                  }
                  // AA = [Axx  axz
                  //       axz' azz]
                  else if (i < nVar && j == nVar)
                  {
                     AA[i * (nVar + 1) + j] = axz[i];
                  }
                  else if (i == nVar && j < nVar)
                  {
                     AA[i * (nVar + 1) + j] = axz[j];
                  }
                  else
                  {
                     AA[i * (nVar + 1) + j] = azz;
                  }
               }
            }
            // bb = [-bx'; -bz]
            // bx = delx - GG'*(dellamyi./diaglamyi)
            // bz = delz - a'*(dellamyi./diaglamyi)
            for (int i = 0; i < nVar; i++)
            {
               sum = 0.0;
               for (int j = 0; j < nCon; j++)
               {
                  sum += GG[j * nVar + i] * (dellamyi[j] / diaglamyi[j]);
               }
               bb[i] = -(delx[i] + sum);
            }
            sum = 0.0;
            for (int i = 0; i < nCon; i++)
            {
               sum += a[i] * (dellamyi[i] / diaglamyi[i]);
            }
            bb[nVar] = -(delz - sum);
            // ----------------------------------------------------------------------------
            //solut = AA\bb --> solve linear system of equations using LAPACK
            int info;
            int nLAP = nVar + 1;
            int nrhs = 1;
            int lda = nLAP;
            int ldb = nLAP;
            int* ipiv = new int[nLAP];
            dgesv_(&nLAP, &nrhs, AA, &lda, ipiv, bb, &ldb, &info);
            delete[] ipiv;
            for (int i = 0; i < (nVar + 1); i++)
            {
               solut[i] = bb[i];
            }
            // ----------------------------------------------------------------------------
            //dx = solut(1:nVar);
            for (int i = 0; i < nVar; i++)
            {
               dx[i] = solut[i];
            }
            dz = solut[nVar];
            //dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
            for (int i = 0; i < nCon; i++)
            {
               sum = 0.0;
               for (int j = 0; j < nVar; j++)
               {
                  sum += GG[i * nVar + j] * dx[j];
               }
               dlam[i] = sum / diaglamyi[i] - dz * (a[i] / diaglamyi[i]) + dellamyi[i] / diaglamyi[i];
            }
         }

         for (int i = 0; i < nCon; i++)
         {
            dy[i] = -dely[i] / diagy[i] + dlam[i] / diagy[i];
            dmu[i] = -mu[i] + epsvecm[i] / y[i] - (mu[i] * dy[i]) / y[i];
            ds[i] = -s[i] + epsvecm[i] / lam[i] - (s[i] * dlam[i]) / lam[i];
            // xx = [y z lam xsi eta mu zet s]
            // dxx = [dy dz dlam dxsi deta dmu dzet ds]
            xx[i] = y[i];
            xx[nCon + 1 + i] = lam[i];
            xx[2 * nCon + 1 + 2 * nVar + i] = mu[i];
            xx[3 * nCon + 2 * nVar + 2 + i] = s[i];
            dxx[i] = dy[i];
            dxx[nCon + 1 + i] = dlam[i];
            dxx[2 * nCon + 1 + 2 * nVar + i] = dmu[i];
            dxx[3 * nCon + 2 * nVar + 2 + i] = ds[i];
         }
         xx[nCon] = z;
         xx[3 * nCon + 2 * nVar + 1] = zet;
         dxx[nCon] = dz;
         dxx[3 * nCon + 2 * nVar + 1] = dzet;
         for (int i = 0; i < nVar; i++)
         {
            // NaN-Avoidance
            if (std::fabs(x[i] - alfa[i]) < machineEpsilon)
            {
               if (std::fabs(beta[i] - x[i]) < machineEpsilon)
               {
                  dxsi[i] = -xsi[i] + epsvecn[i] / machineEpsilon - (xsi[i] * dx[i]) / machineEpsilon;
                  deta[i] = -eta[i] + epsvecn[i] / machineEpsilon + (eta[i] * dx[i]) / machineEpsilon;
               }
               else
               {
                  dxsi[i] = -xsi[i] + epsvecn[i] / machineEpsilon - (xsi[i] * dx[i]) / machineEpsilon;
                  deta[i] = -eta[i] + epsvecn[i] / (beta[i] - x[i]) + (eta[i] * dx[i]) / (beta[i] - x[i]);
               }
            }
            else if (std::fabs(beta[i] - x[i]) < machineEpsilon)
            {
               dxsi[i] = -xsi[i] + epsvecn[i] / (x[i] - alfa[i]) - (xsi[i] * dx[i]) / (x[i] - alfa[i]);
               deta[i] = -eta[i] + epsvecn[i] / machineEpsilon + (eta[i] * dx[i]) / machineEpsilon;
            }
            else
            {
               dxsi[i] = -xsi[i] + epsvecn[i] / (x[i] - alfa[i]) - (xsi[i] * dx[i]) / (x[i] - alfa[i]);
               deta[i] = -eta[i] + epsvecn[i] / (beta[i] - x[i]) + (eta[i] * dx[i]) / (beta[i] - x[i]);
            }
            xx[nCon + 1 + nCon + i] = xsi[i];
            xx[nCon + 1 + nCon + nVar + i] = eta[i];
            dxx[nCon + 1 + nCon + i] = dxsi[i];
            dxx[nCon + 1 + nCon + nVar + i] = deta[i];
         }
         dzet = -zet + epsi / z - zet * dz / z;

         stmxx = 0.0;
         for (int i = 0; i < (4 * nCon + 2 * nVar + 2); i++)
         {
            stepxx[i] = -1.01*dxx[i] / xx[i];
            if (stepxx[i] > stmxx)
            {
               stmxx = stepxx[i];
            }
         }
         stmalfa = 0.0;
         stmbeta = 0.0;
         for (int i = 0; i < nVar; i++)
         {
            
            //NaN-Avoidance
            if (std::fabs(x[i] - alfa[i]) < machineEpsilon)
            {
               stepalfa[i] = -1.01*dx[i] / machineEpsilon;
            }
            else 
            {
               stepalfa[i] = -1.01*dx[i] / (x[i] - alfa[i]);
            }
            if (std::fabs(beta[i] - x[i]) < machineEpsilon)
            {
               stepbeta[i] = 1.01*dx[i] / machineEpsilon;
            }
            else
            {
               stepbeta[i] = 1.01*dx[i] / (beta[i] - x[i]);
            }
            // --------------
            if (stepalfa[i] > stmalfa)
            {
               stmalfa = stepalfa[i];
            }
            if (stepbeta[i] > stmbeta)
            {
               stmbeta = stepbeta[i];
            }
         }
         stmalbe = std::max(stmalfa, stmbeta);
         stmalbexx = std::max(stmalbe, stmxx);
         stminv = std::max(stmalbexx, 1.0);
         steg = 1.0 / stminv;

         for (int i = 0; i < nVar; i++)
         {
            xold[i] = x[i];
            xsiold[i] = xsi[i];
            etaold[i] = eta[i];
         }
         for (int i = 0; i < nCon; i++)
         {
            yold[i] = y[i];
            lamold[i] = lam[i];
            muold[i] = mu[i];
            sold[i] = s[i];
         }
         zold = z;
         zetold = zet;

         itto = 0;
         resinew = 2.0 * residunorm;
         while (resinew > residunorm && itto < 50)
         {
            itto++;

            for (int i = 0; i < nCon; ++i)
            {
               y[i] = yold[i] + steg * dy[i];
               if (std::fabs(y[i]) < machineEpsilon)
               {
                  y[i] = machineEpsilon;
               }
               else
               {
                  while (std::fabs(y[i]) < machineEpsilon)
                  {
                     y[i] *= 10;
                  }
               }
               lam[i] = lamold[i] + steg * dlam[i];
               if (std::fabs(lam[i]) < machineEpsilon)
               {
                  lam[i] = machineEpsilon;
               }
               else
               {
                  while (std::fabs(lam[i]) < machineEpsilon)
                  {
                     lam[i] *= 10;
                  }
               }
               mu[i] = muold[i] + steg * dmu[i];
               s[i] = sold[i] + steg * ds[i];
            }

            for (int i = 0; i < nVar; ++i)
            {
               x[i] = xold[i] + steg * dx[i];
               xsi[i] = xsiold[i] + steg * dxsi[i];
               eta[i] = etaold[i] + steg * deta[i];
               if (std::isnan(eta[i]))
               {
                  eta[i] = machineEpsilon;
               }
               ux1[i] = upp[i] - x[i];
               if (ux1[i] == 0)
               {
                  ux1[i] = machineEpsilon;
               }
               else
               {
                  while (std::fabs(ux1[i]) <= machineEpsilon)
                  {
                     ux1[i] *= 10;
                  }
               }
               xl1[i] = x[i] - low[i];
               if (xl1[i] == 0)
               {
                  xl1[i] = machineEpsilon;
               }
               else
               {
                  while (std::fabs(xl1[i]) <= machineEpsilon)
                  {
                     xl1[i] *= 10;
                  }
               }
               // plam & qlam
               plam[i] = p0[i];
               qlam[i] = q0[i];
               for (int j = 0; j < nCon; j++)
               {
                  plam[i] += P[j * nVar + i] * lam[j];
                  qlam[i] += Q[j * nVar + i] * lam[j];
               }
               dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
               rex[i] = dpsidx[i] - xsi[i] + eta[i];
               rez -= a[i] * lam[i];
               rexsi[i] = xsi[i] * (x[i] - alfa[i]) - epsvecn[i];              
               if (std::fabs(x[i] - alfa[i]) < machineEpsilon || std::isnan(rexsi[i]))
               {
                  rexsi[i] = machineEpsilon;
               }
               reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
               if (std::fabs(beta[i] - x[i]) < machineEpsilon || std::isnan(reeta[i]))
               {
                  reeta[i] = machineEpsilon;
               }
               residu1[i] = rex[i];
               residu2[nCon + i] = rexsi[i];
               residu2[nCon + nVar + i] = reeta[i];
            }
            z = zold + steg * dz;
            if (z == 0)
            {
               z = machineEpsilon;
            }
            else
            {
               while (std::fabs(z) <= machineEpsilon)
               {
                  z *= 10;
               }
            }
            zet = zetold + steg * dzet;

            // gvec = P/ux + Q/xl
            for (int i = 0; i < nCon; i++)
            {
               gvec[i] = 0.0;
               for (int j = 0; j < nVar; j++)
               {
                  gvec[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
               }
               rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
               relam[i] = gvec[i] - a[i] * z - y[i] + s[i] - b[i];
               remu[i] = mu[i] * y[i] - epsvecm[i];
               res[i] = lam[i] * s[i] - epsvecm[i];

               residu2[i] = relam[i];
               residu1[nVar + i] = rey[i];
               residu2[nCon + 2 * nVar + i] = remu[i];
               residu2[2 * nVar + 2 * nCon + 1 + i] = res[i];
            }
            rez = a0 - zet;
            rezet = zet * z - epsi;
            residu1[nVar + nCon] = rez;
            residu2[2 * nVar + 2 * nCon] = rezet;

            // Concatenate the residuals
            for (int i = 0; i < (nVar + nCon + 1); i++)
            {
               residu[i] = residu1[i];
            }
            for (int i = 0; i < (2 * nVar + 3 * nCon + 1); i++)
            {
               residu[nVar + nCon + 1 + i] = residu2[i];
            }
            // New residual
            sum = 0.0;
            for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
            {
               sum += residu[i] * residu[i];
            }
            // Norm of the residual
            resinew = std::sqrt(sum);
            steg = steg / 2.0;
         }
         residunorm = resinew;
         residumax = 0.0;
         for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
         {
            residumax = std::max(residumax, std::abs(residu[i]));
         }
         steg = steg * 2.0;
      }


      if (ittt > 198)
      {
         printf("Warning: Maximum number of iterations reached in subsolv.\n");
      }
      epsi = 0.1 * epsi;
   }
   //results.close();
}


/**
 * This function solves the MMA subproblem:
 * minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
 *         + SUM[ ci*yi + 0.5*di*(yi)^2 ],
 *
 * subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
 *           alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
 *
 * Input: m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b.
 * Output: x, y, z, slack variables and Lagrange multiplers.
*/
void MMA::subsolv()
{
   //std::ofstream results;
   //results.open("sub.dat", std::ios::app);

   // ittt = 0;
   // itto = 0;
   // epsi = 1;
   // itera = 0;
   // z = 1;
   // zet = 1;

   // for (int i = 0; i < nVar; i++)
   // {
   //    epsvecn[i] = epsi;
   //    x[i] = 0.5 * (alfa[i] + beta[i]);
   //    xsi[i] = 1.0/(x[i] - alfa[i]);
   //    xsi[i] = std::max(xsi[i], 1.0);
   //    eta[i] = 1.0/(beta[i] - x[i]);
   //    eta[i] = std::max(eta[i], 1.0);
   // }
   // for (int i = 0; i < nCon; i++)
   // {
   //    epsvecm[i] = epsi;
   //    y[i] = 1.0;
   //    lam[i] = 1.0;
   //    mu[i] = std::max(1.0, 0.5 * c[i]);
   //    s[i] = 1.0;
   // }

   // //for (int Ik = 0; true; Ik++)
   // while (epsi > epsimin)
   // {

   //    rez = a0 - zet;
   //    for (int i = 0; i < nVar; i++)
   //    {
   //       epsvecn[i] = epsi;
   //       ux1[i] = upp[i] - x[i];
   //       if (ux1[i] == 0)
   //       {
   //          ux1[i] = machineEpsilon;
   //       }
   //       else
   //       {
   //          while (std::fabs(ux1[i]) <= machineEpsilon)
   //          {
   //             ux1[i] *= 10;
   //          }
   //       }
   //       xl1[i] = x[i] - low[i];
   //       if (xl1[i] == 0)
   //       {
   //          xl1[i] = machineEpsilon;
   //       }
   //       else
   //       {
   //          while (std::fabs(xl1[i]) <= machineEpsilon)
   //          {
   //             xl1[i] *= 10;
   //          }
   //       }
   //       // plam = P' * lam, qlam = Q' * lam
   //       plam[i] = p0[i];
   //       qlam[i] = q0[i];
   //       for (int j = 0; j < nCon; j++)
   //       {
   //          plam[i] += P[j * nVar + i] * lam[j];
   //          qlam[i] += Q[j * nVar + i] * lam[j];
   //       }
   //       dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
   //       rex[i] = dpsidx[i] - xsi[i] + eta[i];
   //       rez -= a[i] * lam[i];
   //       rexsi[i] = xsi[i] * (x[i] - alfa[i]) - epsvecn[i];
   //       if (std::fabs(x[i]-alfa[i]) == 0)
   //       {
   //          rexsi[i] = xsi[i] * machineEpsilon - epsvecn[i];
   //       }
   //       reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
   //       if (std::fabs(beta[i] - x[i]) == 0)
   //       {
   //          reeta[i] = eta[i] * machineEpsilon - epsvecn[i];;
   //       }
   //       residu1[i] = rex[i];
   //       residu2[nCon + i] = rexsi[i];
   //       residu2[nCon + nVar + i] = reeta[i];
   //    }
   //    for (int i = 0; i < nCon; i++)
   //    {
   //       epsvecm[i] = epsi;
   //       // gvec = P/ux + Q/xl
   //       for (int j = 0; j < nVar; j++)
   //       {
   //          gvec[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
   //       }
   //       rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
   //       relam[i] = gvec[i] - a[i] * z - y[i] + s[i] - b[i];
   //       remu[i] = mu[i] * y[i] - epsvecm[i];
   //       res[i] = lam[i] * s[i] - epsvecm[i];
   //       residu2[i] = relam[i];
   //       residu1[nVar + i] = rey[i];
   //       residu2[nCon + 2 * nVar + i] = remu[i];
   //       residu2[2 * nVar + 2 * nCon + 1 + i] = res[i];
   //    }
   //    rezet = zet * z - epsi;
   //    residu1[nVar + nCon] = rez;
   //    residu2[2 * nVar + 2 * nCon] = rezet;

   //    // Concatenate the residuals
   //    for (int i = 0; i < (nVar + nCon + 1); i++)
   //    {
   //       residu[i] = residu1[i];
   //    }
   //    for (int i = 0; i < (2 * nVar + 3 * nCon + 1); i++)
   //    {
   //       residu[nVar + nCon + 1 + i] = residu2[i];
   //    }

   //    //Get vector product and maximum absolute value
   //    residunorm = 0.0;
   //    residumax = 0.0;
   //    for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
   //    {
   //       residunorm += residu[i] * residu[i];
   //       residumax = std::max(residumax, std::abs(residu[i]));
   //    }
   //    // Norm of the residual
   //    residunorm = std::sqrt(residunorm);

   //    ittt = 0;

   //    while (residumax > 0.9 * epsi && ittt < 200)
   //    {
   //       ittt++;
   //       for (int i = 0; i < nVar; i++)
   //       {
   //          ux1[i] = upp[i] - x[i];
   //          if (ux1[i] == 0)
   //          {
   //             ux1[i] = machineEpsilon;
   //          }
   //          else
   //          {
   //             while (std::fabs(ux1[i]) <= machineEpsilon)
   //             {
   //                ux1[i] *= 10;
   //             }
   //          }
   //          xl1[i] = x[i] - low[i];
   //          if (xl1[i] == 0)
   //          {
   //             xl1[i] = machineEpsilon;
   //          }
   //          else
   //          {
   //             while (std::fabs(xl1[i]) <= machineEpsilon)
   //             {
   //                xl1[i] *= 10;
   //             }
   //          }
   //          // plam = P' * lam, qlam = Q' * lam
   //          plam[i] = p0[i];
   //          qlam[i] = q0[i];
   //          for (int j = 0; j < nCon; j++)
   //          {
   //             plam[i] += P[j * nVar + i] * lam[j];
   //             qlam[i] += Q[j * nVar + i] * lam[j];
   //          }
   //          dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
   //          // NaN-Avoidance
   //          if (std::fabs(x[i] - alfa[i]) < machineEpsilon)
   //          {
   //             if (std::fabs(beta[i] - x[i]) < machineEpsilon)
   //             {
   //                delx[i] = dpsidx[i];
   //                diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / machineEpsilon + eta[i] / machineEpsilon;
   //             }
   //             else
   //             {
   //                delx[i] = dpsidx[i] - epsvecn[i] / machineEpsilon + epsvecn[i] / (beta[i] - x[i]);
   //                diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / (x[i] - alfa[i]) + eta[i] / (beta[i] - x[i]);
   //             }
   //          }
   //          else if (std::fabs(beta[i] - x[i]) < machineEpsilon)
   //          {
   //             delx[i] = dpsidx[i] - epsvecn[i] / (x[i] - alfa[i]) + epsvecn[i] / machineEpsilon;
   //             diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / (x[i] - alfa[i]) + eta[i] / machineEpsilon;
   //          }
   //          else
   //          {
   //             delx[i] = dpsidx[i] - epsvecn[i] / (x[i] - alfa[i]) + epsvecn[i] / (beta[i] - x[i]);
   //             diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + xsi[i] / (x[i] - alfa[i]) + eta[i] / (beta[i] - x[i]);
   //          }
   //       }

   //       delz = a0 - epsi / z;
   //       for (int i = 0; i < nCon; i++)
   //       {
   //          gvec[i] = 0.0;
   //          // gvec = P/ux + Q/xl
   //          for (int j = 0; j < nVar; j++)
   //          {
   //             gvec[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
   //             Puxinv[i * nVar + j] = P[i * nVar + j] / (ux1[j] * ux1[j]);
   //             Qxlinv[i * nVar + j] = Q[i * nVar + j] / (xl1[j] * xl1[j]);

   //             GG[i * nVar + j] = Puxinv[i * nVar + j] - Qxlinv[i * nVar + j];
   //          }

   //          dely[i] = c[i] + d[i] * y[i] - lam[i] - epsvecm[i] / y[i];
   //          delz -= a[i] * lam[i];
   //          dellam[i] = gvec[i] - a[i] * z - y[i] - b[i] + epsvecm[i] / lam[i];
   //          diagy[i] = d[i] + mu[i] / y[i];
   //          //diagyinv[i] = 1.0 / diagy[i];
   //          diaglam[i] = s[i] / lam[i];
   //          diaglamyi[i] = diaglam[i] + 1.0 / diagy[i];
   //       }

   //       if (nCon < nVar)
   //       {
   //          // blam = dellam + dely./diagy - GG*(delx./diagx);
   //          // bb = [blam; delz];
   //          for (int j = 0; j < nCon; j++)
   //          {
   //             sum = 0.0;
   //             for (int i = 0; i < nVar; i++)
   //             {
   //                sum += GG[j * nVar + i] * (delx[i] / diagx[i]);
   //             }
   //             blam[j] = dellam[j] + dely[j] / diagy[j] - sum;
   //             bb1[j] = blam[j];
   //          }
   //          bb1[nCon] = delz;            

   //          // Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             // Axx = GG*spdiags(diagxinv,0,n,n);
   //             for (int k = 0; k < nVar; k++)
   //             {
   //                Axx[i * nVar + k] = GG[k * nCon + i] / diagx[k];
   //             }
   //          }
   //          // Alam = spdiags(diaglamyi,0,m,m) + Axx*GG';
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             for (int j = 0; j < nCon; j++)
   //             {
   //                Alam[i * nCon + j] = 0.0;
   //                for (int k = 0; k < nVar; k++)
   //                {
   //                   Alam[i * nCon + j] += Axx[i * nVar + k] * GG[j * nVar + k];
   //                }
   //                if (i == j)
   //                {
   //                   Alam[i * nCon + j] += diaglamyi[i];
   //                }
   //             }
   //          }
   //          // AA1 = [Alam     a
   //          //       a'    -zet/z];
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             for (int j = 0; j < nCon; j++)
   //             {
   //                AA1[i * nCon + j] = Alam[i * nCon + j];
   //             }
   //             AA1[nCon * nCon + i] = a[i];
   //             AA1[nCon * nCon + nCon + i] = a[i];
   //          }
   //          AA1[nCon * nCon + nCon + nCon] = -zet / z;                              

   //          // ----------------------------------------------------------------------------
   //          //solut = AA1\bb1 --> solve linear system of equations using LAPACK
   //          int info;
   //          int nLAP = nCon + 1;
   //          int nrhs = 1;
   //          int lda = nLAP;
   //          int ldb = nLAP;
   //          int* ipiv = new int[nLAP];
   //          dgesv_(&nLAP, &nrhs, AA1, &lda, ipiv, bb1, &ldb, &info);
   //          delete[] ipiv;
   //          for (int i = 0; i < (nCon + 1); i++)
   //          {
   //             solut[i] = bb1[i];
   //          }
   //          // ----------------------------------------------------------------------------
   //          //dlam = solut(1:nCon);
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             dlam[i] = solut[i];
   //          }
   //          dz = solut[nCon];
   //          //dx = -(GG'*dlam)./diagx - delx./diagx;
   //          for (int i = 0; i < nVar; i++)
   //          {
   //             sum = 0.0;
   //             for (int j = 0; j < nCon; j++)
   //             {
   //                sum += GG[j * nVar + i] * dlam[j];
   //             }
   //             dx[i] = -sum / diagx[i] - delx[i] / diagx[i];           
   //          }
   //       }
   //       else
   //       {
   //          sum = 0.0;
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             //diaglamyiinv[i] = 1.0 / diaglamyi[i];
   //             dellamyi[i] = dellam[i] + dely[i] / diagy[i];
   //             // azz = zet/z + a'*(a./diaglamyi)
   //             sum += a[i] * (a[i] / diaglamyi[i]);
   //          }
   //          azz = zet / z + sum;
   //          // Axx = spdiags(diagx,0,nVar,nVar) + GG'*spdiags(diaglamyiinv,0,nCon,nCon)*GG;
   //          // AA = [Axx      axz
   //          //       axz'     azz];
   //          for (int i = 0; i < nVar; i++)
   //          {
   //             // Axx =  GG'*spdiags(diaglamyiinv,0,nCon,nCon);
   //             for (int k = 0; k < nCon; k++)
   //             {
   //                Axx[i * nCon + k] = GG[k * nVar + i] / diaglamyi[k];
   //             }
   //             sum = 0.0;
   //             // axz = -GG'*(a./diaglamyi)
   //             for (int j = 0; j < nCon; j++)
   //             {
   //                sum -= GG[j * nVar + i] * (a[j] / diaglamyi[j]);
   //             }
   //             axz[i] = sum;
   //          }
   //          //Assemble matrix AA
   //          for (int i = 0; i < (nVar + 1); i++)
   //          {
   //             for (int j = 0; j < (nVar + 1); j++)
   //             {
   //                // AA = [Axx  .
   //                //       .    .]
   //                AA[i * (nVar + 1) + j] = 0.0;
   //                if (i < nVar && j < nVar)
   //                {
   //                   // Axx = Axx*GG
   //                   for (int k = 0; k < nCon; k++)
   //                   {
   //                      AA[i * (nVar + 1) + j] += Axx[i * nCon + k] * GG[k * nVar + j];
   //                   }
   //                   // Axx = Axx + spdiags(diagx,0,nVar,nVar)
   //                   if (i == j)
   //                   {
   //                      AA[i * (nVar + 1) + j] += diagx[j];
   //                   }
   //                }
   //                // AA = [Axx  axz
   //                //       axz' azz]
   //                else if (i < nVar && j == nVar)
   //                {
   //                   AA[i * (nVar + 1) + j] = axz[i];
   //                }
   //                else if (i == nVar && j < nVar)
   //                {
   //                   AA[i * (nVar + 1) + j] = axz[j];
   //                }
   //                else
   //                {
   //                   AA[i * (nVar + 1) + j] = azz;
   //                }
   //             }
   //          }
   //          // bb = [-bx'; -bz]
   //          // bx = delx - GG'*(dellamyi./diaglamyi)
   //          // bz = delz - a'*(dellamyi./diaglamyi)
   //          for (int i = 0; i < nVar; i++)
   //          {
   //             sum = 0.0;
   //             for (int j = 0; j < nCon; j++)
   //             {
   //                sum += GG[j * nVar + i] * (dellamyi[j] / diaglamyi[j]);
   //             }
   //             bb[i] = -(delx[i] + sum);
   //          }
   //          sum = 0.0;
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             sum += a[i] * (dellamyi[i] / diaglamyi[i]);
   //          }
   //          bb[nVar] = -(delz - sum);
   //          // ----------------------------------------------------------------------------
   //          //solut = AA\bb --> solve linear system of equations using LAPACK
   //          int info;
   //          int nLAP = nVar + 1;
   //          int nrhs = 1;
   //          int lda = nLAP;
   //          int ldb = nLAP;
   //          int* ipiv = new int[nLAP];
   //          dgesv_(&nLAP, &nrhs, AA, &lda, ipiv, bb, &ldb, &info);
   //          delete[] ipiv;
   //          for (int i = 0; i < (nVar + 1); i++)
   //          {
   //             solut[i] = bb[i];
   //          }
   //          // ----------------------------------------------------------------------------
   //          //dx = solut(1:nVar);
   //          for (int i = 0; i < nVar; i++)
   //          {
   //             dx[i] = solut[i];
   //          }
   //          dz = solut[nVar];
   //          //dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             sum = 0.0;
   //             for (int j = 0; j < nVar; j++)
   //             {
   //                sum += GG[i * nVar + j] * dx[j];
   //             }
   //             dlam[i] = sum / diaglamyi[i] - dz * (a[i] / diaglamyi[i]) + dellamyi[i] / diaglamyi[i];
   //          }
   //       }

   //       for (int i = 0; i < nCon; i++)
   //       {
   //          dy[i] = -dely[i] / diagy[i] + dlam[i] / diagy[i];
   //          dmu[i] = -mu[i] + epsvecm[i] / y[i] - (mu[i] * dy[i]) / y[i];
   //          ds[i] = -s[i] + epsvecm[i] / lam[i] - (s[i] * dlam[i]) / lam[i];
   //          // xx = [y z lam xsi eta mu zet s]
   //          // dxx = [dy dz dlam dxsi deta dmu dzet ds]
   //          xx[i] = y[i];
   //          xx[nCon + 1 + i] = lam[i];
   //          xx[2 * nCon + 1 + 2 * nVar + i] = mu[i];
   //          xx[3 * nCon + 2 * nVar + 2 + i] = s[i];
   //          dxx[i] = dy[i];
   //          dxx[nCon + 1 + i] = dlam[i];
   //          dxx[2 * nCon + 1 + 2 * nVar + i] = dmu[i];
   //          dxx[3 * nCon + 2 * nVar + 2 + i] = ds[i];
   //       }
   //       xx[nCon] = z;
   //       xx[3 * nCon + 2 * nVar + 1] = zet;
   //       dxx[nCon] = dz;
   //       dxx[3 * nCon + 2 * nVar + 1] = dzet;
   //       for (int i = 0; i < nVar; i++)
   //       {
   //          // NaN-Avoidance
   //          if (std::fabs(x[i] - alfa[i]) < machineEpsilon)
   //          {
   //             if (std::fabs(beta[i] - x[i]) < machineEpsilon)
   //             {
   //                dxsi[i] = -xsi[i] + epsvecn[i] / machineEpsilon - (xsi[i] * dx[i]) / machineEpsilon;
   //                deta[i] = -eta[i] + epsvecn[i] / machineEpsilon + (eta[i] * dx[i]) / machineEpsilon;
   //             }
   //             else
   //             {
   //                dxsi[i] = -xsi[i] + epsvecn[i] / machineEpsilon - (xsi[i] * dx[i]) / machineEpsilon;
   //                deta[i] = -eta[i] + epsvecn[i] / (beta[i] - x[i]) + (eta[i] * dx[i]) / (beta[i] - x[i]);
   //             }
   //          }
   //          else if (std::fabs(beta[i] - x[i]) < machineEpsilon)
   //          {
   //             dxsi[i] = -xsi[i] + epsvecn[i] / (x[i] - alfa[i]) - (xsi[i] * dx[i]) / (x[i] - alfa[i]);
   //             deta[i] = -eta[i] + epsvecn[i] / machineEpsilon + (eta[i] * dx[i]) / machineEpsilon;
   //          }
   //          else
   //          {
   //             dxsi[i] = -xsi[i] + epsvecn[i] / (x[i] - alfa[i]) - (xsi[i] * dx[i]) / (x[i] - alfa[i]);
   //             deta[i] = -eta[i] + epsvecn[i] / (beta[i] - x[i]) + (eta[i] * dx[i]) / (beta[i] - x[i]);
   //          }
   //          xx[nCon + 1 + nCon + i] = xsi[i];
   //          xx[nCon + 1 + nCon + nVar + i] = eta[i];
   //          dxx[nCon + 1 + nCon + i] = dxsi[i];
   //          dxx[nCon + 1 + nCon + nVar + i] = deta[i];
   //       }
   //       dzet = -zet + epsi / z - zet * dz / z;

   //       stmxx = 0.0;
   //       for (int i = 0; i < (4 * nCon + 2 * nVar + 2); i++)
   //       {
   //          stepxx[i] = -1.01*dxx[i] / xx[i];
   //          if (stepxx[i] > stmxx)
   //          {
   //             stmxx = stepxx[i];
   //          }
   //       }
   //       stmalfa = 0.0;
   //       stmbeta = 0.0;
   //       for (int i = 0; i < nVar; i++)
   //       {
            
   //          //NaN-Avoidance
   //          if (std::fabs(x[i] - alfa[i]) < machineEpsilon)
   //          {
   //             stepalfa[i] = -1.01*dx[i] / machineEpsilon;
   //          }
   //          else 
   //          {
   //             stepalfa[i] = -1.01*dx[i] / (x[i] - alfa[i]);
   //          }
   //          if (std::fabs(beta[i] - x[i]) < machineEpsilon)
   //          {
   //             stepbeta[i] = 1.01*dx[i] / machineEpsilon;
   //          }
   //          else
   //          {
   //             stepbeta[i] = 1.01*dx[i] / (beta[i] - x[i]);
   //          }
   //          // --------------
   //          if (stepalfa[i] > stmalfa)
   //          {
   //             stmalfa = stepalfa[i];
   //          }
   //          if (stepbeta[i] > stmbeta)
   //          {
   //             stmbeta = stepbeta[i];
   //          }
   //       }
   //       stmalbe = std::max(stmalfa, stmbeta);
   //       stmalbexx = std::max(stmalbe, stmxx);
   //       stminv = std::max(stmalbexx, 1.0);
   //       steg = 1.0 / stminv;

   //       for (int i = 0; i < nVar; i++)
   //       {
   //          xold[i] = x[i];
   //          xsiold[i] = xsi[i];
   //          etaold[i] = eta[i];
   //       }
   //       for (int i = 0; i < nCon; i++)
   //       {
   //          yold[i] = y[i];
   //          lamold[i] = lam[i];
   //          muold[i] = mu[i];
   //          sold[i] = s[i];
   //       }
   //       zold = z;
   //       zetold = zet;

   //       itto = 0;
   //       resinew = 2.0 * residunorm;
   //       while (resinew > residunorm && itto < 50)
   //       {
   //          itto++;

   //          for (int i = 0; i < nCon; ++i)
   //          {
   //             y[i] = yold[i] + steg * dy[i];
   //             if (std::fabs(y[i]) < machineEpsilon)
   //             {
   //                y[i] = machineEpsilon;
   //             }
   //             else
   //             {
   //                while (std::fabs(y[i]) < machineEpsilon)
   //                {
   //                   y[i] *= 10;
   //                }
   //             }
   //             lam[i] = lamold[i] + steg * dlam[i];
   //             if (std::fabs(lam[i]) < machineEpsilon)
   //             {
   //                lam[i] = machineEpsilon;
   //             }
   //             else
   //             {
   //                while (std::fabs(lam[i]) < machineEpsilon)
   //                {
   //                   lam[i] *= 10;
   //                }
   //             }
   //             mu[i] = muold[i] + steg * dmu[i];
   //             s[i] = sold[i] + steg * ds[i];
   //          }

   //          for (int i = 0; i < nVar; ++i)
   //          {
   //             x[i] = xold[i] + steg * dx[i];
   //             xsi[i] = xsiold[i] + steg * dxsi[i];
   //             eta[i] = etaold[i] + steg * deta[i];
   //             if (std::isnan(eta[i]))
   //             {
   //                eta[i] = machineEpsilon;
   //             }
   //             ux1[i] = upp[i] - x[i];
   //             if (ux1[i] == 0)
   //             {
   //                ux1[i] = machineEpsilon;
   //             }
   //             else
   //             {
   //                while (std::fabs(ux1[i]) <= machineEpsilon)
   //                {
   //                   ux1[i] *= 10;
   //                }
   //             }
   //             xl1[i] = x[i] - low[i];
   //             if (xl1[i] == 0)
   //             {
   //                xl1[i] = machineEpsilon;
   //             }
   //             else
   //             {
   //                while (std::fabs(xl1[i]) <= machineEpsilon)
   //                {
   //                   xl1[i] *= 10;
   //                }
   //             }
   //             // plam & qlam
   //             plam[i] = p0[i];
   //             qlam[i] = q0[i];
   //             for (int j = 0; j < nCon; j++)
   //             {
   //                plam[i] += P[j * nVar + i] * lam[j];
   //                qlam[i] += Q[j * nVar + i] * lam[j];
   //             }
   //             dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
   //             rex[i] = dpsidx[i] - xsi[i] + eta[i];
   //             rez -= a[i] * lam[i];
   //             rexsi[i] = xsi[i] * (x[i] - alfa[i]) - epsvecn[i];              
   //             if (std::fabs(x[i] - alfa[i]) < machineEpsilon || std::isnan(rexsi[i]))
   //             {
   //                rexsi[i] = machineEpsilon;
   //             }
   //             reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
   //             if (std::fabs(beta[i] - x[i]) < machineEpsilon || std::isnan(reeta[i]))
   //             {
   //                reeta[i] = machineEpsilon;
   //             }
   //             residu1[i] = rex[i];
   //             residu2[nCon + i] = rexsi[i];
   //             residu2[nCon + nVar + i] = reeta[i];
   //          }
   //          z = zold + steg * dz;
   //          if (z == 0)
   //          {
   //             z = machineEpsilon;
   //          }
   //          else
   //          {
   //             while (std::fabs(z) <= machineEpsilon)
   //             {
   //                z *= 10;
   //             }
   //          }
   //          zet = zetold + steg * dzet;

   //          // gvec = P/ux + Q/xl
   //          for (int i = 0; i < nCon; i++)
   //          {
   //             gvec[i] = 0.0;
   //             for (int j = 0; j < nVar; j++)
   //             {
   //                gvec[i] += P[i * nVar + j] / ux1[j] + Q[i * nVar + j] / xl1[j];
   //             }
   //             rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
   //             relam[i] = gvec[i] - a[i] * z - y[i] + s[i] - b[i];
   //             remu[i] = mu[i] * y[i] - epsvecm[i];
   //             res[i] = lam[i] * s[i] - epsvecm[i];

   //             residu2[i] = relam[i];
   //             residu1[nVar + i] = rey[i];
   //             residu2[nCon + 2 * nVar + i] = remu[i];
   //             residu2[2 * nVar + 2 * nCon + 1 + i] = res[i];
   //          }
   //          rez = a0 - zet;
   //          rezet = zet * z - epsi;
   //          residu1[nVar + nCon] = rez;
   //          residu2[2 * nVar + 2 * nCon] = rezet;

   //          // Concatenate the residuals
   //          for (int i = 0; i < (nVar + nCon + 1); i++)
   //          {
   //             residu[i] = residu1[i];
   //          }
   //          for (int i = 0; i < (2 * nVar + 3 * nCon + 1); i++)
   //          {
   //             residu[nVar + nCon + 1 + i] = residu2[i];
   //          }
   //          // New residual
   //          sum = 0.0;
   //          for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
   //          {
   //             sum += residu[i] * residu[i];
   //          }
   //          // Norm of the residual
   //          resinew = std::sqrt(sum);
   //          steg = steg / 2.0;
   //       }
   //       residunorm = resinew;
   //       residumax = 0.0;
   //       for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
   //       {
   //          residumax = std::max(residumax, std::abs(residu[i]));
   //       }
   //       steg = steg * 2.0;
   //    }


   //    if (ittt > 198)
   //    {
   //       printf("Warning: Maximum number of iterations reached in subsolv.\n");
   //    }
   //    epsi = 0.1 * epsi;
   // }
   // //results.close();
}

void MMA::kktcheck(double* y, double* const dfdx, double* const gx, double* const dgdx, double* x)
{
   //std::ofstream kkt;
   //kkt.open("KKT.dat", std::ios::app);

   for (int i = 0; i < nVar; i++)
   {
      sum1[i] = 0.0;
      for (int j = 0; j < nCon; j++)
      {
         sum1[i] += dgdx[j * nVar + i] * lam[j];
      }
   }
   for (int i = 0; i < nVar; i++)
   {
      rex[i] = dfdx[i] + sum1[i] - xsi[i] + eta[i];
      rexsi[i] = xsi[i] * (x[i] - xmin[i]);      
      reeta[i] = eta[i] * (xmax[i] - x[i]);
   }
   for (int i = 0; i < nCon; i++)
   {
      rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
      relam[i] = gx[i] - a[i] * z - y[i] + s[i];
      remu[i] = mu[i] * y[i];
      res[i] = lam[i] * s[i];
   }
   sum = 0.0;
   for (int i = 0; i < nCon; i++)
   {
      sum += a[i] * lam[i];
   }
   rez = a0 - zet - sum;
   rezet = zet * z;
   //---------------------------------------------------------------------

   for (int i = 0; i < nVar; i++)
   {
      residu1[i] = rex[i];
      residu2[nCon + i] = rexsi[i];
      residu2[nCon + nVar + i] = reeta[i];
   }

   for (int i = 0; i < nCon; i++)
   {
      residu2[i] = relam[i];
      residu1[nVar + i] = rey[i];
      residu2[nCon + 2 * nVar + i] = remu[i];
      residu2[2 * nVar + 2 * nCon + 1 + i] = res[i];
   }
   residu1[nVar + nCon] = rez;
   residu2[2 * nVar + 2 * nCon] = rezet;

   // Concatenate the residuals
   for (int i = 0; i < (nVar + nCon + 1); i++)
   {
      residu[i] = residu1[i];
   }
   for (int i = 0; i < (2 * nVar + 3 * nCon + 1); i++)
   {
      residu[nVar + nCon + 1 + i] = residu2[i];
   }
   //Get vector product and maximum absolute value
   residunorm = 0.0;
   residumax = 0.0;
   for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
   {
      residunorm += residu[i] * residu[i];
      residumax = std::max(residumax, std::abs(residu[i]));
   }
   // Norm of the residual
   residunorm = std::sqrt(residunorm);
   kktnorm = residunorm;
   
   //kkt.close();
}

void MMA::setRestart(double* xval, int iter, std::string name, std::string path)
{
   if (iter < 3)
   {
      mfem::mfem_error("MMA::setRestart - iter must be greater than 2");
   }
   
   std::string restartFile = path + name;
   std::ifstream input(restartFile);
   input >> iter;
   printf("iter = %d\n", iter);
   for (int i = 0; i < nVar; i++)
   {
      input >> xval[i];
      printf("xval[%d] = %f\n", i, xval[i]);
   }
   for (int i = 0; i < nVar; i++)
   {
      input >> xo1[i];
      printf("xo1[%d] = %f\n", i, xo1[i]);
   }
   for (int i = 0; i < nVar; i++)
   {
      input >> xo2[i];
      printf("xo2[%d] = %f\n", i, xo2[i]);
   }
   for (int i = 0; i < nVar; i++)
   {
      input >> upp[i];
      printf("upp[%d] = %f\n", i, upp[i]);
   }
   for (int i = 0; i < nVar; i++)
   {
      input >> low[i];
      printf("low[%d] = %f\n", i, low[i]);
   }
   input.close();
}

void MMA::outputRestart(double* xval, int iter, std::string name, std::string path)
{
   std::string restartFile = path + name;
   std::ofstream mma;
   mma.open(restartFile);
   //print results
   mma << iter << "\n";
   for (int i = 0; i < nVar; i++)
   {
      mma << xval[i] << "\n";
   }
   for (int i = 0; i < nVar; i++)
   {
      mma << xo1[i] << "\n";
   }
   for (int i = 0; i < nVar; i++)
   {
      mma << xo2[i] << "\n";
   }
   for (int i = 0; i < nVar; i++)
   {
      mma << upp[i] << "\n";
   }
   for (int i = 0; i < nVar; i++)
   {
      mma << low[i] << "\n";
   }
   mma.close();
}

double* MMA::getLow()
{
   return low;
}
double* MMA::getUpp()
{
   return upp;
}

double MMA::getKKT()
{
   return kktnorm;
}

} // end mma namespace
