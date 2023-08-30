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

//Constructors
// ----------------------------------------------------------------------
MMA::MMA(int nVar, int nCon, double *xval, double *xxmin, double *xxmax, enum SubProblemType type )
{
   this->setGlobals(nVar, nCon);
   this->setMMA(nVar, nCon);
   // this->setSubProb(nVar, nCon);
   this->initializeGlobals(nVar, nCon, xval, xxmin, xxmax);

   this->subProblemFactory(nVar, nCon, type);
}

MMA::SubProblemBase::SubProblemBase(MMA* mma) : mma_ptr(mma){}

MMA::SubProblemClassic::SubProblemClassic(MMA* mma, int nVar, int nCon) : SubProblemBase(mma)
{
   this->setSubProb(nVar, nCon);
}

MMA::SubProblemClassicMPI::SubProblemClassicMPI(MMA* mma, int nVar, int nCon) : SubProblemBase(mma)
{
   this->setSubProb(nVar, nCon);
}



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

   // passing in: nVar, nCon, epsimin, low, upp, alfa, beta, p0, q0, P,Q,a0,a,b,c,d
   mSubProblem->Perform();

   // Update design variables
   for (int i = 0; i < nVar; i++)
   {
      xo2[i] = xo1[i];
      xo1[i] = xval[i];
      xval[i] = x[i];
   }   
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
   mma_ptr->z = 1;
   //z = 1;
   mma_ptr->zet = 1;

   for (int i = 0; i < mma_ptr->nVar; i++)
   {
      epsvecn[i] = epsi;
      mma_ptr->x[i] = 0.5 * (mma_ptr->alfa[i] + mma_ptr->beta[i]);
      mma_ptr->xsi[i] = 1.0/(mma_ptr->x[i] - mma_ptr->alfa[i]);
      mma_ptr->xsi[i] = std::max(mma_ptr->xsi[i], 1.0);
      mma_ptr->eta[i] = 1.0/(mma_ptr->beta[i] - mma_ptr->x[i]);
      mma_ptr->eta[i] = std::max(mma_ptr->eta[i], 1.0);
   }
   for (int i = 0; i < mma_ptr->nCon; i++)
   {
      epsvecm[i] = epsi;
      mma_ptr->y[i] = 1.0;
      mma_ptr->lam[i] = 1.0;
      mma_ptr->mu[i] = std::max(1.0, 0.5 * mma_ptr->c[i]);
      mma_ptr->s[i] = 1.0;
   }

   //for (int Ik = 0; true; Ik++)
   while (epsi > mma_ptr->epsimin)
   {

      rez = mma_ptr->a0 - mma_ptr->zet;
      for (int i = 0; i < mma_ptr->nVar; i++)
      {
         epsvecn[i] = epsi;
         ux1[i] = mma_ptr->upp[i] - mma_ptr->x[i];
         if (ux1[i] == 0)
         {
            ux1[i] = mma_ptr->machineEpsilon;
         }
         else
         {
            while (std::fabs(ux1[i]) <= mma_ptr->machineEpsilon)
            {
               ux1[i] *= 10;
            }
         }
         xl1[i] = mma_ptr->x[i] - mma_ptr->low[i];
         if (xl1[i] == 0)
         {
            xl1[i] = mma_ptr->machineEpsilon;
         }
         else
         {
            while (std::fabs(xl1[i]) <= mma_ptr->machineEpsilon)
            {
               xl1[i] *= 10;
            }
         }
         // plam = P' * lam, qlam = Q' * lam
         plam[i] = mma_ptr->p0[i];
         qlam[i] = mma_ptr->q0[i];
         for (int j = 0; j < mma_ptr->nCon; j++)
         {
            plam[i] += mma_ptr->P[j * mma_ptr->nVar + i] * mma_ptr->lam[j];
            qlam[i] += mma_ptr->Q[j * mma_ptr->nVar + i] * mma_ptr->lam[j];
         }
         dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
         rex[i] = dpsidx[i] - mma_ptr->xsi[i] + mma_ptr->eta[i];
         rez -= mma_ptr->a[i] * mma_ptr->lam[i];
         rexsi[i] = mma_ptr->xsi[i] * (mma_ptr->x[i] - mma_ptr->alfa[i]) - epsvecn[i];
         if (std::fabs(mma_ptr->x[i]-mma_ptr->alfa[i]) == 0)
         {
            rexsi[i] = mma_ptr->xsi[i] * mma_ptr->machineEpsilon - epsvecn[i];
         }
         reeta[i] = mma_ptr->eta[i] * (mma_ptr->beta[i] - mma_ptr->x[i]) - epsvecn[i];
         if (std::fabs(mma_ptr->beta[i] - mma_ptr->x[i]) == 0)
         {
            reeta[i] = mma_ptr->eta[i] * mma_ptr->machineEpsilon - epsvecn[i];;
         }
         residu1[i] = rex[i];
         residu2[mma_ptr->nCon + i] = rexsi[i];
         residu2[mma_ptr->nCon + mma_ptr->nVar + i] = reeta[i];
      }
      for (int i = 0; i < mma_ptr->nCon; i++)
      {
         epsvecm[i] = epsi;
         // gvec = P/ux + Q/xl
         for (int j = 0; j < mma_ptr->nVar; j++)
         {
            gvec[i] += mma_ptr->P[i * mma_ptr->nVar + j] / ux1[j] + mma_ptr->Q[i * mma_ptr->nVar + j] / xl1[j];
         }
         rey[i] = mma_ptr->c[i] + mma_ptr->d[i] * mma_ptr->y[i] - mma_ptr->mu[i] - mma_ptr->lam[i];
         relam[i] = gvec[i] - mma_ptr->a[i] * mma_ptr->z - mma_ptr->y[i] + mma_ptr->s[i] - mma_ptr->b[i];
         remu[i] = mma_ptr->mu[i] * mma_ptr->y[i] - epsvecm[i];
         res[i] = mma_ptr->lam[i] * mma_ptr->s[i] - epsvecm[i];
         residu2[i] = relam[i];
         residu1[mma_ptr->nVar + i] = rey[i];
         residu2[mma_ptr->nCon + 2 * mma_ptr->nVar + i] = remu[i];
         residu2[2 * mma_ptr->nVar + 2 * mma_ptr->nCon + 1 + i] = res[i];
      }
      rezet = mma_ptr->zet * mma_ptr->z - epsi;
      residu1[mma_ptr->nVar + mma_ptr->nCon] = rez;
      residu2[2 * mma_ptr->nVar + 2 * mma_ptr->nCon] = rezet;

      // Concatenate the residuals
      for (int i = 0; i < (mma_ptr->nVar + mma_ptr->nCon + 1); i++)
      {
         residu[i] = residu1[i];
      }
      for (int i = 0; i < (2 * mma_ptr->nVar + 3 * mma_ptr->nCon + 1); i++)
      {
         residu[mma_ptr->nVar + mma_ptr->nCon + 1 + i] = residu2[i];
      }

      //Get vector product and maximum absolute value
      residunorm = 0.0;
      residumax = 0.0;
      for (int i = 0; i < (3 * mma_ptr->nVar + 4 * mma_ptr->nCon + 2); i++)
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
         for (int i = 0; i < mma_ptr->nVar; i++)
         {
            ux1[i] = mma_ptr->upp[i] - mma_ptr->x[i];
            if (ux1[i] == 0)
            {
               ux1[i] = mma_ptr->machineEpsilon;
            }
            else
            {
               while (std::fabs(ux1[i]) <= mma_ptr->machineEpsilon)
               {
                  ux1[i] *= 10;
               }
            }
            xl1[i] = mma_ptr->x[i] - mma_ptr->low[i];
            if (xl1[i] == 0)
            {
               xl1[i] = mma_ptr->machineEpsilon;
            }
            else
            {
               while (std::fabs(xl1[i]) <= mma_ptr->machineEpsilon)
               {
                  xl1[i] *= 10;
               }
            }
            // plam = P' * lam, qlam = Q' * lam
            plam[i] = mma_ptr->p0[i];
            qlam[i] = mma_ptr->q0[i];
            for (int j = 0; j < mma_ptr->nCon; j++)
            {
               plam[i] += mma_ptr->P[j * mma_ptr->nVar + i] * mma_ptr->lam[j];
               qlam[i] += mma_ptr->Q[j * mma_ptr->nVar + i] * mma_ptr->lam[j];
            }
            dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
            // NaN-Avoidance
            if (std::fabs(mma_ptr->x[i] - mma_ptr->alfa[i]) < mma_ptr->machineEpsilon)
            {
               if (std::fabs(mma_ptr->beta[i] - mma_ptr->x[i]) < mma_ptr->machineEpsilon)
               {
                  delx[i] = dpsidx[i];
                  diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + mma_ptr->xsi[i] / mma_ptr->machineEpsilon + mma_ptr->eta[i] / mma_ptr->machineEpsilon;
               }
               else
               {
                  delx[i] = dpsidx[i] - epsvecn[i] / mma_ptr->machineEpsilon + epsvecn[i] / (mma_ptr->beta[i] - mma_ptr->x[i]);
                  diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + mma_ptr->xsi[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]) + mma_ptr->eta[i] / (mma_ptr->beta[i] - mma_ptr->x[i]);
               }
            }
            else if (std::fabs(mma_ptr->beta[i] - mma_ptr->x[i]) < mma_ptr->machineEpsilon)
            {
               delx[i] = dpsidx[i] - epsvecn[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]) + epsvecn[i] / mma_ptr->machineEpsilon;
               diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + mma_ptr->xsi[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]) + mma_ptr->eta[i] / mma_ptr->machineEpsilon;
            }
            else
            {
               delx[i] = dpsidx[i] - epsvecn[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]) + epsvecn[i] / (mma_ptr->beta[i] - mma_ptr->x[i]);
               diagx[i] = 2 * (plam[i] / (ux1[i] * ux1[i] * ux1[i]) + qlam[i] / (xl1[i] * xl1[i] * xl1[i])) + mma_ptr->xsi[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]) + mma_ptr->eta[i] / (mma_ptr->beta[i] - mma_ptr->x[i]);
            }
         }

         delz = mma_ptr->a0 - epsi / mma_ptr->z;
         for (int i = 0; i < mma_ptr->nCon; i++)
         {
            gvec[i] = 0.0;
            // gvec = P/ux + Q/xl
            for (int j = 0; j < mma_ptr->nVar; j++)
            {
               gvec[i] += mma_ptr->P[i * mma_ptr->nVar + j] / ux1[j] + mma_ptr->Q[i * mma_ptr->nVar + j] / xl1[j];
               Puxinv[i * mma_ptr->nVar + j] = mma_ptr->P[i * mma_ptr->nVar + j] / (ux1[j] * ux1[j]);
               Qxlinv[i * mma_ptr->nVar + j] = mma_ptr->Q[i * mma_ptr->nVar + j] / (xl1[j] * xl1[j]);

               GG[i * mma_ptr->nVar + j] = Puxinv[i * mma_ptr->nVar + j] - Qxlinv[i * mma_ptr->nVar + j];
            }

            dely[i] = mma_ptr->c[i] + mma_ptr->d[i] * mma_ptr->y[i] - mma_ptr->lam[i] - epsvecm[i] / mma_ptr->y[i];
            delz -= mma_ptr->a[i] * mma_ptr->lam[i];
            dellam[i] = gvec[i] - mma_ptr->a[i] * mma_ptr->z - mma_ptr->y[i] - mma_ptr->b[i] + epsvecm[i] / mma_ptr->lam[i];
            diagy[i] = mma_ptr->d[i] + mma_ptr->mu[i] / mma_ptr->y[i];
            //diagyinv[i] = 1.0 / diagy[i];
            diaglam[i] = mma_ptr->s[i] / mma_ptr->lam[i];
            diaglamyi[i] = diaglam[i] + 1.0 / diagy[i];
         }

         if (mma_ptr->nCon < mma_ptr->nVar)
         {
            // blam = dellam + dely./diagy - GG*(delx./diagx);
            // bb = [blam; delz];
            for (int j = 0; j < mma_ptr->nCon; j++)
            {
               sum = 0.0;
               for (int i = 0; i < mma_ptr->nVar; i++)
               {
                  sum += GG[j * mma_ptr->nVar + i] * (delx[i] / diagx[i]);
               }
               blam[j] = dellam[j] + dely[j] / diagy[j] - sum;
               bb1[j] = blam[j];
            }
            bb1[mma_ptr->nCon] = delz;            

            // Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               // Axx = GG*spdiags(diagxinv,0,n,n);
               for (int k = 0; k < mma_ptr->nVar; k++)
               {
                  Axx[i * mma_ptr->nVar + k] = GG[k * mma_ptr->nCon + i] / diagx[k];
               }
            }
            // Alam = spdiags(diaglamyi,0,m,m) + Axx*GG';
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               for (int j = 0; j < mma_ptr->nCon; j++)
               {
                  Alam[i * mma_ptr->nCon + j] = 0.0;
                  for (int k = 0; k < mma_ptr->nVar; k++)
                  {
                     Alam[i * mma_ptr->nCon + j] += Axx[i * mma_ptr->nVar + k] * GG[j * mma_ptr->nVar + k];
                  }
                  if (i == j)
                  {
                     Alam[i * mma_ptr->nCon + j] += diaglamyi[i];
                  }
               }
            }
            // AA1 = [Alam     a
            //       a'    -zet/z];
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               for (int j = 0; j < mma_ptr->nCon; j++)
               {
                  AA1[i * mma_ptr->nCon + j] = Alam[i * mma_ptr->nCon + j];
               }
               AA1[mma_ptr->nCon * mma_ptr->nCon + i] = mma_ptr->a[i];
               AA1[mma_ptr->nCon * mma_ptr->nCon + mma_ptr->nCon + i] = mma_ptr->a[i];
            }
            AA1[mma_ptr->nCon * mma_ptr->nCon + mma_ptr->nCon + mma_ptr->nCon] = -mma_ptr->zet / mma_ptr->z;           

            // ----------------------------------------------------------------------------
            //solut = AA1\bb1 --> solve linear system of equations using LAPACK
            int info;
            int nLAP = mma_ptr->nCon + 1;
            int nrhs = 1;
            int lda = nLAP;
            int ldb = nLAP;
            int* ipiv = new int[nLAP];
            dgesv_(&nLAP, &nrhs, AA1, &lda, ipiv, bb1, &ldb, &info);
            delete[] ipiv;
            for (int i = 0; i < (mma_ptr->nCon + 1); i++)
            {
               solut[i] = bb1[i];
            }
            // ----------------------------------------------------------------------------
            //dlam = solut(1:nCon);
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               dlam[i] = solut[i];
            }
            dz = solut[mma_ptr->nCon];
            //dx = -(GG'*dlam)./diagx - delx./diagx;
            for (int i = 0; i < mma_ptr->nVar; i++)
            {
               sum = 0.0;
               for (int j = 0; j < mma_ptr->nCon; j++)
               {
                  sum += GG[j * mma_ptr->nVar + i] * dlam[j];
               }
               dx[i] = -sum / diagx[i] - delx[i] / diagx[i];           
            }
         }
         else
         {
            sum = 0.0;
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               //diaglamyiinv[i] = 1.0 / diaglamyi[i];
               dellamyi[i] = dellam[i] + dely[i] / diagy[i];
               // azz = zet/z + a'*(a./diaglamyi)
               sum += mma_ptr->a[i] * (mma_ptr->a[i] / diaglamyi[i]);
            }
            azz = mma_ptr->zet / mma_ptr->z + sum;
            // Axx = spdiags(diagx,0,nVar,nVar) + GG'*spdiags(diaglamyiinv,0,nCon,nCon)*GG;
            // AA = [Axx      axz
            //       axz'     azz];
            for (int i = 0; i < mma_ptr->nVar; i++)
            {
               // Axx =  GG'*spdiags(diaglamyiinv,0,nCon,nCon);
               for (int k = 0; k < mma_ptr->nCon; k++)
               {
                  Axx[i * mma_ptr->nCon + k] = GG[k * mma_ptr->nVar + i] / diaglamyi[k];
               }
               sum = 0.0;
               // axz = -GG'*(a./diaglamyi)
               for (int j = 0; j < mma_ptr->nCon; j++)
               {
                  sum -= GG[j * mma_ptr->nVar + i] * (mma_ptr->a[j] / diaglamyi[j]);
               }
               axz[i] = sum;
            }
            //Assemble matrix AA
            for (int i = 0; i < (mma_ptr->nVar + 1); i++)
            {
               for (int j = 0; j < (mma_ptr->nVar + 1); j++)
               {
                  // AA = [Axx  .
                  //       .    .]
                  AA[i * (mma_ptr->nVar + 1) + j] = 0.0;
                  if (i < mma_ptr->nVar && j < mma_ptr->nVar)
                  {
                     // Axx = Axx*GG
                     for (int k = 0; k < mma_ptr->nCon; k++)
                     {
                        AA[i * (mma_ptr->nVar + 1) + j] += Axx[i * mma_ptr->nCon + k] * GG[k * mma_ptr->nVar + j];
                     }
                     // Axx = Axx + spdiags(diagx,0,nVar,nVar)
                     if (i == j)
                     {
                        AA[i * (mma_ptr->nVar + 1) + j] += diagx[j];
                     }
                  }
                  // AA = [Axx  axz
                  //       axz' azz]
                  else if (i < mma_ptr->nVar && j == mma_ptr->nVar)
                  {
                     AA[i * (mma_ptr->nVar + 1) + j] = axz[i];
                  }
                  else if (i == mma_ptr->nVar && j < mma_ptr->nVar)
                  {
                     AA[i * (mma_ptr->nVar + 1) + j] = axz[j];
                  }
                  else
                  {
                     AA[i * (mma_ptr->nVar + 1) + j] = azz;
                  }
               }
            }
            // bb = [-bx'; -bz]
            // bx = delx - GG'*(dellamyi./diaglamyi)
            // bz = delz - a'*(dellamyi./diaglamyi)
            for (int i = 0; i < mma_ptr->nVar; i++)
            {
               sum = 0.0;
               for (int j = 0; j < mma_ptr->nCon; j++)
               {
                  sum += GG[j * mma_ptr->nVar + i] * (dellamyi[j] / diaglamyi[j]);
               }
               bb[i] = -(delx[i] + sum);
            }
            sum = 0.0;
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               sum += mma_ptr->a[i] * (dellamyi[i] / diaglamyi[i]);
            }
            bb[mma_ptr->nVar] = -(delz - sum);
            // ----------------------------------------------------------------------------
            //solut = AA\bb --> solve linear system of equations using LAPACK
            int info;
            int nLAP = mma_ptr->nVar + 1;
            int nrhs = 1;
            int lda = nLAP;
            int ldb = nLAP;
            int* ipiv = new int[nLAP];
            dgesv_(&nLAP, &nrhs, AA, &lda, ipiv, bb, &ldb, &info);
            delete[] ipiv;
            for (int i = 0; i < (mma_ptr->nVar + 1); i++)
            {
               solut[i] = bb[i];
            }
            // ----------------------------------------------------------------------------
            //dx = solut(1:nVar);
            for (int i = 0; i < mma_ptr->nVar; i++)
            {
               dx[i] = solut[i];
            }
            dz = solut[mma_ptr->nVar];
            //dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               sum = 0.0;
               for (int j = 0; j < mma_ptr->nVar; j++)
               {
                  sum += GG[i * mma_ptr->nVar + j] * dx[j];
               }
               dlam[i] = sum / diaglamyi[i] - dz * (mma_ptr->a[i] / diaglamyi[i]) + dellamyi[i] / diaglamyi[i];
            }
         }

         for (int i = 0; i < mma_ptr->nCon; i++)
         {
            dy[i] = -dely[i] / diagy[i] + dlam[i] / diagy[i];
            dmu[i] = -mma_ptr->mu[i] + epsvecm[i] / mma_ptr->y[i] - (mma_ptr->mu[i] * dy[i]) / mma_ptr->y[i];
            ds[i] = -mma_ptr->s[i] + epsvecm[i] / mma_ptr->lam[i] - (mma_ptr->s[i] * dlam[i]) / mma_ptr->lam[i];
            // xx = [y z lam xsi eta mu zet s]
            // dxx = [dy dz dlam dxsi deta dmu dzet ds]
            xx[i] = mma_ptr->y[i];
            xx[mma_ptr->nCon + 1 + i] = mma_ptr->lam[i];
            xx[2 * mma_ptr->nCon + 1 + 2 * mma_ptr->nVar + i] = mma_ptr->mu[i];
            xx[3 * mma_ptr->nCon + 2 * mma_ptr->nVar + 2 + i] = mma_ptr->s[i];
            dxx[i] = dy[i];
            dxx[mma_ptr->nCon + 1 + i] = dlam[i];
            dxx[2 * mma_ptr->nCon + 1 + 2 * mma_ptr->nVar + i] = dmu[i];
            dxx[3 * mma_ptr->nCon + 2 * mma_ptr->nVar + 2 + i] = ds[i];
         }
         xx[mma_ptr->nCon] = mma_ptr->z;
         xx[3 * mma_ptr->nCon + 2 * mma_ptr->nVar + 1] = mma_ptr->zet;
         dxx[mma_ptr->nCon] = dz;
         dxx[3 * mma_ptr->nCon + 2 * mma_ptr->nVar + 1] = dzet;
         for (int i = 0; i < mma_ptr->nVar; i++)
         {
            // NaN-Avoidance
            if (std::fabs(mma_ptr->x[i] - mma_ptr->alfa[i]) < mma_ptr->machineEpsilon)
            {
               if (std::fabs(mma_ptr->beta[i] - mma_ptr->x[i]) < mma_ptr->machineEpsilon)
               {
                  dxsi[i] = -mma_ptr->xsi[i] + epsvecn[i] / mma_ptr->machineEpsilon - (mma_ptr->xsi[i] * dx[i]) / mma_ptr->machineEpsilon;
                  deta[i] = -mma_ptr->eta[i] + epsvecn[i] / mma_ptr->machineEpsilon + (mma_ptr->eta[i] * dx[i]) / mma_ptr->machineEpsilon;
               }
               else
               {
                  dxsi[i] = -mma_ptr->xsi[i] + epsvecn[i] / mma_ptr->machineEpsilon - (mma_ptr->xsi[i] * dx[i]) / mma_ptr->machineEpsilon;
                  deta[i] = -mma_ptr->eta[i] + epsvecn[i] / (mma_ptr->beta[i] - mma_ptr->x[i]) + (mma_ptr->eta[i] * dx[i]) / (mma_ptr->beta[i] - mma_ptr->x[i]);
               }
            }
            else if (std::fabs(mma_ptr->beta[i] - mma_ptr->x[i]) < mma_ptr->machineEpsilon)
            {
               dxsi[i] = -mma_ptr->xsi[i] + epsvecn[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]) - (mma_ptr->xsi[i] * dx[i]) / (mma_ptr->x[i] - mma_ptr->alfa[i]);
               deta[i] = -mma_ptr->eta[i] + epsvecn[i] / mma_ptr->machineEpsilon + (mma_ptr->eta[i] * dx[i]) / mma_ptr->machineEpsilon;
            }
            else
            {
               dxsi[i] = -mma_ptr->xsi[i] + epsvecn[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]) - (mma_ptr->xsi[i] * dx[i]) / (mma_ptr->x[i] - mma_ptr->alfa[i]);
               deta[i] = -mma_ptr->eta[i] + epsvecn[i] / (mma_ptr->beta[i] - mma_ptr->x[i]) + (mma_ptr->eta[i] * dx[i]) / (mma_ptr->beta[i] - mma_ptr->x[i]);
            }
            xx[mma_ptr->nCon + 1 + mma_ptr->nCon + i] = mma_ptr->xsi[i];
            xx[mma_ptr->nCon + 1 + mma_ptr->nCon + mma_ptr->nVar + i] = mma_ptr->eta[i];
            dxx[mma_ptr->nCon + 1 + mma_ptr->nCon + i] = dxsi[i];
            dxx[mma_ptr->nCon + 1 + mma_ptr->nCon + mma_ptr->nVar + i] = deta[i];
         }
         dzet = -mma_ptr->zet + epsi / mma_ptr->z - mma_ptr->zet * dz / mma_ptr->z;

         stmxx = 0.0;
         for (int i = 0; i < (4 * mma_ptr->nCon + 2 * mma_ptr->nVar + 2); i++)
         {
            stepxx[i] = -1.01*dxx[i] / xx[i];
            if (stepxx[i] > stmxx)
            {
               stmxx = stepxx[i];
            }
         }
         stmalfa = 0.0;
         stmbeta = 0.0;
         for (int i = 0; i < mma_ptr->nVar; i++)
         {
            
            //NaN-Avoidance
            if (std::fabs(mma_ptr->x[i] - mma_ptr->alfa[i]) < mma_ptr->machineEpsilon)
            {
               stepalfa[i] = -1.01*dx[i] / mma_ptr->machineEpsilon;
            }
            else 
            {
               stepalfa[i] = -1.01*dx[i] / (mma_ptr->x[i] - mma_ptr->alfa[i]);
            }
            if (std::fabs(mma_ptr->beta[i] - mma_ptr->x[i]) < mma_ptr->machineEpsilon)
            {
               stepbeta[i] = 1.01*dx[i] / mma_ptr->machineEpsilon;
            }
            else
            {
               stepbeta[i] = 1.01*dx[i] / (mma_ptr->beta[i] - mma_ptr->x[i]);
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

         for (int i = 0; i < mma_ptr->nVar; i++)
         {
            xold[i] = mma_ptr->x[i];
            xsiold[i] = mma_ptr->xsi[i];
            etaold[i] = mma_ptr->eta[i];
         }
         for (int i = 0; i < mma_ptr->nCon; i++)
         {
            yold[i] = mma_ptr->y[i];
            lamold[i] = mma_ptr->lam[i];
            muold[i] = mma_ptr->mu[i];
            sold[i] = mma_ptr->s[i];
         }
         zold = mma_ptr->z;
         zetold = mma_ptr->zet;

         itto = 0;
         resinew = 2.0 * residunorm;
         while (resinew > residunorm && itto < 50)
         {
            itto++;

            for (int i = 0; i < mma_ptr->nCon; ++i)
            {
               mma_ptr->y[i] = yold[i] + steg * dy[i];
               if (mma_ptr->y[i] == 0)
               {
                  mma_ptr->y[i] = mma_ptr->machineEpsilon;
               }
               else
               {
                  while (std::fabs(mma_ptr->y[i]) < mma_ptr->machineEpsilon)
                  {
                     mma_ptr->y[i] *= 10;
                  }
               }
               mma_ptr->lam[i] = lamold[i] + steg * dlam[i];
               if (mma_ptr->lam[i] == 0)
               {
                  mma_ptr->lam[i] = mma_ptr->machineEpsilon;
               }
               else
               {
                  while (std::fabs(mma_ptr->lam[i]) < mma_ptr->machineEpsilon)
                  {
                     mma_ptr->lam[i] *= 10;
                  }
               }
               mma_ptr->mu[i] = muold[i] + steg * dmu[i];
               mma_ptr->s[i] = sold[i] + steg * ds[i];
            }

            for (int i = 0; i < mma_ptr->nVar; ++i)
            {
               mma_ptr->x[i] = xold[i] + steg * dx[i];
               mma_ptr->xsi[i] = xsiold[i] + steg * dxsi[i];
               mma_ptr->eta[i] = etaold[i] + steg * deta[i];
               if (std::isnan(mma_ptr->eta[i]))
               {
                  mma_ptr->eta[i] = mma_ptr->machineEpsilon;
               }
               ux1[i] = mma_ptr->upp[i] - mma_ptr->x[i];
               if (ux1[i] == 0)
               {
                  ux1[i] = mma_ptr->machineEpsilon;
               }
               else
               {
                  while (std::fabs(ux1[i]) <= mma_ptr->machineEpsilon)
                  {
                     ux1[i] *= 10;
                  }
               }
               xl1[i] = mma_ptr->x[i] - mma_ptr->low[i];
               if (xl1[i] == 0)
               {
                  xl1[i] = mma_ptr->machineEpsilon;
               }
               else
               {
                  while (std::fabs(xl1[i]) <= mma_ptr->machineEpsilon)
                  {
                     xl1[i] *= 10;
                  }
               }
               // plam & qlam
               plam[i] = mma_ptr->p0[i];
               qlam[i] = mma_ptr->q0[i];
               for (int j = 0; j < mma_ptr->nCon; j++)
               {
                  plam[i] += mma_ptr->P[j * mma_ptr->nVar + i] * mma_ptr->lam[j];
                  qlam[i] += mma_ptr->Q[j * mma_ptr->nVar + i] * mma_ptr->lam[j];
               }
               dpsidx[i] = plam[i] / (ux1[i] * ux1[i]) - qlam[i] / (xl1[i] * xl1[i]);
               rex[i] = dpsidx[i] - mma_ptr->xsi[i] + mma_ptr->eta[i];
               rez -= mma_ptr->a[i] * mma_ptr->lam[i];
               rexsi[i] = mma_ptr->xsi[i] * (mma_ptr->x[i] - mma_ptr->alfa[i]) - epsvecn[i];              
               if (std::fabs(mma_ptr->x[i] - mma_ptr->alfa[i]) < mma_ptr->machineEpsilon || std::isnan(rexsi[i]))
               {
                  rexsi[i] = mma_ptr->xsi[i] * mma_ptr->machineEpsilon - epsvecn[i]; 
               }
               reeta[i] = mma_ptr->eta[i] * (mma_ptr->beta[i] - mma_ptr->x[i]) - epsvecn[i];
               if (std::fabs(mma_ptr->beta[i] - mma_ptr->x[i]) < mma_ptr->machineEpsilon || std::isnan(reeta[i]))
               {
                  reeta[i] = mma_ptr->eta[i] * mma_ptr->machineEpsilon - epsvecn[i];
               }
               residu1[i] = rex[i];
               residu2[mma_ptr->nCon + i] = rexsi[i];
               residu2[mma_ptr->nCon + mma_ptr->nVar + i] = reeta[i];
            }
            mma_ptr->z = zold + steg * dz;
            if (mma_ptr->z == 0)
            {
               mma_ptr->z = mma_ptr->machineEpsilon;
            }
            else
            {
               while (std::fabs(mma_ptr->z) <= mma_ptr->machineEpsilon)
               {
                  mma_ptr->z *= 10;
               }
            }
            mma_ptr->zet = zetold + steg * dzet;

            // gvec = P/ux + Q/xl
            for (int i = 0; i < mma_ptr->nCon; i++)
            {
               gvec[i] = 0.0;
               for (int j = 0; j < mma_ptr->nVar; j++)
               {
                  gvec[i] += mma_ptr->P[i * mma_ptr->nVar + j] / ux1[j] + mma_ptr->Q[i * mma_ptr->nVar + j] / xl1[j];
               }
               rey[i] = mma_ptr->c[i] + mma_ptr->d[i] * mma_ptr->y[i] - mma_ptr->mu[i] - mma_ptr->lam[i];
               relam[i] = gvec[i] - mma_ptr->a[i] * mma_ptr->z - mma_ptr->y[i] + mma_ptr->s[i] - mma_ptr->b[i];
               remu[i] = mma_ptr->mu[i] * mma_ptr->y[i] - epsvecm[i];
               res[i] = mma_ptr->lam[i] * mma_ptr->s[i] - epsvecm[i];

               residu2[i] = relam[i];
               residu1[mma_ptr->nVar + i] = rey[i];
               residu2[mma_ptr->nCon + 2 * mma_ptr->nVar + i] = remu[i];
               residu2[2 * mma_ptr->nVar + 2 * mma_ptr->nCon + 1 + i] = res[i];
            }
            rez = mma_ptr->a0 - mma_ptr->zet;
            rezet = mma_ptr->zet * mma_ptr->z - epsi;
            residu1[mma_ptr->nVar + mma_ptr->nCon] = rez;
            residu2[2 * mma_ptr->nVar + 2 * mma_ptr->nCon] = rezet;

            // Concatenate the residuals
            for (int i = 0; i < (mma_ptr->nVar + mma_ptr->nCon + 1); i++)
            {
               residu[i] = residu1[i];
            }
            for (int i = 0; i < (2 * mma_ptr->nVar + 3 * mma_ptr->nCon + 1); i++)
            {
               residu[mma_ptr->nVar + mma_ptr->nCon + 1 + i] = residu2[i];
            }
            // New residual
            sum = 0.0;
            for (int i = 0; i < (3 * mma_ptr->nVar + 4 * mma_ptr->nCon + 2); i++)
            {
               sum += residu[i] * residu[i];
            }
            // Norm of the residual
            resinew = std::sqrt(sum);
            steg = steg / 2.0;
         }
         residunorm = resinew;
         residumax = 0.0;
         for (int i = 0; i < (3 * mma_ptr->nVar + 4 * mma_ptr->nCon + 2); i++)
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
   // should return x, y, z, lam, xsi, eta, mu, zet, s
   //kktnorm = kktcheck(mma_ptr->y, dfdx, gx, dgdx, mma_ptr->x);
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

double MMA::SubProblemClassic::kktcheck(double* y, double* const dfdx, double* const gx, double* const dgdx, double* x)
{
   //std::ofstream kkt;
   //kkt.open("KKT.dat", std::ios::app);

   for (int i = 0; i < mma_ptr->nVar; i++)
   {
      sum1[i] = 0.0;
      for (int j = 0; j < mma_ptr->nCon; j++)
      {
         sum1[i] += dgdx[j * mma_ptr->nVar + i] * mma_ptr->lam[j];
      }
   }
   for (int i = 0; i < mma_ptr->nVar; i++)
   {
      rex[i] = dfdx[i] + sum1[i] - mma_ptr->xsi[i] + mma_ptr->eta[i];
      rexsi[i] = mma_ptr->xsi[i] * (mma_ptr->x[i] - mma_ptr->xmin[i]);      
      reeta[i] = mma_ptr->eta[i] * (mma_ptr->xmax[i] - mma_ptr->x[i]);
   }
   for (int i = 0; i < mma_ptr->nCon; i++)
   {
      rey[i] = mma_ptr->c[i] + mma_ptr->d[i] * mma_ptr->y[i] - mma_ptr->mu[i] - mma_ptr->lam[i];
      relam[i] = gx[i] - mma_ptr->a[i] * mma_ptr->z - mma_ptr->y[i] + mma_ptr->s[i];
      remu[i] = mma_ptr->mu[i] * mma_ptr->y[i];
      res[i] = mma_ptr->lam[i] * mma_ptr->s[i];
   }
   sum = 0.0;
   for (int i = 0; i < mma_ptr->nCon; i++)
   {
      sum += mma_ptr->a[i] * mma_ptr->lam[i];
   }
   rez = mma_ptr->a0 - mma_ptr->zet - sum;
   rezet = mma_ptr->zet * mma_ptr->z;
   //---------------------------------------------------------------------

   for (int i = 0; i < mma_ptr->nVar; i++)
   {
      residu1[i] = rex[i];
      residu2[mma_ptr->nCon + i] = rexsi[i];
      residu2[mma_ptr->nCon + mma_ptr->nVar + i] = reeta[i];
   }

   for (int i = 0; i < mma_ptr->nCon; i++)
   {
      residu2[i] = relam[i];
      residu1[mma_ptr->nVar + i] = rey[i];
      residu2[mma_ptr->nCon + 2 * mma_ptr->nVar + i] = remu[i];
      residu2[2 * mma_ptr->nVar + 2 * mma_ptr->nCon + 1 + i] = res[i];
   }
   residu1[mma_ptr->nVar + mma_ptr->nCon] = rez;
   residu2[2 * mma_ptr->nVar + 2 * mma_ptr->nCon] = rezet;

   // Concatenate the residuals
   for (int i = 0; i < (mma_ptr->nVar + mma_ptr->nCon + 1); i++)
   {
      residu[i] = residu1[i];
   }
   for (int i = 0; i < (2 * mma_ptr->nVar + 3 * mma_ptr->nCon + 1); i++)
   {
      residu[mma_ptr->nVar + mma_ptr->nCon + 1 + i] = residu2[i];
   }
   //Get vector product and maximum absolute value
   residunorm = 0.0;
   residumax = 0.0;
   for (int i = 0; i < (3 * mma_ptr->nVar + 4 * mma_ptr->nCon + 2); i++)
   {
      residunorm += residu[i] * residu[i];
      residumax = std::max(residumax, std::abs(residu[i]));
   }
   // Norm of the residual
   residunorm = std::sqrt(residunorm);
   mma_ptr->kktnorm = residunorm;

   return mma_ptr->kktnorm;
   
   //kkt.close();
}

void MMA::SubProblemClassic::setSubProb(int nVar, int nCon)
{
   epsi = 1.0;
   epsvecn = new double[nVar];
   epsvecm = new double[nCon];
   ittt = itto = itera = 0;
   ux1 = new double[nVar];
   xl1 = new double[nVar];
   plam = new double[nVar];
   qlam = new double[nVar];
   gvec = new double[nCon];
   dpsidx = new double[nVar];
   rex = new double[nVar];
   rey = new double[nCon];
   relam = new double[nCon];
   rexsi = new double[nVar];
   reeta = new double[nVar];
   remu = new double[nCon];
   res = new double[nCon];
   residu1 = new double[nVar + nCon + 1];
   residu2 = new double[3 * nCon + 2 * nVar + 1];
   residu = new double[3 * nVar + 4 * nCon + 2];
   GG = new double[nVar * nCon];
   Puxinv = new double[nVar * nCon];
   Qxlinv = new double[nVar * nCon];
   delx = new double[nVar];
   dely = new double[nCon];
   dellam = new double[nCon];
   dellamyi = new double[nCon];
   diagx = new double[nVar];
   diagy = new double[nCon];
   diaglam = new double[nCon];
   diaglamyi = new double[nCon];
   blam = new double[nCon];
   bb = new double[nVar + 1];
   bb1 = new double[nCon + 1];
   Alam = new double[nCon * nCon];
   AA = new double[(nVar + 1) * (nVar + 1)];
   AA1 = new double[(nCon + 1) * (nCon + 1)];
   solut = new double[nVar + nCon + 1];
   dlam = new double[nCon];
   dx = new double[nVar];
   dy = new double[nCon];
   dxsi = new double[nVar];
   deta = new double[nVar];
   dmu = new double[nCon];
   Axx = new double[nVar * nCon];
   axz = new double[nVar];
   ds = new double[nCon];
   xx = new double[4 * nCon + 2 * nVar + 2];
   dxx = new double[4 * nCon + 2 * nVar + 2];
   stepxx = new double[4 * nCon + 2 * nVar + 2];
   sum = 0;
   sum1 = new double[nVar];
   stepalfa = new double[nVar];
   stepbeta = new double[nVar];
   xold = new double[nVar];
   yold = new double[nCon];
   lamold = new double[nCon];
   xsiold = new double[nVar];
   etaold = new double[nVar];
   muold = new double[nCon];
   sold = new double[nCon];
}

void MMA::SubProblemClassic::freeSubProb()
{
   delete[] sum1;
   delete[] epsvecn;
   delete[] epsvecm;
   delete[] ux1;
   delete[] xl1;
   delete[] plam;
   delete[] qlam;
   delete[] gvec;
   delete[] dpsidx;
   delete[] rex;
   delete[] rey;
   delete[] relam;
   delete[] rexsi;
   delete[] reeta;
   delete[] remu;
   delete[] res;
   delete[] residu1;
   delete[] residu2;
   delete[] residu;
   delete[] GG;
   delete[] Puxinv;
   delete[] Qxlinv;
   delete[] delx;
   delete[] dely;
   delete[] dellam;
   delete[] dellamyi;
   delete[] diagx;
   delete[] diagy;
   delete[] diaglam;
   delete[] diaglamyi;
   delete[] blam;
   delete[] bb;
   delete[] bb1;
   delete[] Alam;
   delete[] AA;
   delete[] AA1;
   delete[] solut;
   delete[] dlam;
   delete[] dx;
   delete[] dy;
   delete[] dxsi;
   delete[] deta;
   delete[] dmu;
   delete[] Axx;
   delete[] axz;
   delete[] ds;
   delete[] xx;
   delete[] dxx;
   delete[] stepxx;
   delete[] stepalfa;
   delete[] stepbeta;
   delete[] xold;
   delete[] yold;
   delete[] lamold;
   delete[] xsiold;
   delete[] etaold;
   delete[] muold;
   delete[] sold;
}

void MMA::SubProblemClassicMPI::Perform()
{
   //do stuff
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

void MMA::setGlobals(int nVar, int nCon)
{
   this->nVar = nVar;
   this->nCon = nCon;
   xmin = new double[nVar];
   xmax = new double[nVar];
   low = new double[nVar];
   upp = new double[nVar];
   xo1 = new double[nVar];
   xo2 = new double[nVar];
   x = new double[nVar];
   y = new double[nCon];
   c = new double[nCon];
   d = new double[nCon];
   a = new double[nCon];
   lam = new double[nCon];
   xsi = new double[nVar];
   eta = new double[nVar];
   mu = new double[nCon];
   s = new double[nCon];
   q0 = new double[nVar];
   p0 = new double[nVar];
   P = new double[nCon * nVar];
   Q = new double[nCon * nVar];
   alfa = new double[nVar];
   beta = new double[nVar];
   z = zet = 1.0;
   kktnorm = 10;
   machineEpsilon = 1e-10;

   isInitialized = true;
} 

void MMA::freeGlobals()
{
   delete[] xmin;
   delete[] xmax;
   delete[] low;
   delete[] upp;
   delete[] xo1;
   delete[] xo2;
   delete[] x;
   delete[] y;
   delete[] c;
   delete[] d;
   delete[] a;
   delete[] lam;
   delete[] xsi;
   delete[] eta;
   delete[] mu;
   delete[] s;
   delete[] q0;
   delete[] p0;
   delete[] P;
   delete[] Q;
   isInitialized = false;
}

void MMA::setMMA(int nVar, int nCon)
{
   epsimin = 1e-7;
   raa0 = 0.00001;
   move = 0.5;
   albefa = 0.1;
   asyinit = 0.5;
   asyincr = 1.2;
   asydecr = 0.7;
   xmamieps = 1e-5;
   factor = new double[nVar];
   lowmin = lowmax = uppmin = uppmax = zz = 0.0;
   xmami = new double[nVar];
   pq0 = new double[nVar];
   p = new double[nVar * nCon];
   q = new double[nVar * nCon];
   pq = new double[nVar * nCon];
   b = new double[nCon];
   PQ = new double[nCon * nVar];
   ux1 = new double[nVar];
   xl1 = new double[nVar];
}

void MMA::freeMMA()
{
   delete[] factor;
   delete[] xmami;
   delete[] pq0;
   delete[] p;
   delete[] q;
   delete[] pq;
   delete[] b;
   delete[] PQ;
   delete[] ux1;
   delete[] xl1;
}

void MMA::subProblemFactory(int nVar, int nCon, enum SubProblemType type)
{
   switch (type)
   {
      case CLASSIC:
         mSubProblem = new SubProblemClassic(this, nVar, nCon);
         break;
      case MPIClassic:
         mSubProblem = new SubProblemClassicMPI(this, nVar, nCon);
         //reinterpret_cast<SubProblemClassicMPI*>(mSubProblem)->Set.);
         break;
      default:
         mfem::mfem_error("MMA::subProblemFactory() invalid subproblem type");
         break;
   }
}

}
// end mma namespace
