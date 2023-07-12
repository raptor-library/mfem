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

namespace mma
{

MMA::MMA(int n, int m, double * x, int sx)
{

}

//MMA::MMA(MPI_Comm Comm,int n, int m, double * x, int sx) {}

void MMA::mmasub(int nVar, int nCon, int iter, double* xval, double* xmin,
                 double* xmax, double* xold1, double* xold2, double fval,
                 double* dfdx, double* gx, double* dgdx, double* low,
                 double* upp, double a0, double* a, double* c, double* d,
                 double* xmma, double* ymma, double* zmma, double* lam,
                 double* xsi, double* eta, double* mu, double& zet, double* s)
{
   double epsimin = 1e-7;
   double raa0 = 0.00001;
   double move = 0.5;
   double albefa = 0.1;
   double asyinit = 0.5;
   double asyincr = 1.2;
   double asydecr = 0.7;
   double xmamieps = 1e-5;

   double* eeen = new double[nVar];
   double* eeem = new double[nCon];
   double* zeron = new double[nVar];

   double* factor = new double[nVar];
   double lowmin, lowmax, uppmin, uppmax, z, zzz, zzz1, zzz2;
   double* xmami = new double(nVar);
   // double* xmamiinv = new double(nVar);
   double* ux1 = new double(nVar);
   double* ux2 = new double(nVar);
   double* xl1 = new double(nVar);
   double* xl2 = new double(nVar);
   double* uxinv = new double(nVar);
   double* xlinv = new double(nVar);
   double* p0 = new double(nVar);
   double* q0 = new double(nVar);
   double* p = new double(nCon * nVar);
   double* q = new double(nCon * nVar);
   double* b = new double(nCon);
   double* alfa = new double(nVar);
   double* beta = new double(nVar);
   double* pq0 = new double(nVar);

   double* P = new double(nCon * nVar);
   double* Q = new double(nCon * nVar);

   for (int i = 0; i < nVar; i++)
   {
      eeen[i] = 1.0;
      zeron[i] = 0.0;
      factor[i] = asyincr;
      p0[i] = 0.0;
      q0[i] = 0.0;
   }

   for (int i = 0; i < nCon; i++)
   {
      eeem[i] = 1.0;
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
      //Determine sign
      for (int i = 0; i < nVar; i++)
      {
         z = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]);
         if (z > 0.0)
         {
            factor[i] = asyincr;
         }
         else if (z < 0.0)
         {
            factor[i] = asydecr;
         }
      }
      //Find new asymptote
      for (int i = 0; i < nVar; ++i)
      {
         low[i] = xval[i] - factor[i] * (xold1[i] - low[i]);
         upp[i] = xval[i] + factor[i] * (upp[i] - xold1[i]);

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

   // Calculation of bounds alfa and beta according to:
   // alfa = max{xmin, low + 0.1(xval-low), xval-0.5(xmax-xmin)}
   // beta = min{xmax, upp - 0.1(upp-xval), xval+0.5(xmax-xmin)}

   for (int i = 0; i < nVar; i++)
   {
      zzz1 = low[i] + albefa * (xval[i] - low[i]);
      zzz2 = xval[i] - move * (xmax[i] - xmin[i]);
      zzz = std::max(zzz1, zzz2);
      alfa[i] = std::max(zzz, xmin[i]);
      zzz1 = upp[i] - albefa * (upp[i] - xval[i]);
      zzz2 = xval[i] + move * (xmax[i] - xmin[i]);
      zzz = std::min(zzz1, zzz2);
      beta[i] = std::min(zzz, xmax[i]);
   }

   for (int i = 0; i < nVar; i++)
   {
      xmami[i] = std::max(xmax[i] - xmin[i], xmamieps);
   }

   // Calculations of p0, q0, P, Q, and b
   for (int i = 0; i < nVar; i++)
   {
      ux1[i] = upp[i] - xval[i];
      ux2[i] = ux1[i] * ux1[i];
      xl1[i] = xval[i] - low[i];
      xl2[i] = xl1[i] * xl1[i];
      uxinv[i] = eeen[i] / ux1[i];
      xlinv[i] = eeen[i] / xl1[i];
   }

   for (int i = 0; i < nVar; ++i)
   {
      p0[i] = std::max(dfdx[i], 0.0);
      q0[i] = std::max(-dfdx[i], 0.0);
      pq0[i] = 0.001 * (p0[i] + q0[i]) + raa0 / xmami[i];
      p0[i] += pq0[i];
      q0[i] += pq0[i];
      p0[i] *= ux2[i];
      q0[i] *= xl2[i];
   }

   for (int i = 0; i < nCon; ++i)
   {
      for (int j = 0; j < nVar; ++j)
      {
         p[i * nVar + j] = std::max(dgdx[i * nVar + j], 0.0);
         q[i * nVar + j] = std::max(-dgdx[i * nVar + j], 0.0);
      }
   }
   for (int i = 0; i < nCon; ++i)
   {
      for (int j = 0; j < nVar; ++j)
      {
         P[i * nVar + j] = p[i * nVar + j] * ux2[j];
         Q[i * nVar + j] = q[i * nVar + j] * xl2[j];
      }
   }

   for (int i = 0; i < nCon; ++i)
   {
      for (int j = 0; j < nVar; ++j)
      {
         b[i] += P[i * nVar + j] * uxinv[j] + Q[i * nVar + j] * xlinv[j] - gx[i];
      }
   }

   subsolv(nVar, nCon, epsimin, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d,
           xmma, ymma, zmma, lam, xsi, eta, mu, &zet, s);
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
 * Output: xmma, ymma, zmma, slack variables and Lagrange multiplers.
*/
void MMA::subsolv(int nVar, int nCon, double epsimin, double* low, double* upp,
                  double* alfa, double* beta, double* p0,
                  double* q0, double* P, double* Q,
                  double a0, double* a, double* b, double* c,
                  double* d, double* xmma, double* ymma,
                  double* zmma, double* lamma, double* xsimma,
                  double* etamma, double* mumma, double* zetmma, double* smma)
{
   double* een = new double(nVar);
   double* eem = new double(nCon);
   double epsi = 1.0;
   double* epsvecn = new double(nVar);
   double* epsvecm = new double(nCon);
   double* x = new double(nVar);
   double* y = new double(nCon);
   double z = 1.0;
   double* lam = new double(nCon);
   double* xsi = new double(nVar);
   double* eta = new double(nVar);
   double* mu = new double(nCon);
   double zet = 1.0;
   double* s = new double(nCon);
   int ittt, itto, itera = 0;

   double* ux1 = new double(nVar);
   double* ux2 = new double(nVar);
   double* ux3 = new double(nVar);
   double* xl1 = new double(nVar);
   double* xl2 = new double(nVar);
   double* xl3 = new double(nVar);
   double* uxinv1 = new double(nVar);
   double* xlinv1 = new double(nVar);
   double* uxinv2 = new double(nVar);
   double* xlinv2 = new double(nVar);
   double* plam = new double(nVar);
   double* qlam = new double(nVar);
   double* gvec = new double(nCon);
   double* dpsidx = new double(nVar);
   double* rex = new double(nVar);
   double* rey = new double(nCon);
   double rez;
   double* relam = new double(nCon);
   double* rexsi = new double(nVar);
   double* reeta = new double(nVar);
   double* remu = new double(nCon);
   double rezet;
   double* res = new double(nCon);
   double* residu1 = new double(nVar + nCon + 1);
   double* residu2 = new double(3 * nCon + 2 * nVar + 1);
   double* residu = new double(3 * nVar + 4* nCon + 2);
   double residunorm, residumax, resinew;
   double* GG = new double(nVar * nCon);
   double* Puxinv = new double(nVar * nCon);
   double* Qxlinv = new double(nVar * nCon);
   double* delx = new double(nVar);
   double* dely = new double(nCon);
   double delz;
   double* dellam = new double(nCon);
   double* dellamyi = new double(nCon);
   double* diagx = new double(nVar);
   double* diagxinv = new double(nVar);
   double* diagy = new double(nCon);
   double* diagyinv = new double(nCon);
   double* diaglam = new double(nCon);
   double* diaglamyi = new double(nCon);
   double* diaglamyiinv = new double(nCon);
   double* blam = new double(nCon);
   double* bb = new double(nCon + 1);
   double* Alam = new double(nCon * nCon);
   double* AA = new double((nVar + nCon + 1) * (nVar + nCon + 1));
   double* solut = new double(nVar + nCon + 1);
   double* dlam = new double(nCon);
   double dz;
   double* dx = new double(nVar);
   double* dy = new double(nCon);
   double* dxsi = new double(nVar);
   double* deta = new double(nVar);
   double* dmu = new double(nCon);
   double dzet;
   double* azz = new double(nCon);
   double* axz = new double(nVar);
   double* ds = new double(nCon);
   double* xx = new double(4 * nCon + 2 * nVar + 2);
   double* dxx = new double(4 * nCon + 2 * nVar + 2);
   double* stepxx = new double(4 * nCon + 2 * nVar + 2);
   double stmxx;
   double* stepalfa = new double(nVar);
   double stmalfa;
   double* stepbeta = new double(nVar);
   double stmbeta;
   double stmalbe;
   double stmalbexx;
   double stminv;
   double steg;

   double* xold = new double(nVar);
   double* yold = new double(nCon);
   double zold;
   double* lamold = new double(nCon);
   double* xsiold = new double(nVar);
   double* etaold = new double(nVar);
   double* muold = new double(nCon);
   double zetold;
   double* sold = new double(nCon);

   for (int i = 0; i < nVar; i++)
   {
      een[i] = 1.0;
      epsvecn[i] = epsi;
      x[i] = 0.5 * (alfa[i] + beta[i]);
      xsi[i] = een[i]/(x[i] - alfa[i]);
      xsi[i] = std::max(xsi[i], een[i]);
      eta[i] = een[i]/(beta[i] - x[i]);
      eta[i] = std::max(eta[i], een[i]);
   }
   for (int i = 0; i < nCon; i++)
   {
      eem[i] = 1.0;
      epsvecm[i] = epsi;
      y[i] = 1.0;
      lam[i] = 1.0;
      mu[i] = std::max(eem[i], 0.5 * c[i]);
      s[i] = eem[i];
   }

   while (epsi > epsimin)
   {
      for (int i = 0; i < nVar; i++)
      {
         ux1[i] = upp[i] - x[i];
         xl1[i] = x[i] - low[i];
         ux2[i] = ux1[i]*ux1[i];
         xl2[i] = xl1[i]*xl1[i];
         uxinv1[i] = een[i]/ux1[i];
         xlinv1[i] = een[i]/xl1[i];
      }
      // Matrix operations
      for (int i = 0; i < nCon; i++)
      {
         plam[i] = p0[i];
         qlam[i] = q0[i];
         for (int j = 0; j < nVar; j++)
         {
            plam[i] += P[i * nVar + j] * lam[j];
            qlam[i] += Q[i * nVar + j] * lam[j];
            gvec[j] = P[i * nVar + j] * uxinv1[j] + Q[i * nVar + j] * xlinv1[j];
         }
      }

      rez = a0 - zet;
      for (int i = 0; i < nVar; i++)
      {
         dpsidx[i] = plam[i] / ux2[i] - qlam[i] / xl2[i];
         rex[i] = dpsidx[i] - xsi[i] + eta[i];
         rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
         rez -= a[i] * lam[i];
         relam[i] = gvec[i] - a[i] * z - y[i] + s[i] - b[i];
         rexsi[i] = xsi[i] * (x[i] - alfa[i]) - epsvecn[i];
         reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
      }

      for (int i = 0; i < nCon; i++)
      {
         remu[i] = mu[i] * y[i] - epsvecm[i];
         res[i] = lam[i] * s[i] - epsvecm[i];
      }
      rezet = zet * z - epsi;

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

      ittt = 0;

      while (residumax > 0.9 * epsi && ittt < 200)
      {
         ittt++;
         for (int i = 0; i < nVar; i++)
         {
            ux1[i] = upp[i] - x[i];
            xl1[i] = x[i] - low[i];
            ux2[i] = ux1[i];
            xl2[i] = xl1[i];
            ux3[i] = ux1[i];
            xl3[i] = xl1[i];
            uxinv1[i] = een[i]/ux1[i];
            xlinv1[i] = een[i]/xl1[i];
            uxinv2[i] = een[i]/ux2[i];
            xlinv2[i] = een[i]/xl2[i];
         }
         for (int i = 0; i < nCon; i++)
         {
            plam[i] = p0[i];
            qlam[i] = q0[i];
            for (int j = 0; j < nVar; j++)
            {
               plam[i] += P[i * nVar + j] * lam[j];
               qlam[i] += Q[i * nVar + j] * lam[j];
               gvec[j] = P[i * nVar + j] * uxinv1[j] + Q[i * nVar + j] * xlinv1[j];
            }
         }
         for (int i = 0; i < (nVar * nCon); i++)
         {
            Puxinv[i] = P[i] * uxinv2[i];
            Qxlinv[i] = Q[i] * xlinv2[i];
         }
         for (int i = 0; i < nCon; i++)
         {
            for (int j = 0; j < nVar; j++)
            {
               GG[i * nVar + j] = Puxinv[i * nVar + j] + Qxlinv[i * nVar + j];
            }
         }
         for (int i = 0; i < nVar; i++)
         {
            dpsidx[i] = plam[i] / ux2[i] - qlam[i] / xl2[i];
            delx[i] = dpsidx[i] - epsvecn[i] / (x[i] - alfa[i]) + epsvecn[nVar + i] /
                      (beta[i] - x[i]);
            dely[i] = c[i] + d[i] * y[i] - lam[i] - epsvecn[i] / y[i];
         }
         delz = a0 - epsi/z;
         for (int i = 0; i < nCon; i++)
         {
            delz -= a[i] * lam[i];
            dellam[i] = gvec[i] - a[i] * z - y[i] - b[i] + epsvecm[i] / lam[i];
            diagy[i] = d[i] + mu[i] / y[i];
            diagyinv[i] = een[i]/diagy[i];
            diaglam[i] = s[i] / lam[i];
            diaglamyi[i] = diaglam[i] + diagyinv[i];
         }
         for (int i = 0; i < nVar; i++)
         {
            diagx[i] = 2 * (plam[i] / ux3[i] + qlam[i] / xl3[i]) + xsi[i] /
                       (x[i] - alfa[i]) + eta[i] / (beta[i] - x[i]);
            diagxinv[i] = een[i]/diagx[i];
         }

         if (nCon < nVar)
         {
            //To-Do: Double check this in debug
            // blam = dellam + dely./diagy - GG*(delx./diagx);
            // bb = [blam; delz];
            double sum = 0.0;
            for (int j = 0; j < nCon; j++)
            {
               sum = 0.0;
               for (int i = 0; i < nVar; i++)
               {
                  sum -= GG[j * nVar + i] * (delx[i] / diagx[i]);
               }
               blam[j] = dellam[j] + dely[j] / diagy[j] - sum;
               bb[j] = blam[j];
            }
            bb[nCon] = delz;

            // Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
            // AA = [Alam     a
            //       a'    -zet/z];
            for (int i = 0; i < nCon; i++)
            {
               for (int j = 0; j < nVar; j++)
               {
                  AA[i * nCon + j] = diaglamyi[i] * een[i] + GG[i * nCon + j] * diagxinv[j] * GG[i
                                                                                                 * nCon + j];
               }
               AA[nCon * nCon + i] = a[i];
               AA[3 * nCon + i] = a[i];
            }
            AA[4 * nCon] = -zet / z;
            // ----------------------------------------------------------------------------
            //solut = AA\bb --> solve linear system of equations using Gaussian elimination
            for (int i = 0; i < (nCon - 1); i++)
            {
               for (int j = i + 1; j < nCon; j++)
               {
                  double ratio = AA[j * nCon + i] / AA[i * nCon + i];
                  for (int k = i; k < nCon; k++)
                  {
                     AA[j * nCon + k] -= ratio * AA[i * nCon + k];
                  }
                  bb[j] -= ratio * bb[i];
               }
            }
            // back substitution
            for (int i = nCon - 1; i >= 0; i--)
            {
               double sum = 0.0;
               for (int j = i + 1; j < nCon; j++)
               {
                  sum += AA[i * nCon + j] * solut[j];
               }
               solut[i] = (bb[i] - sum) / AA[i * nCon + i];
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
               double sum = 0.0;
               for (int j = 0; j < nCon; j++)
               {
                  sum -= GG[j * nVar + i] * dlam[j];
               }
               dx[i] = -sum / diagx[i] - delx[i] / diagx[i];
            }
         }
         else
         {
            for (int i = 0; i < nCon; i++)
            {
               diaglamyiinv[i] = eem[i] / diaglamyi[i];
               dellamyi[i] = dellam[i] + dely[i] / diagy[i];
            }
            // Axx = spdiags(diagx,0,nVar,nVar) + GG'*spdiags(diaglamyiinv,0,nCon,nCon)*GG;
            // AA = [Axx      axz
            //       axz'     azz];
            for (int i = 0; i < nVar; i++)
            {
               for (int j = 0; j < nCon; j++)
               {
                  AA[i * nCon + j] = diagx[i] * een[i] + GG[j * nCon + i] * diaglamyiinv[j] * GG[j
                                                                                                 * nCon + i];
               }
            }
            // axz = -GG'*(a./diaglamyi)
            for (int i = 0; i < nVar; i++)
            {
               double sum = 0.0;
               for (int j = 0; j < nCon; j++)
               {
                  sum -= GG[j * nVar + i] * (a[j] / diaglamyi[j]);
               }
               axz[i] = sum;
               AA[nVar * nCon + i] = axz[i];
               AA[nVar * nCon + nVar + i] = axz[i];
            }
            for (int i = 0; i < nCon; i++)
            {
               azz[i] = zet / z + a[i] * (a[i] / diaglamyi[i]);
               AA[nVar * nCon + 2 * nVar + i] = azz[i];
            }
            // bb = [-bx'; -bz]
            // bx = delx - GG'*(dellamyi./diaglamyi)
            // bz = delz - a'*(dellamyi./diaglamyi)
            double sum1 = 0.0;
            for (int i = 0; i < nVar; i++)
            {
               double sum = 0.0;
               for (int j = 0; j < nCon; j++)
               {
                  sum -= GG[j * nVar + i] * (dellamyi[j] / diaglamyi[j]);
               }
               bb[i] = -(delx[i] + sum);
               sum1 += a[i] * (dellamyi[i] / diaglamyi[i]);
            }
            bb[nVar] = -(delz - sum1);
            // ----------------------------------------------------------------------------
            //solut = AA\bb --> solve linear system of equations using Gaussian elimination
            for (int i = 0; i < (nVar - 1); i++)
            {
               for (int j = i + 1; j < nVar; j++)
               {
                  double ratio = AA[j * nVar + i] / AA[i * nVar + i];
                  for (int k = i; k < nVar; k++)
                  {
                     AA[j * nVar + k] -= ratio * AA[i * nVar + k];
                  }
                  bb[j] -= ratio * bb[i];
               }
            }
            // back substitution
            for (int i = nVar - 1; i >= 0; i--)
            {
               double sum = 0.0;
               for (int j = i + 1; j < nVar; j++)
               {
                  sum += AA[i * nVar + j] * solut[j];
               }
               solut[i] = (bb[i] - sum) / AA[i * nVar + i];
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
               double sum = 0.0;
               for (int j = 0; j < nVar; j++)
               {
                  sum += GG[i * nVar + j] * dx[j];
               }
               dlam[i] = sum / diaglamyi[i] - dz * (a[i] / diaglamyi[i]) + dellamyi[i] /
                         diaglamyi[i];
            }
         }

         for (int i = 0; i < nCon; i++)
         {
            dy[i] = -dely[i] / diagy[i] + dlam[i] / diagy[i];
            dmu[i] = -mu[i] + epsvecm[i] / y[i] - (mu[i] * dy[i]) / y[i];
            ds[i] = -s[i] + epsvecm[i] / lam[i] - (s[i] * dlam[i]) / lam[i];
         }
         for (int i = 0; i < nVar; i++)
         {
            dxsi[i] = -xsi[i] + epsvecn[i] / (x[i] - alfa[i]) - (xsi[i] * dx[i]) /
                      (x[i] - alfa[i]);
            deta[i] = -eta[i] + epsvecn[i] / (beta[i] - x[i]) + (eta[i] * dx[i]) /
                      (beta[i] - x[i]);

         }
         dzet = -zet + epsi / z - zet * dz / z;

         // Fill arrays:
         // xx = [y z lam xsi eta mu zet s]
         // dxx = [dy dz dlam dxsi deta dmu dzet ds]
         for (int i = 0; i < nCon; i++)
         {
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
            xx[nCon + 1 + nCon + i] = xsi[i];
            xx[nCon + 1 + nCon + nVar + i] = eta[i];
            dxx[nCon + 1 + nCon + i] = dxsi[i];
            dxx[nCon + 1 + nCon + nVar + i] = deta[i];
         }
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
            stepalfa[i] = -1.01*dx[i] / (x[i] - alfa[i]);
            if (stepalfa[i] > stmalfa)
            {
               stmalfa = stepalfa[i];
            }
            stepbeta[i] = 1.01*dx[i] / (beta[i] - x[i]);
            if (stepbeta[i] > stmbeta)
            {
               stmbeta = stepbeta[i];
            }
         }
         stmalbe = std::max(stmalfa, stmbeta);
         stmalbexx = std::max(stmalbe, stmxx);
         stminv = std::max(stmalbexx, 1.0);
         steg = 1.0 / stminv;

         xold = x;
         yold = y;
         zold = z;
         lamold = lam;
         xsiold = xsi;
         etaold = eta;
         muold = mu;
         zetold = zet;
         sold = s;

         itto = 0;
         resinew = 2.0 * residunorm;

         while (resinew > residunorm && itto < 50)
         {
            itto++;
            for (int i = 0; i < nVar; ++i)
            {
               x[i] = xold[i] + steg * dx[i];
               xsi[i] = xsiold[i] + steg * dxsi[i];
               eta[i] = etaold[i] + steg * deta[i];

               ux1[i] = upp[i] - x[i];
               xl1[i] = x[i] - low[i];
               ux2[i] = ux1[i];
               xl2[i] = xl1[i];
               uxinv1[i] = een[i] / ux1[i];
               xlinv1[i] = een[i] / xl1[i];

            }
            for (int i = 0; i < nCon; ++i)
            {
               y[i] = yold[i] + steg * dy[i];
               lam[i] = lamold[i] + steg * dlam[i];
               mu[i] = muold[i] + steg * dmu[i];
               s[i] = sold[i] + steg * ds[i];
            }
            z = zold + steg * dz;
            zet = zetold + steg * dzet;

            // plam & qlam
            for (int i = 0; i < nCon; i++)
            {
               plam[i] = p0[i];
               qlam[i] = q0[i];
               for (int j = 0; j < nVar; j++)
               {
                  plam[i] += P[i * nVar + j] * lam[j];
                  qlam[i] += Q[i * nVar + j] * lam[j];
                  gvec[j] = P[i * nVar + j] * uxinv1[j] + Q[i * nVar + j] * xlinv1[j];
               }
            }

            // dpsidx, rex, rey, relam, rexsi, reeta
            rez = a0 - zet;
            for (int i = 0; i < nVar; i++)
            {
               dpsidx[i] = plam[i] / ux2[i] - qlam[i] / xl2[i];
               rex[i] = dpsidx[i] - xsi[i] + eta[i];
               rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
               rez -= a[i] * lam[i];
               relam[i] = gvec[i] - a[i] * z - y[i] + s[i] - b[i];
               rexsi[i] = xsi[i] * (x[i] - alfa[i]) - epsvecn[i];
               reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
            }

            for (int i = 0; i < nCon; i++)
            {
               remu[i] = mu[i] * y[i] - epsvecm[i];
               res[i] = lam[i] * s[i] - epsvecm[i];
            }
            rezet = zet * z - epsi;

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
            // New residual
            double sum = 0.0;
            for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
            {
               sum += residu[i] * residu[i];
               residumax = std::max(residumax, std::abs(residu[i]));
            }
            // Norm of the residual
            resinew = std::sqrt(sum);

            steg = steg / 2.0;
         }
         residunorm = resinew;
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
   // Update new values
   for (int i = 0; i < nVar; i++)
   {
      xmma[i] = x[i];
      xsimma[i] = xsi[i];
      etamma[i] = eta[i];
   }
   for (int i = 0; i < nCon; i++)
   {
      ymma[i] = y[i];
      lamma[i] = lam[i];
      mumma[i] = mu[i];
      smma[i] = s[i];
   }
   *zmma = z;
   *zetmma = zet;
}

void MMA::kktcheck(int nCon, int nVar, double* x, double* y, double z,
                   double* lam, double* xsi, double* eta,
                   double* mu, double zet, double* s,
                   double* xmin, double* xmax,
                   double* dfdx, double* gx, double* dgdx,
                   double a0, double* a, const double* c, double* d,
                   double* residu, double* residunorm, double* residumax)
{

   double* rex = new double(nVar);
   double* rey = new double(nCon);
   double rez, rezet;
   double* relam = new double(nCon);
   double* rexsi = new double(nVar);
   double* reeta = new double(nVar);
   double* remu = new double(nCon);
   double* res = new double(nCon);
   double* residu1 = new double(nVar + nCon + 1);
   double* residu2 = new double(2 * nVar + 3 * nCon + 2);


   double sum = 0.0;
   for (int j = 0; j < nVar; j++)
   {
      sum += dgdx[j] * lam[j];
   }
   for (int i = 0; i < nVar; i++)
   {
      rex[i] = dfdx[i] + sum - xsi[i] + eta[i];
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
   *residunorm = 0.0;
   *residumax = 0.0;
   for (int i = 0; i < (3 * nVar + 4 * nCon + 2); i++)
   {
      *residunorm += residu[i] * residu[i];
      *residumax = std::max(*residumax, std::abs(residu[i]));
   }
   // Norm of the residual
   *residunorm = std::sqrt(*residunorm);
}

/*void MMA::KKTresidual(double* xval, double* dfdx, double* gx, double* dgdx,
                      double* xmin, double* xmax, double* norm2, double* normInf)
{
   *norm2 = 1.0;
   *normInf = 1.0;
}
*/


void MMA::Restart(double* xo1, double* xo2, double* xo3)
{
   printf("Restart not implemented yet.\n");
}

void MMA::SetAsymptotes(double init, double decrease, double increase)
{

}

} // end mma namespace
