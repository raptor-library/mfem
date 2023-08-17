#ifndef __MTOP_MMA__HPP
#define __MTOP_MMA__HPP

#include "mfem.hpp"
#include <iostream>

namespace mma
{

class MMA
{
private:
   // Local vectors
   double *a, *c, *d;
   double a0, zet, zetmma;
   int z;

   // Lagrange multipliers:
   double *lam, *xsi, *eta, *mu, *s;

   // Global: Asymptotes, bounds, objective approx., constraint approx.
   double *low, *upp, *alfa, *beta, *p0, *q0, *P, *Q, *xmma, *ymma, *zmma, *lamma, *xsimma, *etamma, *mumma, *smma;

   // Global: MMA-specific
   double epsimin, raa0, move, albefa, asyinit, asyincr, asydecr, xmamieps, lowmin, lowmax, uppmin, uppmax, zz;
   double *factor, *xmami, *pq0, *p, *q, *pq, *b, *PQ;

   // Global: Subproblem
   int ittt, itto, itera;
   double epsi, machineEpsilon, rez, rezet, delz, dz, dzet, azz, stmxx, stmalfa, stmbeta, stmalbe, sum, stmalbexx, stminv, steg, zold, zetold, residunorm, residumax, resinew;
   double *sum1, *epsvecn, *epsvecm, *x, *y, *ux1, *ux2, *ux3, *xl1, *xl2, *xl3, *uxinv1, *xlinv1, *plam, *qlam, *gvec, *dpsidx, *rex, *rey, *relam,
          *rexsi, *reeta, *remu, *res, *residu1, *residu2, *residu, *GG, *Puxinv, *Qxlinv, *delx, *dely, *dellam, *dellamyi, *diagx, *diagxinv,
          *diagy, *diagyinv, *diaglam, *diaglamyi, *diaglamyiinv, *blam, *bb, *Alam, *AA, *solut, *dlam, *dx, *dy, *dxsi, *deta, *dmu, *Axx,
          *axz, *ds, *xx, *dxx, *stepxx, *stepalfa, *stepbeta, *xold, *yold, *lamold, *xsiold, *etaold, *muold, *sold;

   // Global: Old design variables
   double *xo1, *xo2;

   //Convergence studies
   double kktnorm;

   bool isInitialized = false;

   void setGlobals(int nVar, int nCon)
   {
      low = new double[nVar];
      upp = new double[nVar];
      xo1 = new double[nVar];
      xo2 = new double[nVar];
      c = new double[nCon];
      d = new double[nCon];
      a = new double[nCon];
      xmma = new double[nVar];
      ymma = new double[nCon];
      zmma = new double[nCon];
      lam = new double[nCon];
      xsi = new double[nVar];
      eta = new double[nVar];
      mu = new double[nCon];
      s = new double[nCon];
      lamma = new double[nCon];
      xsimma = new double[nVar];
      etamma = new double[nVar];
      mumma = new double[nCon];
      smma = new double[nCon];
      q0 = new double[nVar];
      p0 = new double[nVar];
      P = new double[nCon * nVar];
      Q = new double[nCon * nVar];
      alfa = new double[nVar];
      beta = new double[nVar];
      z = zet = 1.0;
      kktnorm = 10;
      printf("Initialized ");

      isInitialized = true;
   }

   void setMMA(int nVar, int nCon)
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
   }
   void setSubProb(int nVar, int nCon)
   {
         //Allocate all subproblem variables
         epsi = 1.0;
         machineEpsilon = 1e-10;
         epsvecn = new double[nVar];
         epsvecm = new double[nCon];
         x = new double[nVar];
         y = new double[nCon];
         ittt = itto = itera = 0;
         ux1 = new double[nVar];
         ux2 = new double[nVar];
         ux3 = new double[nVar];
         xl1 = new double[nVar];
         xl2 = new double[nVar];
         xl3 = new double[nVar];
         uxinv1 = new double[nVar];
         xlinv1 = new double[nVar];
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
         diagxinv = new double[nVar];
         diagy = new double[nCon];
         diagyinv = new double[nCon];
         diaglam = new double[nCon];
         diaglamyi = new double[nCon];
         diaglamyiinv = new double[nCon];
         blam = new double[nCon];
         bb = new double[nVar + 1];
         Alam = new double[nCon * nCon];
         AA = new double[(nVar + 1) * (nVar + 1)];
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

public:
   // Construct using defaults subproblem penalization
   MMA(int nVar, int nCon, int iter)
   {
      printf("0 ");
      this->setGlobals(nVar, nCon);
      printf("1 ");
      this->setMMA(nVar, nCon);
      printf("2 ");
      this->setSubProb(nVar, nCon);
      printf("3 ");
      this->initializeGlobals(nVar, nCon, iter);
      printf("4\n");
   }

   //MMA(Comm, int nVar, int nCon, int iter)
   //{
   //   this->setGlobals(nVar, nCon);
   //   this->setMMA(nVar, nCon);
   //   this->setSubProb(nVar, nCon);
   //   this->initializeGlobals(nVar, nCon, iter);
   //}

   ~MMA()
   {
       this->freeGlobals();
       this->freeMMA();
       this->freeSubProb();
   }

   void freeGlobals()
   {
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
      delete[] lamma;
      delete[] xsimma;
      delete[] etamma;
      delete[] mumma;
      delete[] smma;
      delete[] q0;
      delete[] p0;
      delete[] P;
      delete[] Q;

      isInitialized = false;
   }

   void freeMMA()
   {
      delete[] factor;
      delete[] xmami;
      delete[] pq0;
      delete[] p;
      delete[] q;
      delete[] pq;
      delete[] b;
      delete[] PQ;
   }
   void freeSubProb()
   {
      delete[] sum1;
      delete[] epsvecn;
      delete[] epsvecm;
      delete[] x;
      delete[] y;
      delete[] ux1;
      delete[] ux2;
      delete[] ux3;
      delete[] xl1;
      delete[] xl2;
      delete[] xl3;
      delete[] uxinv1;
      delete[] xlinv1;
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
      delete[] diagxinv;
      delete[] diagy;
      delete[] diagyinv;
      delete[] diaglam;
      delete[] diaglamyi;
      delete[] diaglamyiinv;
      delete[] blam;
      delete[] bb;
      delete[] Alam;
      delete[] AA;
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

   void initializeGlobals(int nVar, int nCon, int iter)
   {
      if (iter == 0)
      {
         for (int i = 0; i < nVar; i++)
         {
            //MMAmain.xval[i] = 0.0;
            xo1[i] = 0.0;
            xo2[i] = 0.0;
            //xmin[i] = -2.0;
            //xmax[i] = 2.0;
            low[i] = -2.0;
            upp[i] = 2.0;
         }
      }
      else
      {
         /*
         std::ifstream input("Restart.dat");
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
         */
      }

      for (int i = 0; i < nCon; i++)
      {
         a[i] = 0.0;
         c[i] = 1000.0;
         d[i] = 1.0;
      }
      double a0 = 1.0;
      if (!isInitialized)
      {
         mfem::mfem_error("MMA::initialize() member data not initialized, call setGlobals first");
      }
      
      

   }

   void Update(int nVar, int nCon, int iter, double* xval, double* xmin,
                 double* xmax, double* fval, double* dfdx, double* gx, double* dgdx);

   // Set and solve a subproblem: return new xval
   void mmasub(int nVar, int nCon, int iter, double* xval, double* xmin, double* xmax, double* fval, double* dfdx, double* gx, double* dgdx);

   void kktcheck(int nCon, int nVar, double* x, double* y, 
                 double* xmin, double* xmax,
                 double* dfdx, double* gx, double* dgdx);

   void subsolv(int nVar, int nCon, double epsimin, double* b);


   // Options
   // Return necessary data for possible restart
   void Restart(double* xval, int length, int iter);

};

}

#endif