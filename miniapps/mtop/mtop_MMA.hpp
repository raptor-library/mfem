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
   double *a = nullptr, *c = nullptr, *d;
   double a0, zet, zetmma, z, zmma;
   int nCon, nVar;

   // Lagrange multipliers:
   double *lam, *xsi, *eta, *mu, *s;

   // Global: Asymptotes, bounds, objective approx., constraint approx.
   double *xmin, *xmax, *low, *upp, *alfa, *beta, *p0, *q0, *P, *Q, *xmma, *ymma, *lamma, *xsimma, *etamma, *mumma, *smma;

   // Global: MMA-specific
   double epsimin, raa0, move, albefa, asyinit, asyincr, asydecr, xmamieps, lowmin, lowmax, uppmin, uppmax, zz;
   double *factor, *xmami, *pq0, *p, *q, *pq, *b, *PQ;

   // Global: Subproblem
   int ittt, itto, itera;
   double epsi, machineEpsilon, rez, rezet, delz, dz, dzet, azz, stmxx, stmalfa, stmbeta, stmalbe, sum, stmalbexx, stminv, steg, zold, zetold, residunorm, residumax, resinew;
   double *sum1, *epsvecn, *epsvecm, *x, *y, *ux1, *xl1, *plam, *qlam, *gvec, *dpsidx, *rex, *rey, *relam, *rexsi, *reeta, *remu, *res, *residu1, *residu2, *residu, *GG, 
          *Puxinv, *Qxlinv, *delx, *dely, *dellam, *dellamyi, *diagx, *diagy, *diaglam, *diaglamyi, *blam, *bb, *bb1, *Alam, *AA, *AA1, *solut, *dlam, *dx, *dy, *dxsi, *deta, 
          *dmu, *Axx, *axz, *ds, *xx, *dxx, *stepxx, *stepalfa, *stepbeta, *xold, *yold, *lamold, *xsiold, *etaold, *muold, *sold;

   // Global: Old design variables
   double *xo1, *xo2;

   //Convergence study
   double kktnorm;

   bool isInitialized = false;

   void setGlobals(int nVar, int nCon)
   {
      xmin = new double[nVar];
      xmax = new double[nVar];
      low = new double[nVar];
      upp = new double[nVar];
      xo1 = new double[nVar];
      xo2 = new double[nVar];
      c = new double[nCon];
      d = new double[nCon];
      a = new double[nCon];
      xmma = new double[nVar];
      ymma = new double[nCon];
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

   void sanityCheck()
   {
      //isFinite

      //mfem_error(a == nullptr, "a is nullprt");
   }

public:
   // Construct using defaults subproblem penalization
   MMA(int nVar, int nCon, int iter, double *xval, double *xxmin, double *xxmax)
   {
      this->setGlobals(nVar, nCon);
      this->setMMA(nVar, nCon);
      this->setSubProb(nVar, nCon);
      this->initializeGlobals(nVar, nCon, iter, xval, xxmin, xxmax);
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

   void initializeGlobals(int nVariables, int nConstraints, int iter, double *xval, double *xxmin, double *xxmax)
   {
      nVar = nVariables;
      nCon = nConstraints;
      xmin = xxmin;
      xmax = xxmax;

      if (iter == 0)
      {
         for (int i = 0; i < nVar; i++)
         {
            xo1[i] = 0.0;
            xo2[i] = 0.0;
            low[i] = xmin[i];
            upp[i] = xmax[i];
         }
      }
      else
      {
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
      }

      for (int i = 0; i < nCon; i++)
      {
         a[i] = 0.0;
         c[i] = 1000.0;
         d[i] = 1.0;
      }
      a0 = 1.0;
      
      if (!isInitialized)
      {
         mfem::mfem_error("MMA::initialize() member data not initialized, call setGlobals first");
      }
      
      

   }

   void Update(int iter, double* fval, double* dfdx, double* gx, double* dgdx, double* xval);

   // Set and solve a subproblem: return new xval
   void mmasub(int iter, double* fval, double* dfdx, double* gx, double* dgdx, double* xval);

   void kktcheck(double* y, double* dfdx, double* gx, double* dgdx, double* x);

   void subsolv();

   double* getLow();
   double* getUpp();
   double getKKT();

   // Options
   // Return necessary data for possible restart
   void Restart(double* xval, int iter);

};

}

#endif