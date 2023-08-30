#ifndef __MTOP_MMA__HPP
#define __MTOP_MMA__HPP

#include "mfem.hpp"
#include <iostream>

namespace mma
{

class MMA
{
private:

   enum SubProblemType
   {
      CLASSIC,
      MPIClassic,
      ENDENUM
   };

   // Private: MMA-specific
   double raa0, move, albefa, asyinit, asyincr, asydecr, xmamieps, lowmin, lowmax,
          uppmin, uppmax, zz;
   double *factor, *xmami, *pq0, *p, *q, *pq, *PQ;
   double *ux1, *xl1;
   // Global: Old design variables
   double *xo1, *xo2;
   //Convergence study
   double kktnorm;
   bool isInitialized = false;

   void setGlobals(int nVar, int nCon);
   void setMMA(int nVar, int nCon);
   void freeGlobals();
   void freeMMA();
   void SubProblemFactory(int nVar, int nCon, enum SubProblemType type = CLASSIC);

public:

   // Local vectors
   double *a, *b, *c, *d;
   double a0, machineEpsilon, epsimin;
   double z, zet;
   int nCon, nVar;

   // Global: Asymptotes, bounds, objective approx., constraint approx.
   double *xmin, *xmax, *low, *upp, *p0, *q0, *P, *Q;
   double *x, *y, *alfa, *beta, *xsi, *eta, *lam, *mu, *s;


   // Construct using subproblem
   MMA(int nVar, int nCon, double *xval, double *xxmin, double *xxmax,
       enum SubProblemType type = CLASSIC );

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
   }

   // SUBPROBLEM BASE CLASS
   class SubProblemBase
   {
   protected:
      MMA* mma_ptr;

   public:
      SubProblemBase(MMA* mma);
      virtual ~SubProblemBase()
      {
         this->mma_ptr = nullptr;
      };
      virtual void Perform() = 0;
   };

   // SubProblem according to Svanberg, implemented in serial
   class SubProblemClassic : public SubProblemBase
   {
   private:
      int ittt, itto, itera;
      double epsi, rez, rezet, delz, dz, dzet, azz, stmxx, stmalfa, stmbeta, stmalbe,
             sum, stmalbexx, stminv, steg, zold, zetold,
             residunorm, residumax, resinew;
      double *sum1, *ux1, *xl1, *plam, *qlam, *gvec, *residu, *GG, *Puxinv, *Qxlinv,
             *delx, *dely, *dellam, *dellamyi, *diagx, *diagy, *diaglam, *diaglamyi, *bb,
             *bb1, *Alam, *AA, *AA1, *dlam, *dx, *dy, *dxsi, *deta, *dmu, *Axx, *axz, *ds,
             *xx, *dxx, *stepxx, *stepalfa, *stepbeta, *xold, *yold, *lamold, *xsiold,
             *etaold, *muold, *sold;

      void setSubProb(int nVar, int nCon);
      void freeSubProb();
      double kktcheck(double* y, double* const dfdx, double* const gx,
                      double* const dgdx, double* x);

   public:

      SubProblemClassic(MMA* mma, int nVar, int nCon);
      ~SubProblemClassic()
      {
         this->freeSubProb();
      };
      void Perform()  override;
   };

   // SubProblem according to Svanberg, implemented in parallel
   class SubProblemClassicMPI : public SubProblemBase
   {
   private:
      void setSubProb(int nVar, int nCon);
      void freeSubProb();

   public:
      SubProblemClassicMPI(MMA* mma, int nVar, int nCon);
      ~SubProblemClassicMPI()
      {
         this->freeSubProb();
      };
      void Perform() override;
   };

   //instance of class
   SubProblemBase *mSubProblem;


   void initializeGlobals(int nVariables, int nConstraints, double *xval,
                          double *xxmin, double *xxmax)
   {
      nVar = nVariables;
      nCon = nConstraints;
      xmin = xxmin;
      xmax = xxmax;

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

   void Update(int iter, double* const fval, double* const dfdx, double* const gx,
               double* const dgdx, double* xval);

   double* getLow();
   double* getUpp();
   double getKKT();

   // Options
   // Return necessary data for restart
   void setRestart(double* xval, int iter, std::string name,
                   std::string path = "./");
   void outputRestart(double* xval, int iter, std::string name,
                      std::string path = "./");

};

}

#endif