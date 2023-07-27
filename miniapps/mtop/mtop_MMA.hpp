#ifndef __MTOP_MMA__HPP
#define __MTOP_MMA__HPP

#include "mfem.hpp"

namespace mma
{

class MMA
{
public:
   // Construct using defaults subproblem penalization
   MMA(int n, int m, double * x, int sx);

   //MMA(MPI_Comm Comm,int n, int m, double * x, int sx);

   //~MMA()
   //{};


   // Set and solve a subproblem: return new xval
   void mmasub(int nVar, int nCon, int iter, double* xval, double* xmin,
               double* xmax, double* xold1, double* xold2, double* fval,
               double* dfdx, double* gx, double* dgdx, double* low,
               double* upp, double a0, double* a, double* c, double* d,
               double* xmma, double* ymma, double* zmma, double* lam,
               double* xsi, double* eta, double* mu, double& zet, double* s);


   // Return KKT residual norms (norm2 and normInf)
   //void KKTresidual(double* xval, double* dfdx, double* gx, double* dgdx,
   //                 double* xmin, double* xmax, double* norm2,
   //                 double* normInf);

   void kktcheck(int nCon, int nVar, double* x, double* y, double z,
                 double* lam, double* xsi, double* eta,
                 double* mu, double zet, double* s,
                 double* xmin, double* xmax,
                 double* dfdx, double* gx, double* dgdx,
                 double a0, double* a, const double* c, double* d,
                 double* kktnorm);

   void subsolv(int nVar, int nCon, double epsimin, double* low, double* upp,
                double* alfa, double* beta, double* p0,
                double* q0, double* P, double* Q,
                double a0, double* a, double* b, double* c,
                double* d, double* xmma, double* ymma,
                double* zmma, double* lamma, double* xsimma,
                double* etamma, double* mumma, double* zetmma, double* smma);


   // Options
   // Return necessary data for possible restart
   void Restart(double* xo1, double* xo2, double* xo3);

   // Set the aggresivity of the moving asymptotes
   void SetAsymptotes(double init, double decrease, double increase);


private:


   // Local vectors: elastic variables
   double* y;
   int  z;

   // Local vectors: Lagrange multipliers:
   double *lam, *mu, *s;

   // Global: Asymptotes, bounds, objective approx., constraint approx.
   double* L, U, alpha, beta, p0, q0, pij, qij;

   // Local: subproblem constant terms, dual gradient, dual hessian
   double *b, *grad, *Hess;

   // Global: Old design variables
   double* xo1, xo2;

};

}

#endif