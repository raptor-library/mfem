#ifndef DIFFUSION_SOLVER_HPP
#define DIFFUSION_SOLVER_HPP

#include "mfem.hpp"

namespace mfem
{

class DiffusionCoeff:public mfem::Coefficient
{
public:
   DiffusionCoeff()
   {
   }

   virtual
   ~DiffusionCoeff()
   {

   }

   void SetDensity(mfem::Coefficient* coeff)
   {
      dcoeff=coeff;
   }

   virtual double Eval(mfem::ElementTransformation &T,
                       const mfem::IntegrationPoint &ip)
   {
      if (dcoeff == nullptr)
      {
         mfem_error("dcoeff is nullptr!");
      }

      double dens=dcoeff->Eval(T,ip);


      return 1.0* std::pow(dens, 3.0);
   }

private:

   mfem::Coefficient* dcoeff = nullptr;

};

class PDEFilter
{

public:
   PDEFilter(mfem::ParMesh* mesh_, int order_=2)
   {
      pmesh=mesh_;
      int dim=pmesh->Dimension();
      // TODO why not DG space
      //fec = new DG_FECollection(order_, dim, BasisType::GaussLobatto);
      fec = new H1_FECollection(order_,dim);
      fes = new ParFiniteElementSpace(pmesh,fec);

      sol.SetSize(fes->GetTrueVSize()); sol=0.0;
      rhs.SetSize(fes->GetTrueVSize()); rhs=0.0;
      adj.SetSize(fes->GetTrueVSize()); adj=0.0;

      solgf.SetSpace(fes);
      adjgf.SetSpace(fes);

      SetLinearSolver();
   }

   ~PDEFilter()
   {
      delete ls;
      delete prec;

      delete fes;
      delete fec;

      delete b;

      delete a;
   }

   /// Set the Linear Solver
   void SetLinearSolver(double rtol=1e-8, double atol=1e-12, int miter=2000)
   {
      linear_rtol=rtol;
      linear_atol=atol;
      linear_iter=miter;
   }

   void freeMemory()
   {
      delete b;
      b = nullptr;
   }

   /// Solves the forward problem.
   void FSolve();

   void ASolve(mfem::Vector& rhs);

   /// Returns the solution
   mfem::ParGridFunction& GetSolution() {return solgf;}

   void SetDesignGF(
      mfem::ParGridFunction *  designGF )
   {
      designGF_ = designGF;

      delete mDensCoeff;

      mDensCoeff = new mfem::GridFunctionCoefficient(designGF_);
   }

   ParFiniteElementSpace * GetFES()
   {
      return fes;
   }

   /// Returns the solution vector.
   mfem::Vector& GetSol() {return sol;}

   void GetSol(ParGridFunction& sgf)
   {
      sgf.SetSpace(fes); sgf.SetFromTrueDofs(sol);
   }

   /// Returns the adjoint solution vector.
   mfem::Vector& GetAdj() {return adj;}

   void GetAdj(ParGridFunction& agf)
   {
      agf.SetSpace(fes); agf.SetFromTrueDofs(adj);
   }


private:
   mfem::ParMesh* pmesh;

   ParBilinearForm *a = nullptr;
   ParLinearForm *b = nullptr;

   //solution true vector
   mfem::Vector sol;
   mfem::Vector adj;
   mfem::Vector rhs;
   mfem::ParGridFunction solgf;
   mfem::ParGridFunction adjgf;

   mfem::ParGridFunction *  designGF_ = nullptr;
   mfem::Coefficient*       mDensCoeff = nullptr;

   mfem::FiniteElementCollection *fec;
   mfem::ParFiniteElementSpace    *fes;

   //Linear solver parameters
   double linear_rtol;
   double linear_atol;
   int linear_iter;

   int print_level = 1;

   const double alpha = 1.0;

   double filterRad = 5e-4;

   //mfem::HypreBoomerAMG *prec = nullptr; //preconditioner
   mfem::HypreILU *prec = nullptr;
   //mfem::CGSolver *ls = nullptr;  //linear solver
   mfem::GMRESSolver *ls = nullptr;

   mfem::Array<int> ess_tdofv;
};

class Diffusion_Solver
{

public:
   Diffusion_Solver(mfem::ParMesh* mesh_, int order_=2)
   {
      pmesh=mesh_;
      int dim=pmesh->Dimension();
      // TODO why not DG space
      //fec = new DG_FECollection(order_, dim, BasisType::GaussLobatto);
      fec = new H1_FECollection(order_,dim);
      fes = new ParFiniteElementSpace(pmesh,fec);

      sol.SetSize(fes->GetTrueVSize()); sol=0.0;
      rhs.SetSize(fes->GetTrueVSize()); rhs=0.0;
      adj.SetSize(fes->GetTrueVSize()); adj=0.0;

      solgf.SetSpace(fes);
      adjgf.SetSpace(fes);

      SetLinearSolver();
   }

   ~Diffusion_Solver()
   {
      delete ls;
      delete prec;

      delete fes;
      delete fec;

      delete b;

      delete a;
   }

   /// Set the Linear Solver
   void SetLinearSolver(double rtol=1e-8, double atol=1e-12, int miter=2000)
   {
      linear_rtol=rtol;
      linear_atol=atol;
      linear_iter=miter;
   }

   /// Solves the forward problem.
   void FSolve();

   void ASolve(mfem::Vector& rhs);

   /// Adds Dirichlet BC
   void AddDirichletBC(int id, double val)
   {
      bc[id]=mfem::ConstantCoefficient(val);
      AddDirichletBC(id,bc[id]);
   }

   /// Adds Dirichlet BC
   void AddDirichletBC(int id, mfem::Coefficient& val)
   {
      bcc[id]=&val;
   }

   /// Adds Neumann BC
   void AddNeumannBC(int id, double val)
   {
      nc[id]=mfem::ConstantCoefficient(val);
      AddNeumannBC(id,nc[id]);
   }

   /// Adds Neumann BC
   void AddNeumannBC(int id, mfem::Coefficient& val)
   {
      ncc[id]=&val;
   }

   /// Returns the solution
   mfem::ParGridFunction& GetSolution() {return solgf;}

   void SetDesignGF(
      mfem::ParGridFunction *  designGF )
   {
      designGF_ = designGF;

      delete mDensCoeff;

      mDensCoeff = new mfem::GridFunctionCoefficient(designGF_);
   }

   ParFiniteElementSpace * GetFES()
   {
      return fes;
   }

   /// Returns the solution vector.
   mfem::Vector& GetSol() {return sol;}

   void GetSol(ParGridFunction& sgf)
   {
      sgf.SetSpace(fes); sgf.SetFromTrueDofs(sol);
   }

   /// Returns the adjoint solution vector.
   mfem::Vector& GetAdj() {return adj;}

   void GetAdj(ParGridFunction& agf)
   {
      agf.SetSpace(fes); agf.SetFromTrueDofs(adj);
   }

   void Postprocess();

private:
   mfem::ParMesh* pmesh;

   ParBilinearForm *a = nullptr;
   ParLinearForm *b = nullptr;

   //solution true vector
   mfem::Vector sol;
   mfem::Vector adj;
   mfem::Vector rhs;
   mfem::ParGridFunction solgf;
   mfem::ParGridFunction adjgf;

   mfem::ParGridFunction *  designGF_ = nullptr;
   mfem::Coefficient*       mDensCoeff = nullptr;

   mfem::FiniteElementCollection *fec;
   mfem::ParFiniteElementSpace    *fes;

   //Linear solver parameters
   double linear_rtol;
   double linear_atol;
   int linear_iter;

   int print_level = 1;

   const double alpha = 1.0;

   //mfem::HypreBoomerAMG *prec = nullptr; //preconditioner
   mfem::HypreILU *prec = nullptr;
   //mfem::CGSolver *ls = nullptr;  //linear solver
   mfem::GMRESSolver *ls = nullptr;

   // holds DBC in coefficient form
   std::map<int, mfem::Coefficient*> bcc;

   // holds internal DBC
   std::map<int, mfem::ConstantCoefficient> bc;

   // holds NBC in coefficient form
   std::map<int, mfem::Coefficient*> ncc;

   // holds internal NBC
   std::map<int, mfem::ConstantCoefficient> nc;

   mfem::Array<int> ess_tdofv;

   ParaViewDataCollection * mPvdc = nullptr;
};

class ThermalComplianceIntegrator : public NonlinearFormIntegrator
{
public:

   ThermalComplianceIntegrator()
   {};

   void SetFieldsAndMicrostructure(
      ParGridFunction * tempfield_,
      ParGridFunction * desfield_)
   {
      tempGF         =tempfield_;
      designGF       =desfield_;
   }

   virtual
   double GetElementEnergy(const FiniteElement &el, ElementTransformation &Tr,
                           const Vector &elfun);

   virtual
   void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr,
                              const Vector &elfun, Vector &elvect);

   virtual
   void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr,
                            const Vector &elfun, DenseMatrix &elmat)
   {
      {
         mfem::mfem_error("ThermalComplianceIntegrator::AssembleElementGrad is not defined!");
      }
   }

private:
   mfem::ParGridFunction             * tempGF          = nullptr;
   mfem::ParGridFunction             * designGF        = nullptr;

};

class ThermalComplianceIntegrator_1 : public LinearFormIntegrator
{
public:

   ThermalComplianceIntegrator_1()
   {};

   void SetDesignAndTempGF(
      ParGridFunction * desfield_,
      ParGridFunction * tempfield_)
   {
      tempGF         =tempfield_;
      designGF       =desfield_;
   }

   virtual
   double GetElementEnergy(const FiniteElement &el, ElementTransformation &Tr,
                           const Vector &elfun)
   {
      {
         mfem::mfem_error("ThermalComplianceIntegrator_1::AssembleElementGrad is not defined!");
      }
      return 0.0;
   }

   virtual
   void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr,
                               Vector &elvect);

   virtual
   void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr,
                            const Vector &elfun, DenseMatrix &elmat)
   {
      {
         mfem::mfem_error("ThermalComplianceIntegrator_1::AssembleElementGrad is not defined!");
      }
   }

private:
   mfem::ParGridFunction             * tempGF          = nullptr;
   mfem::ParGridFunction             * designGF        = nullptr;

};

class DiffusionAdjointPostIntegrator : public LinearFormIntegrator
{
public:

   DiffusionAdjointPostIntegrator()
   {};

   void SetAdjoint(ParGridFunction* Adjoint_)
   {
      AdjointGF=Adjoint_;
   };

   void SetDesignAndTempGF(
      mfem::ParGridFunction * desfield_,
      mfem::ParGridFunction * temp_ )
   {
      desfield               =desfield_;
      tempGF                = temp_;
   }

   virtual
   double GetElementEnergy(const FiniteElement &el, ElementTransformation &Tr,
                           const Vector &elfun)
   {
      {
         mfem::mfem_error("ADVDiffAdjointPostIntegrator::GetElementEnergy is not defined!");
      }
      return 0.0; //FIX THIS
   }

   virtual
   void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr,
                               Vector &elvect);

   virtual
   void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr,
                            const Vector &elfun, DenseMatrix &elmat)
   {
      {
         mfem::mfem_error("ADVDiffAdjointPostIntegrator::AssembleElementGrad is not defined!");
      }
   }

private:
   mfem::ParGridFunction             * AdjointGF       = nullptr;
   mfem::ParGridFunction             * desfield        = nullptr;
   mfem::ParGridFunction             * tempGF          = nullptr;

};

class VolIntegrator:public NonlinearFormIntegrator
{
public:

   VolIntegrator()
   {};

   void SetDesingField(ParGridFunction * desfield_)
   {
      desfield=desfield_;
   }

   virtual
   double GetElementEnergy(const FiniteElement &el, ElementTransformation &Tr,
                           const Vector &elfun);

   virtual
   void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr,
                              const Vector &elfun, Vector &elvect);

   virtual
   void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr,
                            const Vector &elfun, DenseMatrix &elmat);

private:
   mfem::ParGridFunction             * desfield        = nullptr;
};

class ThermalComplianceQoI
{
public:
   ThermalComplianceQoI()
   {};

   ~ThermalComplianceQoI()
   { delete nf;};

   void SetFESpaceAndField(
      ParFiniteElementSpace       * fes,
      ParGridFunction             * designGF_,
      ParGridFunction             * tempGF_)
   {
      dfes                   = fes;
      designGF               = designGF_;
      tempGF                 = tempGF_;
   };

   double Eval();

   void Grad(Vector& grad);

private:

   ParFiniteElementSpace       * dfes            = nullptr;
   ParNonlinearForm            * nf              = nullptr;
   ThermalComplianceIntegrator * intgr           = nullptr;
   ParGridFunction             * designGF        = nullptr;
   ParGridFunction             * tempGF          = nullptr;
};

class VolumeQoI
{
public:
   VolumeQoI()
   {};

   ~VolumeQoI()
   { delete nf;};

   void SetDesignFES(ParFiniteElementSpace* fes)
   {
      dfes=fes;
   }

   void SetDesField(ParGridFunction* desfield_)
   {
      desfield=desfield_;
   }

   double Eval();

   void Grad(Vector& grad);

private:

   ParFiniteElementSpace       * dfes            = nullptr;
   ParNonlinearForm            * nf              = nullptr;
   VolIntegrator               * intgr           = nullptr;
   ParGridFunction             * desfield        = nullptr;
};




}

#endif
