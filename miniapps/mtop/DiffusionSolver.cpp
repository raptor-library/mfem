#include "DiffusionSolver.hpp"

namespace mfem
{

void PDEFilter::FSolve()
{

   //    Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero.
   ParGridFunction x(fes);
   x = 0.0;

   sol = solgf;

   ess_tdofv.DeleteAll();

   mfem::ConstantCoefficient filterRadius(filterRad);

   //    Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system.
   if (b==nullptr)
   {
      b = new ParLinearForm(fes);

      mfem::Coefficient * LoadCoeff = new GridFunctionCoefficient( designGF_ );
      b->AddDomainIntegrator(new DomainLFIntegrator(*LoadCoeff));
      b->Assemble();
   }

   //    Set up the bilinear form a(.,.) on the finite element space
   if (a==nullptr)
   {
      a = new ParBilinearForm(fes);


      // add diffusion integrators
      a->AddDomainIntegrator(new DiffusionIntegrator(filterRadius));
      a->AddDomainIntegrator(new MassIntegrator());
      a->Assemble();
      a->Finalize();
   }

   HypreParMatrix A;
   HypreParVector B, X;

   //allocate the preconditioner and the linear solver
   if (prec==nullptr)
   {
      //prec = new HypreBoomerAMG();
      prec = new HypreILU ();
      prec->SetLevelOfFill(5);
      prec->SetPrintLevel(print_level);
   }

   if (ls==nullptr)
   {
      //ls = new CGSolver(pmesh->GetComm());
      ls = new GMRESSolver(pmesh->GetComm());
      ls->SetAbsTol(linear_atol);
      ls->SetRelTol(linear_rtol);
      ls->SetMaxIter(linear_iter);
      ls->SetPrintLevel(print_level);
      ls->SetKDim(1000);
      ls->SetPreconditioner(*prec);
   }

   a->FormLinearSystem(ess_tdofv, sol, *b, A, X, B);

   ls->SetOperator(A);
   ls->Mult(B, X);


   sol = X;     // copy solution
   solgf.SetFromTrueDofs(sol);

   delete (a);
   a = nullptr;
}

void PDEFilter::ASolve(mfem::Vector& rhs)
{
   MFEM_ASSERT( ls != nullptr, "Liner solver does not exist");

   adj = 0.0;
   ess_tdofv.DeleteAll();

   mfem::ConstantCoefficient filterRadius(filterRad);

   if (a==nullptr)
   {
      a = new ParBilinearForm(fes);

      a->AddDomainIntegrator(new DiffusionIntegrator(filterRadius));
      a->AddDomainIntegrator(new MassIntegrator());
      a->Assemble();
      a->Finalize();

   }

   HypreParMatrix A;
   HypreParVector B, X;

   //std::cout << "\n\n b= " << b->Size() << "\n\n"<<std::endl;
   a->FormLinearSystem(ess_tdofv, sol, *b, A, X, B);

   mfem::HypreParMatrix* tTransOp = (&A)->Transpose();

   ls->SetOperator(*tTransOp);
   ls->Mult(rhs, adj);


   delete tTransOp;
   adjgf.SetFromTrueDofs(adj);

   delete (a);
   //delete b;

   a = nullptr;
   //b = nullptr;
}

void Diffusion_Solver::FSolve()
{

   //    Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero.
   ParGridFunction x(fes);
   x = 0.0;

   sol = solgf;

   ess_tdofv.DeleteAll();

   // set the boundary conditions
   {
      for (auto it=bcc.begin(); it!=bcc.end(); it++)
      {
         mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
         ess_bdr=0;
         ess_bdr[it->first -1]=1;
         mfem::Array<int> ess_tdof_list;
         fes->GetEssentialTrueDofs(ess_bdr,ess_tdof_list);
         ess_tdofv.Append(ess_tdof_list);
         solgf.ProjectBdrCoefficient(*(it->second),ess_bdr);
      }

      //copy BC values from the grid function to the solution vector
      {
         solgf.GetTrueDofs(rhs);
         for (int ii=0; ii<ess_tdofv.Size(); ii++)
         {
            sol[ess_tdofv[ii]]=rhs[ess_tdofv[ii]];
         }
      }
   }

   //the BC are setup in the solution vector sol
   std::cout<<"BC dofs size="<<ess_tdofv.Size()<<std::endl;

   mfem::DiffusionCoeff * DiffCoeff = new DiffusionCoeff();
   DiffCoeff->SetDensity(mDensCoeff);

   //    Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system.
   if (b==nullptr)
   {
      b = new ParLinearForm(fes);

      mfem::Coefficient * LoadCoeff = new ConstantCoefficient( 1.0 );
      b->AddDomainIntegrator(new DomainLFIntegrator(*LoadCoeff));
      b->Assemble();
   }

   //    Set up the bilinear form a(.,.) on the finite element space
   if (a==nullptr)
   {
      a = new ParBilinearForm(fes);


      // add diffusion integrators
      a->AddDomainIntegrator(new DiffusionIntegrator(*DiffCoeff));
      a->Assemble();
      a->Finalize();
   }

   HypreParMatrix A;
   HypreParVector B, X;

   //allocate the preconditioner and the linear solver
   if (prec==nullptr)
   {
      //prec = new HypreBoomerAMG();
      prec = new HypreILU ();
      prec->SetLevelOfFill(5);
      prec->SetPrintLevel(print_level);
   }

   if (ls==nullptr)
   {
      //ls = new CGSolver(pmesh->GetComm());
      ls = new GMRESSolver(pmesh->GetComm());
      ls->SetAbsTol(linear_atol);
      ls->SetRelTol(linear_rtol);
      ls->SetMaxIter(linear_iter);
      ls->SetPrintLevel(print_level);
      ls->SetKDim(1000);
      ls->SetPreconditioner(*prec);
   }

   a->FormLinearSystem(ess_tdofv, sol, *b, A, X, B);

   ls->SetOperator(A);
   ls->Mult(B, X);


   sol = X;     // copy solution
   solgf.SetFromTrueDofs(sol);

   delete DiffCoeff;
   delete (a);
   a = nullptr;
}

void Diffusion_Solver::ASolve(mfem::Vector& rhs)
{
   MFEM_ASSERT( ls != nullptr, "Liner solver does not exist");

   adj = 0.0;
   ess_tdofv.DeleteAll();
   // set the boundary conditions
   {
      for (auto it=bcc.begin(); it!=bcc.end(); it++)
      {
         mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
         ess_bdr=0;
         ess_bdr[it->first -1]=1;
         mfem::Array<int> ess_tdof_list;
         fes->GetEssentialTrueDofs(ess_bdr,ess_tdof_list);
         ess_tdofv.Append(ess_tdof_list);
         adjgf.ProjectBdrCoefficient(*(it->second),ess_bdr);

      }

      //copy BC values from the grid function to the solution vector
      {
         //adjgf.GetTrueDofs(rhs);
         for (int ii=0; ii<ess_tdofv.Size(); ii++)
         {
            adj[ess_tdofv[ii]]=0.0;
            rhs[ess_tdofv[ii]]=0.0;
         }
      }
   }

   mfem::DiffusionCoeff * DiffCoeff = new DiffusionCoeff();
   DiffCoeff->SetDensity(mDensCoeff);

   if (a==nullptr)
   {
      a = new ParBilinearForm(fes);

      a->AddDomainIntegrator(new DiffusionIntegrator(*DiffCoeff));
      a->Assemble();
      a->Finalize();

   }

   HypreParMatrix A;
   HypreParVector B, X;

   a->FormLinearSystem(ess_tdofv, sol, *b, A, X, B);

   mfem::HypreParMatrix* tTransOp = (&A)->Transpose();

   ls->SetOperator(*tTransOp);
   ls->Mult(rhs, adj);


   delete tTransOp;
   adjgf.SetFromTrueDofs(adj);

   delete DiffCoeff;
   delete (a);
   delete b;

   a = nullptr;
   b = nullptr;
}

void Diffusion_Solver::Postprocess()
{
   if ( true )
   {
      mPvdc = new ParaViewDataCollection("AdvDiff", pmesh);
      mPvdc->SetDataFormat(VTKFormat::BINARY32);
      mPvdc->SetHighOrderOutput(true);
      //mPvdc->SetLevelsOfDetail(mCtk.order);
      mPvdc->SetCycle(0);
      mPvdc->SetTime(0.0);
      mPvdc->RegisterField("Temperature", &solgf);
      mPvdc->RegisterField("density", designGF_);
      mPvdc->Save();

   }
}

double ThermalComplianceIntegrator::GetElementEnergy(
   const FiniteElement &el,
   ElementTransformation &Tr,
   const Vector &elfun)
{
   //integrate the dot product gradTT^T \kappa GradT

   const int dim=el.GetDim();
   {
      const int spaceDim=Tr.GetDimension();
      if (dim!=spaceDim)
      {
         mfem::mfem_error("ThermalComplianceIntegrator::GetElementEnergy is not define on manifold meshes.");
      }
   }

   Vector gradT; gradT.SetSize(dim);
   Vector tempvec;  tempvec.SetSize(dim);

   const IntegrationRule *ir = nullptr;
   int order= 2 * el.GetOrder() + Tr.OrderGrad(&el)+designGF->FESpace()->GetOrder(
                 Tr.ElementNo);
   ir=&IntRules.Get(Tr.GetGeometryType(),order);

   double w;
   double energy=0.0;
   for (int i=0; i<ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      w=Tr.Weight();
      w = ip.weight * w;

      double DesingVal = designGF->GetValue(Tr,ip);
      tempGF     ->GetGradient(Tr,gradT);

      tempvec = 0.0;
      tempvec.Add(std::pow(DesingVal, 3.0),gradT);

      // Mult gradTT^T \kappa GradT
      energy=energy+w*0.5*(gradT*tempvec);
   }

   return energy;

}

//the finite element space is the space of the filtered design
void ThermalComplianceIntegrator::AssembleElementVector(
   const FiniteElement &el,
   ElementTransformation &Tr,
   const Vector &elfun,
   Vector &elvect)
{
   const int dof=el.GetDof();
   const int dim=el.GetDim();
   {
      const int spaceDim=Tr.GetDimension();
      if (dim!=spaceDim)
      {
         mfem::mfem_error("ComplianceNLIntegrator::AssembleElementVector is not define on manifold meshes.");
      }
   }

   elvect.SetSize(dof); elvect=0.0;

   Vector shapef(dof);

   const IntegrationRule *ir = nullptr;
   int order= 2 * el.GetOrder() + Tr.OrderGrad(&el)+2*
              (designGF->FESpace()->GetOrder(Tr.ElementNo));
   ir=&IntRules.Get(Tr.GetGeometryType(),order);

   Vector gradT; gradT.SetSize(dim);
   Vector tempvec;  tempvec.SetSize(dim);

   double w;
   for (int i=0; i<ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      w=Tr.Weight();
      w = ip.weight * w;

      double DesingVal = designGF->GetValue(Tr,ip);
      tempGF     ->GetGradient(Tr,gradT);

      tempvec = 0.0;
      tempvec.Add(3.0*std::pow(DesingVal, 2.0),gradT);

      // Mult gradTT^T \kappa GradT
      double cpl =(gradT*tempvec);

      el.CalcShape(ip,shapef);
      elvect.Add(0.5*cpl*w,shapef);               // shapef from chain rule
   }
}

void ThermalComplianceIntegrator_1::AssembleRHSElementVect(
   const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
{
   const int dof=el.GetDof();
   const int dim=el.GetDim();
   {
      const int spaceDim=Tr.GetDimension();
      if (dim!=spaceDim)
      {
         mfem::mfem_error("ThermalComplianceIntegrator_1::AssembleElementVector is not define on manifold meshes.");
      }
   }

   elvect.SetSize(dof); elvect=0.0;

   Vector shapef(dof);
   DenseMatrix dshape_iso(dof, dim);
   DenseMatrix dshape_xyz(dof, dim);

   const IntegrationRule *ir = nullptr;
   int order= 2 * el.GetOrder() + Tr.OrderGrad(&el)+2*
              (designGF->FESpace()->GetOrder(Tr.ElementNo));
   ir=&IntRules.Get(Tr.GetGeometryType(),order);

   Vector gradT; gradT.SetSize(dim);
   Vector tempvec;  tempvec.SetSize(dim);

   double w;
   for (int i=0; i<ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      w=Tr.Weight();
      w = ip.weight * w;

      el.CalcDShape(ip, dshape_iso);
      el.CalcShape(ip, shapef);
      Mult(dshape_iso, Tr.InverseJacobian(), dshape_xyz);

      double DesingVal = designGF->GetValue(Tr,ip);
      tempGF     ->GetGradient(Tr,gradT);

      tempvec = 0.0;
      tempvec.Add(std::pow(DesingVal, 3.0),gradT);

      Vector dQdT(dof);
      dshape_xyz.Mult(tempvec, dQdT);
      elvect.Add(0.5*2.0*w,
                 dQdT);                              //2.0 because of gradT^T /kappa gradT
   }
}

void DiffusionAdjointPostIntegrator::AssembleRHSElementVect(
   const FiniteElement &el,
   ElementTransformation &Tr,
   Vector &elvect)
{
   const int dof=el.GetDof();
   const int dim=el.GetDim();
   {
      const int spaceDim=Tr.GetDimension();
      if (dim!=spaceDim)
      {
         mfem::mfem_error("ComplianceNLIntegrator::AssembleElementVector is not define on manifold meshes.");
      }
   }

   elvect.SetSize(dof); elvect=0.0;

   Vector shapef(dof);

   const IntegrationRule *ir = nullptr;
   int order= 2 * el.GetOrder() + Tr.OrderGrad(&el)+2*
              (desfield->FESpace()->GetOrder(Tr.ElementNo));
   ir=&IntRules.Get(Tr.GetGeometryType(),order);

   Vector gradT; gradT.SetSize(dim);
   Vector tempvec;  tempvec.SetSize(dim);

   DenseMatrix dshape_iso(dof, dim);
   DenseMatrix dshape_xyz(dof, dim);

   double w;
   for (int i=0; i<ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      w=Tr.Weight();
      w = ip.weight * w;

      el.CalcDShape(ip, dshape_iso);
      el.CalcShape(ip, shapef);
      Mult(dshape_iso, Tr.InverseJacobian(), dshape_xyz);

      tempGF     ->GetGradient(Tr,gradT);
      double DesingVal = desfield->GetValue(Tr,ip);

      tempvec = 0.0;
      tempvec.Add(3.0*std::pow(DesingVal, 2.0),gradT);

      Vector AdjointGrad(dim);
      AdjointGF->GetGradient(Tr,AdjointGrad);

      double ScalarVal = AdjointGrad*tempvec;

      elvect.Add(w*ScalarVal,shapef);
   }
}

double VolIntegrator::GetElementEnergy(
   const FiniteElement &el,
   ElementTransformation &Tr,
   const Vector &elfun)
{
   //const int dof=el.GetDof();
   const int dim=el.GetDim();
   {
      const int spaceDim=Tr.GetDimension();
      if (dim!=spaceDim)
      {
         mfem::mfem_error("VolIntegrator::GetElementEnergy is not define on manifold meshes.");
      }
   }

   const IntegrationRule *ir = nullptr;
   int order= 2 * el.GetOrder() + Tr.OrderGrad(&el)+desfield->FESpace()->GetOrder(
                 Tr.ElementNo);
   ir=&IntRules.Get(Tr.GetGeometryType(),order);

   double w;
   double energy=0.0;
   for (int i=0; i<ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      w=Tr.Weight();
      w = ip.weight * w;

      double Desingval = desfield->GetValue(Tr,ip);

      energy=energy+w*(Desingval);
   }

   return energy;

}

//the finite element space is the space of the filtered design
void VolIntegrator::AssembleElementVector(
   const FiniteElement &el,
   ElementTransformation &Tr,
   const Vector &elfun,
   Vector &elvect)
{
   const int dof=el.GetDof();
   const int dim=el.GetDim();
   {
      const int spaceDim=Tr.GetDimension();
      if (dim!=spaceDim)
      {
         mfem::mfem_error("ComplianceNLIntegrator::AssembleElementVector is not define on manifold meshes.");
      }
   }

   elvect.SetSize(dof); elvect=0.0;
   Vector shapef(dof);

   const IntegrationRule *ir = nullptr;
   int order= 2 * el.GetOrder() + Tr.OrderGrad(&el)+2*
              (desfield->FESpace()->GetOrder(Tr.ElementNo));
   ir=&IntRules.Get(Tr.GetGeometryType(),order);

   double w;
   for (int i=0; i<ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      w=Tr.Weight();
      w = ip.weight * w;

      //double DesingVal = desfield->GetValue(Tr,ip);

      el.CalcShape(ip,shapef);
      elvect.Add(w,shapef);
   }
}

void VolIntegrator::AssembleElementGrad(
   const FiniteElement &el,
   ElementTransformation &Tr,
   const Vector &elfun,
   DenseMatrix &elmat)
{
   {
      mfem::mfem_error("ComplianceNLIntegrator::AssembleElementGrad is not defined!");
   }
}

double ThermalComplianceQoI::Eval()
{
   if (designGF==nullptr)
   {
      MFEM_ABORT("fsolv of dfes in ThermalComplianceQoI should be set before calling the Eval method!");
   }


   if (nf==nullptr)
   {
      nf=new ParNonlinearForm(dfes);
      intgr=new ThermalComplianceIntegrator();
      nf->AddDomainIntegrator(intgr);
   }

   intgr->SetFieldsAndMicrostructure(
      tempGF,
      designGF);

   double rt=nf->GetEnergy(*designGF);

   return rt;
}

void ThermalComplianceQoI::Grad(Vector& grad)
{

   if (designGF==nullptr)
   {
      MFEM_ABORT("fsolv or dfes in ThermalComplianceQoI should be set before calling the Grad method!");
   }

   if (nf==nullptr)
   {
      nf=new ParNonlinearForm(dfes);
      intgr=new ThermalComplianceIntegrator();
      nf->AddDomainIntegrator(intgr);
   }
   intgr->SetFieldsAndMicrostructure(
      tempGF,
      designGF);

   nf->Mult(*designGF,grad);
}

double VolumeQoI::Eval()
{
   if (desfield==nullptr)
   {
      MFEM_ABORT("fsolv of dfes in VolumeQoI should be set before calling the Eval method!");
   }


   if (nf==nullptr)
   {
      nf=new ParNonlinearForm(dfes);
      intgr=new VolIntegrator();
      nf->AddDomainIntegrator(intgr);
   }

   intgr->SetDesingField(desfield);

   double rt=nf->GetEnergy(*desfield);

   return rt;
}

void VolumeQoI::Grad(Vector& grad)
{
   if (desfield==nullptr)
   {
      MFEM_ABORT("fsolv or dfes in VolumeQoI should be set before calling the Grad method!");
   }

   if (nf==nullptr)
   {
      nf=new ParNonlinearForm(dfes);
      intgr=new VolIntegrator();
      nf->AddDomainIntegrator(intgr);
   }

   intgr->SetDesingField(desfield);

   nf->Mult(*desfield,grad);
}

}
