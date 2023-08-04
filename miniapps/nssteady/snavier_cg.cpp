#include "snavier_cg.hpp"

namespace mfem{

/// Constructor
SNavierPicardCGSolver::SNavierPicardCGSolver(ParMesh* mesh,
                                             int vorder,
                                             int porder,
                                             double kin_vis_,
                                             bool verbose):
pmesh(mesh), vorder(vorder), porder(porder), verbose(verbose)
{
   // mesh
   dim=pmesh->Dimension();

   // FE collection and spaces for velocity and pressure
   vfec=new H1_FECollection(vorder,dim);
   vfes=new ParFiniteElementSpace(pmesh,vfec,dim);
   pfec=new H1_FECollection(porder);
   pfes=new ParFiniteElementSpace(pmesh,pfec,1);

   // determine spaces dimension
   vdim = vfes->GetTrueVSize();
   pdim = pfes->GetTrueVSize(); 
   
   // initialize vectors of essential attributes
   vel_ess_attr.SetSize(pmesh->bdr_attributes.Max());      vel_ess_attr=0;
   vel_ess_attr_x.SetSize(pmesh->bdr_attributes.Max());    vel_ess_attr_x=0;
   vel_ess_attr_y.SetSize(pmesh->bdr_attributes.Max());    vel_ess_attr_y=0;
   vel_ess_attr_z.SetSize(pmesh->bdr_attributes.Max());    vel_ess_attr_z=0;

   // initialize GridFunctions
   v_gf.SetSpace(vfes);  v_gf=0.0;
   vk_gf.SetSpace(vfes); vk_gf=0.0;
   z_gf.SetSpace(vfes);  z_gf=0.0;
   p_gf.SetSpace(pfes);  p_gf=0.0;
   pk_gf.SetSpace(pfes); pk_gf=0.0;

   // initialize vectors
   v    = new Vector(vdim); *v = 0.0;
   vk   = new Vector(vdim); *vk = 0.0; 
   z    = new Vector(vdim); *z = 0.0; 
   p    = new Vector(pdim); *p = 0.0; 
   pk   = new Vector(pdim); *pk = 0.0; 
   fv   = new Vector(vdim); *fv = 0.0; 
   fp   = new Vector(pdim); *fp = 0.0; 
   rhs1 = new Vector(vdim); *rhs1 = 0.0; 
   rhs2 = new Vector(pdim); *rhs2 = 0.0;
   rhs3 = new Vector(vdim); *rhs3 = 0.0; 
   tmp  = new Vector(vdim); *tmp = 0.0; 
   
   // setup GridFunctionCoefficients
   vk_vc = new VectorGridFunctionCoefficient(&vk_gf);
   pk_c  = new GridFunctionCoefficient(&pk_gf);

   // set kinematic viscosity
   kin_vis.constant = kin_vis_;

   // Error computation setup
   err_v = err_p = 0;
   norm_v = norm_p = 0;
   int order_quad = std::max(2, 2*vorder+1);
   for (int i=0; i < Geometry::NumGeom; ++i)
   {
      irs[i] = &(IntRules.Get(i, order_quad));
   }

}



/// Public Interface

void SNavierPicardCGSolver::AddVelDirichletBC(VectorCoefficient *coeff, Array<int> &attr)
{
   vel_dbcs.emplace_back(attr, coeff);

   // Check for duplicate
   for (int i = 0; i < attr.Size(); ++i)
   {
      MFEM_ASSERT(( (vel_ess_attr[i] || vel_ess_attr_x[i] || vel_ess_attr_y[i] || vel_ess_attr_z[i]) && attr[i]) == 0,
                  "Duplicate boundary definition detected.");
      if (attr[i] == 1)
      {
         vel_ess_attr[i] = 1;
      }
   }

   // Output
   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Velocity Dirichlet BC (full) to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }
}

void SNavierPicardCGSolver::AddVelDirichletBC(Coefficient *coeff, Array<int> &attr, int &dir)
{
   // Add bc container to list of componentwise velocity bcs
   vel_dbcs_xyz.emplace_back(attr, coeff, dir);

   // Check for duplicate and add attributes for current bc to global list (for that specific component)
   for (int i = 0; i < attr.Size(); ++i)
   {
      switch (dir) {
            case 0: // x 
               dir_string = "x";
               MFEM_ASSERT(( (vel_ess_attr[i] || vel_ess_attr_x[i]) && attr[i]) == 0,
                           "Duplicate boundary definition for x component detected.");
               if (attr[i] == 1){vel_ess_attr_x[i] = 1;}
               break;
            case 1: // y
               dir_string = "y";
               MFEM_ASSERT(( (vel_ess_attr[i] || vel_ess_attr_y[i]) && attr[i]) == 0,
                           "Duplicate boundary definition for y component detected.");
               if (attr[i] == 1){vel_ess_attr_y[i] = 1;}
               break;
            case 2: // z
               dir_string = "z";
               MFEM_ASSERT(( (vel_ess_attr[i] || vel_ess_attr_z[i]) && attr[i]) == 0,
                           "Duplicate boundary definition for z component detected.");
               if (attr[i] == 1){vel_ess_attr_z[i] = 1;}
               break;
            default:;
         }      
   }

   // Output
   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Velocity Dirichlet BC ( " << dir_string << " component) to attributes: " << std::endl;
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << ", ";
         }
      }
      mfem::out << std::endl;
   }
}

void SNavierPicardCGSolver::AddVelDirichletBC(VectorCoefficient *coeff, int &attr)
{
   // Create array for attributes and mark given mark given mesh boundary
   ess_attr_tmp = 0;
   ess_attr_tmp[ attr - 1] = 1;

   // Call AddVelDirichletBC accepting array of essential attributes
   AddVelDirichletBC(coeff, ess_attr_tmp);
}

void SNavierPicardCGSolver::AddVelDirichletBC(Coefficient *coeff, int &attr, int &dir)
{
   // Create array for attributes and mark given mark given mesh boundary
   ess_attr_tmp = 0;
   ess_attr_tmp[ attr - 1] = 1;

   // Call AddVelDirichletBC accepting array of essential attributes
   AddVelDirichletBC(coeff, ess_attr_tmp, dir);
}

void SNavierPicardCGSolver::AddTractionBC(VectorCoefficient *coeff, Array<int> &attr)
{
   traction_bcs.emplace_back(attr, coeff);

   for (int i = 0; i < attr.Size(); ++i)
   {
      MFEM_ASSERT(( (vel_ess_attr[i] || vel_ess_attr_x[i] || vel_ess_attr_y[i] || vel_ess_attr_z[i]) && attr[i]) == 0,
                  "Trying to enforce traction bc on dirichlet boundary.");
   }

   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Traction (Neumann) BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }
}

void SNavierPicardCGSolver::AddAccelTerm(VectorCoefficient *coeff, Array<int> &attr)
{
   accel_terms.emplace_back(attr, coeff);

   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Acceleration term to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }
}

void SNavierPicardCGSolver::SetFixedPointSolver(SolverParams params)
{
   sParams = params;    
}

void SNavierPicardCGSolver::SetLinearSolvers( SolverParams params1,
                                            SolverParams params2,
                                            SolverParams params3)
{
   s1Params = params1;
   s2Params = params2;
   s3Params = params3;                          
}

void SNavierPicardCGSolver::SetAlpha(double &alpha_, const AlphaType &type_)
{
   alpha0    = alpha_;
   alphaType = type_;
}

void SNavierPicardCGSolver::SetInitialConditionVel(VectorCoefficient &v_in)
{
   // Project coefficient onto velocity ParGridFunction
   v_gf.ProjectCoefficient(v_in);

   // Initialize provisional velocity and velocity at previous iteration
   v_gf.GetTrueDofs(*v);
   *z = *v;
   z_gf.SetFromTrueDofs(*z);
   //*vk = *v;                         // CHECK: do we need to initialize also vk?
   //vk_gf.SetFromTrueDofs(*vk);
}

void SNavierPicardCGSolver::SetInitialConditionPres(Coefficient &p_in)
{
   // Project coefficient onto pressure ParGridFunction
   p_gf.ProjectCoefficient(p_in);

   // Initialize pressure at previous iteration
   p_gf.GetTrueDofs(*p);
   //*pk = *p;                     // CHECK: do we need to initialize also pk?
   //pk_gf.SetFromTrueDofs(*pk);
}

void SNavierPicardCGSolver::Setup()
{
   /// 1. Setup and assemble bilinear forms 
   K_form = new ParBilinearForm(vfes);
   C_form = new ParBilinearForm(vfes);
   B_form = new ParMixedBilinearForm(vfes, pfes);

   K_form->AddDomainIntegrator(new VectorDiffusionIntegrator(kin_vis));
   B_form->AddDomainIntegrator(new VectorDivergenceIntegrator);

   K_form->Assemble();  K_form->Finalize();
   B_form->Assemble();  B_form->Finalize();
  
   K = K_form->ParallelAssemble();
   B = B_form->ParallelAssemble();


   /// 2. Setup and assemble linear form for rhs
   f_form = new ParLinearForm(vfes);
   // Adding forcing terms
   for (auto &accel_term : accel_terms)
   {
      f_form->AddDomainIntegrator( new VectorDomainLFIntegrator(*accel_term.coeff) );
   }
   // Adding traction bcs
   for (auto &traction_bc : traction_bcs)
   {
      f_form->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(*traction_bc.coeff), traction_bc.attr);
   }
   f_form->Assemble(); 
   fv = f_form->ParallelAssemble(); 
   

   /// 3. Apply boundary conditions
   // Extract to list of true dofs
   vfes->GetEssentialTrueDofs(vel_ess_attr_x,vel_ess_tdof_x,0);
   vfes->GetEssentialTrueDofs(vel_ess_attr_y,vel_ess_tdof_y,1);
   vfes->GetEssentialTrueDofs(vel_ess_attr_z,vel_ess_tdof_z,2);
   vfes->GetEssentialTrueDofs(vel_ess_attr, vel_ess_tdof_full);
   vel_ess_tdof.Append(vel_ess_tdof_x);
   vel_ess_tdof.Append(vel_ess_tdof_y);
   vel_ess_tdof.Append(vel_ess_tdof_z);
   vel_ess_tdof.Append(vel_ess_tdof_full);

   // Projection of coeffs (full velocity applied)
   for (auto &vel_dbc : vel_dbcs)
   {
      v_gf.ProjectBdrCoefficient(*vel_dbc.coeff, vel_dbc.attr);
   }
   v_gf.GetTrueDofs(*v);

   // Projection of coeffs (velocity component applied)
   ParGridFunction tmp_gf(vfes);        // temporary velocity gf for projection
   Vector          tmp_vec(vdim);       // temporary velocity vector for projection
   Array<int>      tmp_tdofs;
   for (auto &vel_dbc : vel_dbcs_xyz)
   {
      VectorArrayCoefficient tmp_coeff(dim);                           // Set coefficient with right component
      tmp_coeff.Set(vel_dbc.dir, vel_dbc.coeff, false);
      tmp_gf.ProjectBdrCoefficient(tmp_coeff, vel_dbc.attr);           // Project on dummy gf
      tmp_gf.GetTrueDofs(tmp_vec);

      vfes->GetEssentialTrueDofs(vel_dbc.attr,tmp_tdofs,vel_dbc.dir);  // Update solution dofs
      for(int i=0;i<tmp_tdofs.Size();i++)
      {
         v[tmp_tdofs[i]]=tmp_vec[tmp_tdofs[i]];
      }      
   }

   // Initialize solution gf with vector containing projected coefficients
   v_gf.SetFromTrueDofs(*v);

   /// 4. Apply transformation for essential bcs
   // NOTE: alternatively the following function performs the modification on both matrix and rhs
   //       EliminateRowsCols(const Array<int> &rows_cols, const HypreParVector &x, HypreParVector &b)
   Ke = K->EliminateRowsCols(vel_ess_tdof);  // Remove rows/cols for ess tdofs
   Be  = B->EliminateCols(vel_ess_tdof);
   Bt  = B->Transpose(); 
   Bte = Be->Transpose(); 

   //ModifyRHS(vel_ess_tdof, Ke, *v, *fv);         // Modify rhs for velocity    rhs -= Ke u
   //ModifyRHS(vel_ess_tdof, Bte, *p, *fv);         // Modify rhs for velocity    rhs -= Bt p   Not needed as we don't prescribe EssPres bcs

   ModifyRHS(vel_ess_tdof, Be, *v, *fp, false);  // Modify rhs for pressure ( no pressure ess bcs)   fp = 0 - Be v

   // Update grid function and vector for provisional velocity
   // CHECK: do we need to initialize also vk?
   *z  = *v;
   z_gf.SetFromTrueDofs(*z);

   /// 5. Setup solvers and preconditioners
   //  5.1 Velocity prediction       A = K + alpha*C(uk)
   // solved with CGSolver preconditioned with HypreBoomerAMG (elasticity version)
   invA_pc = new HypreBoomerAMG();
   invA_pc->SetSystemsOptions(dim);
   invA_pc->SetPrintLevel(s1Params.pl);
   invA_pc->SetElasticityOptions(vfes);
   invA = new CGSolver(pmesh->GetComm());
   invA->iterative_mode = true;     // uses 2nd argument of mult as initial guess
   invA->SetPrintLevel(s1Params.pl);
   invA->SetRelTol(s1Params.rtol);
   //invA->SetAbsTol(s1Params.atol);
   invA->SetMaxIter(s1Params.maxIter);

   // 5.2 Pressure correction       S = B K^{-1} Bt
   // solved with CGSolver preconditioned with HypreBoomerAMG 
   // NOTE: test different approaches to deal with Schur Complement:
   // * now using Jacobi, but this may not be a good approximation when involving Brinkman Volume Penalization
   // * alternative may be to use Multigrid to get better approximation
   HypreParVector* Kd = new HypreParVector(pmesh->GetComm(), K->GetGlobalNumRows(), K->GetRowStarts());
   HypreParMatrix* local = new HypreParMatrix(*Bt); // local = Bt
   K->GetDiag(*Kd);
   local->InvScaleRows(*Kd);  // local = Kd^{-1} Bt
   S = ParMult(B, local);     // S = B Kd^{-1} Bt
   delete local; local=nullptr;
   delete Kd; Kd = nullptr;

   invS_pc = new HypreBoomerAMG(*S);
   //invS_pc->SetSystemsOptions(dim);
   invS_pc->SetPrintLevel(s2Params.pl);
   invS = new CGSolver(pmesh->GetComm());
   invS->iterative_mode = true;
   invS->SetOperator(*S);
   //invS->SetPreconditioner(*invS_pc);
   invS->SetPrintLevel(s2Params.pl);
   invS->SetRelTol(s2Params.rtol);
   //invS->SetAbsTol(s2Params.atol);
   invS->SetMaxIter(s2Params.maxIter);

   // 5.3 Velocity correction
   // solved with CGSolver preconditioned with HypreBoomerAMG 
   invK_pc = new HypreBoomerAMG(*K);
   //invK_pc->SetSystemsOptions(dim);
   invK_pc->SetPrintLevel(s3Params.pl);
   invK = new CGSolver(pmesh->GetComm());
   invK->iterative_mode = true;
   invK->SetOperator(*K);
   invK->SetPreconditioner(*invK_pc);
   invK->SetPrintLevel(s3Params.pl);
   invK->SetRelTol(s3Params.rtol);
   //invK->SetAbsTol(s3Params.atol);
   invK->SetMaxIter(s3Params.maxIter);

}

void SNavierPicardCGSolver::FSolve()
{
#ifdef MFEM_DEBUG
   PrintMatricesVectors( "setup", 0); // Export matrices/vectors before step
#endif

   PrintInfo();

   if (pmesh->GetMyRank() == 0)
   {
      mfem::out << std::endl;
      mfem::out << "=========================================================="<< std::endl;
      mfem::out << "======    Picard-aCT Steady Navier-Stokes Solver    ======"<< std::endl;
      mfem::out << "=========================================================="<< std::endl;
   }

   timer.Clear();
   timer.Start();

   // Print header
   mfem::out << std::endl;
   mfem::out << std::setw(7) << "" << std::setw(3) << "It" << std::setw(8)
             << "Res" << std::setw(12) << "AbsTol" << "\n";


   // Export initial solution
   SaveResults( 0 );

   for (iter = 0; iter < sParams.maxIter; iter++)
   {
      // Update parameter alpha
      UpdateAlpha();

      // Solve current iteration.
      Step();

      // Output results
      SaveResults( iter+1 );

      // Compute errors.
      ComputeError();

      // Update solution at previous iterate and gridfunction coefficients.
      UpdateSolution();

      // Print results
      mfem::out << iter << "   " << std::setw(3)
                << std::setprecision(2) << std::scientific << err_v
                << "   " << sParams.atol << "\n";

      // Check convergence.
      if (err_v < sParams.atol)
      {
         out << "Solver converged to steady state solution \n";
         flag = 1;
         break;
      }

   }

   timer.Stop();
}


void SNavierPicardCGSolver::SetupOutput( const char* folderPath, bool visit_, bool paraview_,
                                         DataCollection::Format par_format )
{
   visit    = visit_;
   paraview = paraview_;
   outfolder = folderPath;

   // Creating output directory if not existent
   if (mkdir(folderPath, 0777) == -1) {std::cerr << "Error :  " << strerror(errno) << std::endl;}

   // Create output data collections   
   if( visit )
   {
      visit_dc = new VisItDataCollection("Results-VISit", pmesh);
      visit_dc->SetPrefixPath(folderPath);
      visit_dc->RegisterField("velocity", &v_gf);
      visit_dc->RegisterField("intermediate velocity", &z_gf);
      visit_dc->RegisterField("pressure", &p_gf);
      visit_dc->SetFormat(par_format);
   }

   if( paraview )
   {
      paraview_dc = new ParaViewDataCollection("Results-Paraview", pmesh);
      paraview_dc->SetPrefixPath(folderPath);
      paraview_dc->SetLevelsOfDetail(vorder);
      paraview_dc->SetDataFormat(VTKFormat::BINARY);
      paraview_dc->SetHighOrderOutput(true);
      paraview_dc->SetCycle(0);
      paraview_dc->SetTime(0.0);
      paraview_dc->RegisterField("velocity",&v_gf);
      paraview_dc->RegisterField("intermediate velocity",&z_gf);
      paraview_dc->RegisterField("pressure",&p_gf);
   }   

}

/// Private Interface

void SNavierPicardCGSolver::Step()
{
   /// Assemble convective term with new velocity vk and modify matrix for essential bcs.
   delete C_form; C_form = nullptr;
   C_form = new ParBilinearForm(vfes);
   C_form->AddDomainIntegrator(new VectorConvectionIntegrator(*vk_vc, 1));
   C_form->Assemble(); C_form->Finalize();
   C = C_form->ParallelAssemble();              

   /// Solve.
   // 1: Velocity prediction      ( K + alpha*C(vk) ) z = f - (1-alpha)*C(uk)*k
   // Assemble rhs
   rhs1->Set(1,*fv);                           // rhs1 = fv
   C->Mult(*vk,*tmp);                          // tmp = C*vk
   rhs1->Add(alpha-1, *tmp);                   // rhs1 = fv += (alpha-1)*C*vk      
   
   // Assemble operator
   A = Add(1.0, *K, alpha, *C);                // A = Km + alpha C                   
   A->Add(1, *Ke);                             // A = Km + Ke + alpha C 

#ifdef MFEM_DEBUG
   PrintMatricesVectors( "prestep", iter);  // Export matrices/vectors after assembly of A, before modifications
#endif
  
   Ae = A->EliminateRowsCols(vel_ess_tdof);

   ModifyRHS(vel_ess_tdof, Ae, *z, *rhs1);

   invA->SetOperator(*A);     
   invA_pc->SetOperator(*A);

   invA->Mult(*rhs1,*z);

#ifdef MFEM_DEBUG
   PrintMatricesVectors( "step1", iter);  // Export matrices/vectors after 1st step
#endif

   // 2: Pressure correction                   B*K^-1*B^T p = B*z
   FullMult(B, Be, *z, *rhs2);                 // rhs2 = B z
   rhs2->Add(1, *fp);                          // rhs2 += fp        

   invS->Mult(*rhs2, *p);

#ifdef MFEM_DEBUG
   PrintMatricesVectors( "step2", iter);  // Export matrices/vectors after 2nd step
#endif

   // 3: Velocity correction         K u = K*z - B^T*p = f - (1-alpha)*C(uk)*uk - alpha C(uk) z - B^T*p  
   // NOTE: Could be more efficient storing and reusing rhs1, SparseMatrix -alpha C(uk)
   FullMult(K, Ke, *z, *rhs3);   // rhs3 = K z
   FullMult(Bt, Bte, *p, *tmp);  // tmp  = Bt p  
   rhs3->Add(-1, *tmp);          // rhs3 -= tmp        
   ModifyRHS(vel_ess_tdof, Ke, *v, *rhs3);

   invK->Mult(*rhs3,*v);

#ifdef MFEM_DEBUG
   PrintMatricesVectors( "step3", iter);  // Export matrices/vectors after 3rd step
#endif

   /// Update GridFunctions for solution.
   v_gf.SetFromTrueDofs(*v);
   p_gf.SetFromTrueDofs(*p);
   z_gf.SetFromTrueDofs(*z);

   delete A; A = nullptr;
   delete Ae; Ae = nullptr;
   delete C; C = nullptr;
}

void SNavierPicardCGSolver::ComputeError()
{
   err_v  = v_gf.ComputeL2Error(*vk_vc);
   norm_v = ComputeGlobalLpNorm(2., *vk_vc, *pmesh, irs);
   err_p  = p_gf.ComputeL2Error(*pk_c);
   norm_p = ComputeGlobalLpNorm(2., *pk_c, *pmesh, irs);

   if (verbose)
   {
      out << "|| v - v_k || / || v_k || = " << err_v / norm_v << "\n";
      out << "|| p - p_k || / || p_k || = " << err_p / norm_p << "\n";
   }
}

void SNavierPicardCGSolver::UpdateSolution()
{
   *vk = *v;
   *z  = *v;
   z_gf.SetFromTrueDofs(*z);
   vk_gf.SetFromTrueDofs(*vk);

   *pk = *p;
   pk_gf.SetFromTrueDofs(*pk);
}


void SNavierPicardCGSolver::UpdateAlpha()
{
   if ( alphaType == AlphaType::CONSTANT) { alpha = alpha0;}
   else {  MFEM_ABORT("Error: SNavierPicardCGSolver::UpdateAlpha does not implement"
                       "adaptive update of the segregation parameter yet!");} // NYI!
}


void SNavierPicardCGSolver::ModifyRHS(Array<int> &ess_tdof_list, HypreParMatrix* mat_e,
                                      Vector &sol, Vector &rhs, bool copy_sol)
{
   // Initialize temporary vector for solution
   Vector tmp(sol);
   tmp.SetSubVectorComplement(ess_tdof_list, 0.0);

   // Perform elimination
   mat_e->Mult(-1.0, tmp, 1.0, rhs); // rhs -= mat_e*sol
 
   // Set rhs equal to solution at essential tdofs
   if( copy_sol )
   {
      int idx;
      for (int i = 0; i < ess_tdof_list.Size(); i++)
      {
         idx = ess_tdof_list[i];
         rhs(idx) = sol(idx);
      }
   }

}

void SNavierPicardCGSolver::FullMult(HypreParMatrix* mat, HypreParMatrix* mat_e, Vector &x, Vector &y)
{
   mat->Mult(x, y);        // y =  mat x
   mat_e->AddMult(x, y);   // y += mat_e x
}

void SNavierPicardCGSolver::SaveResults( int iter )
{

   if ( visit ) // Save to GLVis format if visit enabled
   { 
      visit_dc->SetCycle(iter);
      visit_dc->SetTime(iter);
      visit_dc->Save();
   }
   if ( paraview ) // Save to Paraview format if visit enabled
   {
      paraview_dc->SetCycle(iter);
      paraview_dc->SetTime(iter);
      paraview_dc->Save();
   }
}

void SNavierPicardCGSolver::PrintInfo()
{
   int fes_sizeVel = vfes->GlobalVSize();
   int fes_sizePres = pfes->GlobalVSize();

   if (pmesh->GetMyRank() == 0)
   {
      mfem::out << std::endl;
      mfem::out << "NAVIER version: " << SNAVIER_CG_VERSION << std::endl
                << "MFEM version: " << MFEM_VERSION << std::endl
                << "MFEM GIT: " << MFEM_GIT_STRING << std::endl
                << "Velocity #DOFs: " << fes_sizeVel << std::endl
                << "Pressure #DOFs: " << fes_sizePres << std::endl;
   }
}

void SNavierPicardCGSolver::PrintMatricesVectors( const char* id, int num )
{
   // Create folder
   std::string folderName(outfolder);
   folderName += "/MatVecs_iter";
   folderName += std::to_string(num);

   if (mkdir(folderName.c_str(), 0777) == -1) {std::cerr << "Error :  " << strerror(errno) << std::endl;}

   //Create files
   std::ofstream K_file(std::string(folderName) + '/' + "K_" + std::string(id) + ".dat");
   std::ofstream C_file(std::string(folderName) + '/' + "C_" + std::string(id) + ".dat");
   std::ofstream A_file(std::string(folderName) + '/' + "A_" + std::string(id) + ".dat");
   std::ofstream S_file(std::string(folderName) + '/' + "S_" + std::string(id) + ".dat");
   std::ofstream B_file(std::string(folderName) + '/' + "B_" + std::string(id) + ".dat");
   std::ofstream Bt_file(std::string(folderName) + '/' + "Bt_" + std::string(id) + ".dat");

   std::ofstream fv_file(std::string(folderName) + '/' + "fv_" + std::string(id) + ".dat");
   std::ofstream fp_file(std::string(folderName) + '/' + "fp_" + std::string(id) + ".dat");
   std::ofstream rhs1_file(std::string(folderName) + '/' + "rhs1_" + std::string(id) + ".dat");
   std::ofstream rhs2_file(std::string(folderName) + '/' + "rhs2_" + std::string(id) + ".dat");
   std::ofstream rhs3_file(std::string(folderName) + '/' + "rhs3_" + std::string(id) + ".dat");
   std::ofstream dofs_file(std::string(folderName) + '/' + "dofs_" + std::string(id) + ".dat");

   std::ofstream v_file(std::string(folderName) + '/' + "v_" + std::string(id) + ".dat");
   std::ofstream vk_file(std::string(folderName) + '/' + "vk_" + std::string(id) + ".dat");
   std::ofstream p_file(std::string(folderName) + '/' + "p_" + std::string(id) + ".dat");
   std::ofstream pk_file(std::string(folderName) + '/' + "pk_" + std::string(id) + ".dat");
   std::ofstream z_file(std::string(folderName) + '/' + "z_" + std::string(id) + ".dat");

   // Print matrices in matlab format
   K->PrintMatlab(K_file);
   if(C==nullptr)
   {
      C = new HypreParMatrix();
      C->PrintMatlab(C_file);
      delete C; C = nullptr;
   }
   else
   {
      C->PrintMatlab(C_file);
   }

   if(A==nullptr)
   {
      A = new HypreParMatrix();
      A->PrintMatlab(C_file);
      delete A; A = nullptr;
   }
   else
   {
      A->PrintMatlab(A_file);
   }

   S->PrintMatlab(S_file);
   B->PrintMatlab(B_file);
   Bt->PrintMatlab(Bt_file);

   fv->Print(fv_file,1);
   fp->Print(fp_file,1);
   rhs1->Print(rhs1_file,1);
   rhs2->Print(rhs2_file,1);
   rhs3->Print(rhs3_file,1);

   v->Print(v_file,1);
   vk->Print(vk_file,1);
   p->Print(p_file,1);
   pk->Print(pk_file,1);
   z->Print(z_file,1);


   for (int i = 0; i < vel_ess_tdof.Size(); ++i)
   {
      dofs_file << vel_ess_tdof[i] << std::endl;
   }
   dofs_file.close();

}

/// Destructor
SNavierPicardCGSolver::~SNavierPicardCGSolver()
{
   delete vfes;
   delete vfec;
   delete pfes;
   delete pfec;

   delete K_form;
   delete B_form;
   delete C_form; 
   delete f_form;

   delete K;
   delete B;
   delete Bt;
   //delete A;
   //delete C;
   delete S;
   delete Ke;
   delete Be;
   delete Bte;
   //delete Ae;

   delete fv;
   delete fp;
   delete rhs1;
   delete rhs2;
   delete rhs3;
   delete tmp;

   delete vk_vc;
   delete pk_c;

   delete fcoeff;
   delete traction;

   delete invA;     
   delete invK;     
   delete invS;     

   delete invA_pc;  
   delete invK_pc;  
   delete invS_pc;  

   delete paraview_dc;
   delete visit_dc;
}

}