//                       MFEM Example 0 - Parallel Version
//
// Compile with: make ex0p
//
// Sample runs:  mpirun -np 4 ex0p
//               mpirun -np 4 ex0p -m ../data/fichera.mesh
//               mpirun -np 4 ex0p -m ../data/square-disc.mesh -o 2
//
// Description: This example code demonstrates the most basic parallel usage of
//              MFEM to define a simple finite element discretization of the
//              Laplace problem -Delta u = 1 with zero Dirichlet boundary
//              conditions. General 2D/3D serial mesh files and finite element
//              polynomial degrees can be specified by command line options.

#include "mfem.hpp"
#include <iostream>

using namespace std;
using namespace mfem;

template <class T> void set_iters(T&) {}
template <>
void set_iters<RaptorRugeStuben>(RaptorRugeStuben & s) {
   s.SetMaxIter(1);
}

template <class O, class S> void raptorpcg(O &, S &, const Vector &, Vector &) {}
template <>
void raptorpcg<RaptorParMatrix, RaptorRugeStuben>(RaptorParMatrix & A,
                                                  RaptorRugeStuben & M,
                                                  const Vector & B, Vector & X)
{
   RaptorPCG cg(MPI_COMM_WORLD);
   cg.SetTol(1e-12);
   cg.SetMaxIter(2000);
   cg.SetPreconditioner(M);
   cg.SetOperator(A);
   cg.Mult(B, X);
}

template <class OpType>
void run(ParBilinearForm & a, ParGridFunction & x,
         ParLinearForm & b, Array<int> & boundary_dofs,
         ParMesh & mesh, bool raptor_pcg) {

   // 10. Form the linear system A X = B. This includes eliminating boundary
   //     conditions, applying AMR constraints, parallel assembly, etc.
   OpType A;

   Vector B, X;
   a.FormLinearSystem(boundary_dofs, x, b, A, X, B);

   using solver_type =
      typename std::conditional<std::is_same<OpType, HypreParMatrix>::value,
                                HypreBoomerAMG, RaptorRugeStuben>::type;
   solver_type M(A);
   set_iters(M);

   if (raptor_pcg) {
      raptorpcg(A, M, B, X);
   } else {
      // GMRESSolver cg(MPI_COMM_WORLD);
      CGSolver cg(MPI_COMM_WORLD);
      // MINRESSolver cg(MPI_COMM_WORLD);
      cg.SetRelTol(1e-12);
      cg.SetMaxIter(2000);
      cg.SetPrintLevel(1);
      cg.SetPreconditioner(M);
      cg.SetOperator(A);

      cg.Mult(B, X);
   }
   // 12. Recover the solution x as a grid function and save to file. The output
   //     can be viewed using GLVis as follows: "glvis -np <np> -m mesh -g sol"
   a.RecoverFEMSolution(X, b, x);
   x.Save("sol");
   mesh.Save("mesh");
}


int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   Hypre::Init();

   // 2. Parse command line options.
   string mesh_file = "../data/star.mesh";
   int order = 1;
   bool use_raptor = false;
   bool raptor_pcg = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
   args.AddOption(&use_raptor, "-r", "--raptor", "-h", "--hypre", "Use raptor");
   args.AddOption(&raptor_pcg, "-p", "--raptorpcg", "-f", "--mfem", "use mfem cg");
   args.ParseCheck();

   if (Mpi::Root()) {
      if (use_raptor) {
         std::cout << "Using Raptor" << std::endl;
      } else {
         std::cout << "Using Hypre" << std::endl;
      }
   }

   // 3. Read the serial mesh from the given mesh file.
   Mesh serial_mesh(mesh_file);

   // 4. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh once in parallel to increase the resolution.
   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   serial_mesh.Clear(); // the serial mesh is no longer needed
   mesh.UniformRefinement();

   // 5. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   H1_FECollection fec(order, mesh.Dimension());
   ParFiniteElementSpace fespace(&mesh, &fec);
   HYPRE_BigInt total_num_dofs = fespace.GlobalTrueVSize();
   if (Mpi::Root())
   {
      cout << "Number of unknowns: " << total_num_dofs << endl;
   }

   // 6. Extract the list of all the boundary DOFs. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.
   Array<int> boundary_dofs;
   fespace.GetBoundaryTrueDofs(boundary_dofs);

   // 7. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   ParGridFunction x(&fespace);
   x = 0.0;

   // 8. Set up the linear form b(.) corresponding to the right-hand side.
   ConstantCoefficient one(1.0);
   ParLinearForm b(&fespace);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   // 9. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
   ParBilinearForm a(&fespace);
   if (use_raptor)
      a.SetOperatorType(Operator::RAPTOR_ParCSR);
   a.AddDomainIntegrator(new DiffusionIntegrator);
   a.Assemble();

   if (use_raptor)
      run<RaptorParMatrix>(a, x, b, boundary_dofs, mesh, raptor_pcg);
   else
      run<HypreParMatrix>(a, x, b, boundary_dofs, mesh, raptor_pcg);

   return 0;
}
