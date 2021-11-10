#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include "../linalg/dtensor.hpp"

using namespace std;
using namespace mfem;

int problem;
void velocity_function(const Vector &x, Vector &v);
double u0_function(const Vector &x);
double inflow_function(const Vector &x);

Vector bb_min, bb_max;

void AddDGIntegrators(BilinearForm &k, VectorCoefficient &velocity)
{
   double alpha = 1.0;
   double beta = -0.5;
   k.AddDomainIntegrator(new ConvectionIntegrator(velocity, -alpha));
   k.AddDomainIntegrator(new MassIntegrator);
   k.AddInteriorFaceIntegrator(
      new TransposeIntegrator(new DGTraceIntegrator(velocity, alpha, beta)));
   k.AddBdrFaceIntegrator(
      new TransposeIntegrator(new DGTraceIntegrator(velocity, alpha, beta)));
   // k.AddInteriorFaceIntegrator(new DGTraceIntegrator(velocity, alpha, beta));
   // k.AddBdrFaceIntegrator(new DGTraceIntegrator(velocity, alpha, beta));
}

void SaveSolution(const std::string &fname, GridFunction &gf)
{
   ofstream osol(fname);
   osol.precision(16);
   gf.Save(osol);
}

Mesh *oriented_mesh()
{
   static const int dim = 3;
   static const int nv = 12;
   static const int nel = 2;
   Mesh *mesh = new Mesh(dim, nv, nel);
   double x[dim];
   x[0] = 0.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 0.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 1.0;   x[2] = 2.0;
   mesh->AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;   x[2] = 1.0;
   mesh->AddVertex(x);

   //EAST
   x[0] = 2.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 2.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 2.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh->AddVertex(x);
   x[0] = 2.0;   x[1] = 1.0;   x[2] = 1.0;
   mesh->AddVertex(x);
   //WEST
   // x[0] = -1.0;   x[1] = 0.0;   x[2] = 0.0;
   // mesh->AddVertex(x);
   // x[0] = -1.0;   x[1] = 1.0;   x[2] = 0.0;
   // mesh->AddVertex(x);
   // x[0] = -1.0;   x[1] = 0.0;   x[2] = 1.0;
   // mesh->AddVertex(x);
   // x[0] = -1.0;   x[1] = 1.0;   x[2] = 1.0;
   // mesh->AddVertex(x);

   int el[8];
   el[0] = 0;
   el[1] = 1;
   el[2] = 2;
   el[3] = 3;
   el[4] = 4;
   el[5] = 5;
   el[6] = 6;
   el[7] = 7;
   mesh->AddHex(el);
   // ELEM1 WEST
   // orientation 3 WEST/EAST OK
   // el[0] = 8;
   // el[1] = 0;
   // el[2] = 3;
   // el[3] = 9;
   // el[4] = 10;
   // el[5] = 4;
   // el[6] = 7;
   // el[7] = 11;
   // orientation 3 WEST/SOUTH OK
   // el[0] = 0;
   // el[1] = 3;
   // el[2] = 9;
   // el[3] = 8;
   // el[4] = 4;
   // el[5] = 7;
   // el[6] = 11;
   // el[7] = 10;
   // orientation 3 WEST/NORTH OK
   // el[0] = 8;
   // el[1] = 9;
   // el[2] = 0;
   // el[3] = 3;
   // el[4] = 10;
   // el[5] = 11;
   // el[6] = 4;
   // el[7] = 7;
   // orientation 5 WEST/TOP OK
   // el[0] = 10;
   // el[1] = 8;
   // el[2] = 9;
   // el[3] = 11;
   // el[4] = 4;
   // el[5] = 0;
   // el[6] = 3;
   // el[7] = 7;
   // orientation 3 WEST/TOP OK
   // el[0] = 8;
   // el[1] = 9;
   // el[2] = 11;
   // el[3] = 10;
   // el[4] = 0;
   // el[5] = 3;
   // el[6] = 7;
   // el[7] = 4;
   // orientation 3 WEST/BOTTOM OK
   // el[0] = 4;
   // el[1] = 7;
   // el[2] = 3;
   // el[3] = 0;
   // el[4] = 10;
   // el[5] = 11;
   // el[6] = 9;
   // el[7] = 8;

   // ELEM1 EAST
   // orientation 3 EAST/WEST OK
   el[0] = 1;
   el[1] = 8;
   el[2] = 9;
   el[3] = 2;
   el[4] = 5;
   el[5] = 10;
   el[6] = 11;
   el[7] = 6;
   // orientation 1 EAST/WEST OK
   // el[0] = 5;
   // el[1] = 10;
   // el[2] = 8;
   // el[3] = 1;
   // el[4] = 6;
   // el[5] = 11;
   // el[6] = 9;
   // el[7] = 2;
   // orientation 7 EAST/WEST OK
   // el[0] = 6;
   // el[1] = 11;
   // el[2] = 10;
   // el[3] = 5;
   // el[4] = 2;
   // el[5] = 9;
   // el[6] = 8;
   // el[7] = 1;
   // orientation 5 EAST/WEST OK
   // el[0] = 2;
   // el[1] = 9;
   // el[2] = 11;
   // el[3] = 6;
   // el[4] = 1;
   // el[5] = 8;
   // el[6] = 10;
   // el[7] = 5;
   // orientation 3 EAST/EAST OK
   // el[0] = 9;
   // el[1] = 2;
   // el[2] = 1;
   // el[3] = 8;
   // el[4] = 11;
   // el[5] = 6;
   // el[6] = 5;
   // el[7] = 10;
   // orientation 1 EAST/EAST OK
   // el[0] = 8;
   // el[1] = 1;
   // el[2] = 5;
   // el[3] = 10;
   // el[4] = 9;
   // el[5] = 2;
   // el[6] = 6;
   // el[7] = 11;
   // orientation 7 EAST/EAST OK
   // el[0] = 10;
   // el[1] = 5;
   // el[2] = 6;
   // el[3] = 11;
   // el[4] = 8;
   // el[5] = 1;
   // el[6] = 2;
   // el[7] = 9;
   // orientation 5 EAST/EAST OK
   // el[0] = 11;
   // el[1] = 6;
   // el[2] = 2;
   // el[3] = 9;
   // el[4] = 10;
   // el[5] = 5;
   // el[6] = 1;
   // el[7] = 8;
   // orientation 3 EAST/TOP OK
   // el[0] = 9;
   // el[1] = 8;
   // el[2] = 10;
   // el[3] = 11;
   // el[4] = 2;
   // el[5] = 1;
   // el[6] = 5;
   // el[7] = 6;
   // orientation 1 EAST/TOP OK
   // el[0] = 8;
   // el[1] = 10;
   // el[2] = 11;
   // el[3] = 9;
   // el[4] = 1;
   // el[5] = 5;
   // el[6] = 6;
   // el[7] = 2;
   // orientation 7 EAST/TOP OK
   // el[0] = 10;
   // el[1] = 11;
   // el[2] = 9;
   // el[3] = 8;
   // el[4] = 5;
   // el[5] = 6;
   // el[6] = 2;
   // el[7] = 1;
   // orientation 5 EAST/TOP OK
   // el[0] = 11;
   // el[1] = 9;
   // el[2] = 8;
   // el[3] = 10;
   // el[4] = 6;
   // el[5] = 2;
   // el[6] = 1;
   // el[7] = 5;
   // orientation 5 EAST/BOTTOM OK
   // el[0] = 5;
   // el[1] = 1;
   // el[2] = 2;
   // el[3] = 6;
   // el[4] = 10;
   // el[5] = 8;
   // el[6] = 9;
   // el[7] = 11;
   // orientation 7 EAST/BOTTOM OK
   // el[0] = 1;
   // el[1] = 2;
   // el[2] = 6;
   // el[3] = 5;
   // el[4] = 8;
   // el[5] = 9;
   // el[6] = 11;
   // el[7] = 10;
   // orientation 1 EAST/BOTTOM OK
   // el[0] = 2;
   // el[1] = 6;
   // el[2] = 5;
   // el[3] = 1;
   // el[4] = 9;
   // el[5] = 11;
   // el[6] = 10;
   // el[7] = 8;
   // orientation 3 EAST/BOTTOM OK
   // el[0] = 6;
   // el[1] = 5;
   // el[2] = 1;
   // el[3] = 2;
   // el[4] = 11;
   // el[5] = 10;
   // el[6] = 8;
   // el[7] = 9;
   // orientation 3 EAST/SOUTH OK
   // el[0] = 2;
   // el[1] = 1;
   // el[2] = 8;
   // el[3] = 9;
   // el[4] = 6;
   // el[5] = 5;
   // el[6] = 10;
   // el[7] = 11;
   // orientation 5 EAST/SOUTH OK
   // el[0] = 6;
   // el[1] = 2;
   // el[2] = 9;
   // el[3] = 11;
   // el[4] = 5;
   // el[5] = 1;
   // el[6] = 8;
   // el[7] = 10;
   // orientation 7 EAST/SOUTH OK
   // el[0] = 5;
   // el[1] = 6;
   // el[2] = 11;
   // el[3] = 10;
   // el[4] = 1;
   // el[5] = 2;
   // el[6] = 9;
   // el[7] = 8;
   // orientation 1 EAST/SOUTH OK
   // el[0] = 1;
   // el[1] = 5;
   // el[2] = 10;
   // el[3] = 8;
   // el[4] = 2;
   // el[5] = 6;
   // el[6] = 11;
   // el[7] = 9;
   // orientation 3 EAST/NORTH OK
   // el[0] = 8;
   // el[1] = 9;
   // el[2] = 2;
   // el[3] = 1;
   // el[4] = 10;
   // el[5] = 11;
   // el[6] = 6;
   // el[7] = 5;
   // orientation 5 EAST/NORTH OK
   // el[0] = 9;
   // el[1] = 11;
   // el[2] = 6;
   // el[3] = 2;
   // el[4] = 8;
   // el[5] = 10;
   // el[6] = 5;
   // el[7] = 1;
   // orientation 7 EAST/NORTH OK
   // el[0] = 11;
   // el[1] = 10;
   // el[2] = 5;
   // el[3] = 6;
   // el[4] = 9;
   // el[5] = 8;
   // el[6] = 1;
   // el[7] = 2;
   // orientation 1 EAST/NORTH OK
   // el[0] = 10;
   // el[1] = 8;
   // el[2] = 1;
   // el[3] = 5;
   // el[4] = 11;
   // el[5] = 9;
   // el[6] = 2;
   // el[7] = 6;
   mesh->AddHex(el);

   mesh->FinalizeHexMesh(true);
   mesh->GenerateBoundaryElements();
   mesh->Finalize();
   return mesh;
}

Mesh *skewed_mesh_2d()
{
   static const int dim = 2;
   static const int nv = 4;
   static const int nel = 1;
   Mesh *mesh = new Mesh(dim, nv, nel);
   double x[2];
   x[0] = 0.0;   x[1] = 0.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;
   mesh->AddVertex(x);
   x[0] = 2.0;   x[1] = 1.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 2.0;
   mesh->AddVertex(x);
   int el[4];
   el[0] = 0;
   el[1] = 1;
   el[2] = 2;
   el[3] = 3;
   mesh->AddQuad(el);
   mesh->FinalizeQuadMesh(true);
   mesh->GenerateBoundaryElements();
   mesh->Finalize();
   return mesh;
}

Mesh *skewed_mesh_3d()
{
   static const int dim = 3;
   static const int nv = 8;
   static const int nel = 1;
   Mesh *mesh = new Mesh(dim, nv, nel);
   double x[dim];
   x[0] = 0.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh->AddVertex(x);
   x[0] = 0.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh->AddVertex(x);
   x[0] = 1.0;   x[1] = 2.0;   x[2] = 1.0;
   mesh->AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;   x[2] = 1.0;
   mesh->AddVertex(x);

   int el[8];
   el[0] = 0;
   el[1] = 1;
   el[2] = 2;
   el[3] = 3;
   el[4] = 4;
   el[5] = 5;
   el[6] = 6;
   el[7] = 7;
   mesh->AddHex(el);

   mesh->FinalizeHexMesh(true);
   mesh->GenerateBoundaryElements();
   mesh->Finalize();
   return mesh;
}

Mesh nc_2dmesh()
{
   static const int dim = 2;
   static const int nv = 7;
   static const int nel = 2;
   Mesh mesh(dim, nv, nel);
   double x[dim];
   x[0] = 0.5;   x[1] = 0.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;
   mesh.AddVertex(x);
   x[0] = 0.0;   x[1] = 0.5;
   mesh.AddVertex(x);
   x[0] = 0.5;   x[1] = 0.5;
   mesh.AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;
   mesh.AddVertex(x);
   x[0] = 0.5;   x[1] = 1.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 1.0;
   mesh.AddVertex(x);

   int el[4];
   el[0] = 0;
   el[1] = 1;
   el[2] = 6;
   el[3] = 5;
   mesh.AddQuad(el);

   el[0] = 2;
   el[1] = 3;
   el[2] = 5;
   el[3] = 4;
   mesh.AddQuad(el);

   mesh.FinalizeQuadMesh(true);
   mesh.GenerateBoundaryElements();
   mesh.Finalize();
   return mesh;
}

Mesh rotated_2dmesh(int face_perm_1, int face_perm_2)
{
   static const int dim = 2;
   static const int nv = 6;
   static const int nel = 2;
   Mesh mesh(dim, nv, nel);
   double x[dim];
   x[0] = 0.0;   x[1] = 0.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;
   mesh.AddVertex(x);
   x[0] = 2.0;   x[1] = 0.0;
   mesh.AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 1.0;
   mesh.AddVertex(x);
   x[0] = 2.0;   x[1] = 1.0;
   mesh.AddVertex(x);
   int el[4];
   el[0] = 0;
   el[1] = 1;
   el[2] = 4;
   el[3] = 3;
   std::rotate(&el[0], &el[face_perm_1], &el[3] + 1);

   mesh.AddQuad(el);

   el[0] = 1;
   el[1] = 2;
   el[2] = 5;
   el[3] = 4;
   std::rotate(&el[0], &el[face_perm_2], &el[3] + 1);
   mesh.AddQuad(el);

   mesh.FinalizeQuadMesh(true);
   mesh.GenerateBoundaryElements();
   mesh.Finalize();
   return mesh;
}

void rotate_3d_vertices(int *v, int ref_face, int rot)
{
   std::vector<int> face_1, face_2;

   switch (ref_face/2)
   {
      case 0:
         face_1 = {v[0], v[1], v[2], v[3]};
         face_2 = {v[4], v[5], v[6], v[7]};
         break;
      case 1:
         face_1 = {v[1], v[5], v[6], v[2]};
         face_2 = {v[0], v[4], v[7], v[3]};
         break;
      case 2:
         face_1 = {v[4], v[5], v[1], v[0]};
         face_2 = {v[7], v[6], v[2], v[3]};
         break;
   }
   if (ref_face % 2 == 0)
   {
      std::reverse(face_1.begin(), face_1.end());
      std::reverse(face_2.begin(), face_2.end());
      std::swap(face_1, face_2);
   }

   std::rotate(face_1.begin(), face_1.begin() + rot, face_1.end());
   std::rotate(face_2.begin(), face_2.begin() + rot, face_2.end());

   for (int i=0; i<4; ++i)
   {
      v[i] = face_1[i];
      v[i+4] = face_2[i];
   }
}

Mesh rotated_3dmesh(int face_perm_1, int face_perm_2)
{
   static const int dim = 3;
   static const int nv = 12;
   static const int nel = 2;
   Mesh mesh(dim, nv, nel);
   double x[dim];
   x[0] = 0.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh.AddVertex(x);
   x[0] = 2.0;   x[1] = 0.0;   x[2] = 0.0;
   mesh.AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh.AddVertex(x);
   x[0] = 2.0;   x[1] = 1.0;   x[2] = 0.0;
   mesh.AddVertex(x);
   x[0] = 0.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh.AddVertex(x);
   x[0] = 2.0;   x[1] = 0.0;   x[2] = 1.0;
   mesh.AddVertex(x);
   x[0] = 0.0;   x[1] = 1.0;   x[2] = 1.0;
   mesh.AddVertex(x);
   x[0] = 1.0;   x[1] = 1.0;   x[2] = 1.0;
   mesh.AddVertex(x);
   x[0] = 3.0;   x[1] = 1.0;   x[2] = 1.0;
   mesh.AddVertex(x);

   int el[8];

   el[0] = 0;
   el[1] = 1;
   el[2] = 4;
   el[3] = 3;
   el[4] = 6;
   el[5] = 7;
   el[6] = 10;
   el[7] = 9;
   rotate_3d_vertices(el, face_perm_1/4, face_perm_1%4);
   mesh.AddHex(el);

   el[0] = 1;
   el[1] = 2;
   el[2] = 5;
   el[3] = 4;
   el[4] = 7;
   el[5] = 8;
   el[6] = 11;
   el[7] = 10;
   rotate_3d_vertices(el, face_perm_2/4, face_perm_2%4);
   mesh.AddHex(el);

   mesh.FinalizeHexMesh(true);
   mesh.GenerateBoundaryElements();
   mesh.Finalize();
   return mesh;
}

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   MPI_Session mpi;

   // 1. Parse command-line options.
   problem = 0;
   const char *mesh_file = "../data/inline-quad.mesh";
   int ref_levels = 0;
   int order = 3;
   const char *device_config = "cpu";
   bool visualization = true;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&problem, "-p", "--problem",
                  "Problem setup to use. See options in velocity_function().");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   int pi(0), pj(0);
   args.AddOption(&pi, "-pi", "--permi", "Permutation i.");
   args.AddOption(&pj, "-pj", "--permj", "Permutation j.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root())
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (mpi.Root())
   {
      args.PrintOptions(cout);
   }

   Device device(device_config);
   device.Print();

   //Creating custom mesh
   // Mesh *mesh_ptr = oriented_mesh();
   // Mesh *mesh_ptr = skewed_mesh_3d();
   // Mesh& mesh = *mesh_ptr;

   Mesh mesh(mesh_file, 1, 1);
   // Mesh mesh = nc_2dmesh();
   // Mesh mesh = Mesh::MakeCartesian3D(2,1,1,Element::HEXAHEDRON);
   // Not working: (8,0),(10,0),(12,0),(14,0),(16,0),(17,0),(18,0),(19,0),(20,0)
   // (21,0),(22,0),(23,0)
   // (8,1),(10,1),(12,1),
   // Mesh mesh = rotated_2dmesh(pi,pj);
   // Mesh mesh = rotated_2dmesh(pi,pj);
   int dim = mesh.Dimension();

   mesh.EnsureNodes();
   // mesh.SetCurvature(3,true);
   // Array<Refinement> refs(1);
   // refs[0].index = 0;
   // refs[0].ref_type = Refinement::XYZ;
   // mesh.GeneralRefinement(refs,1,4);

   // mesh.UniformRefinement();
   for (int lev = 0; lev < ref_levels; lev++)
   {
      // mesh.UniformRefinement();
      mesh.RandomRefinement(0.6, false, 1, 4);
   }
   mesh.GetBoundingBox(bb_min, bb_max, max(order, 1));

   int* partitioning = mesh.GeneratePartitioning(mpi.WorldSize());
   // ParMesh pmesh(MPI_COMM_WORLD, mesh, partitioning);
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   for (int i = 0; i < mesh.GetNE(); i++)
   {
      partitioning[i] = i*mpi.WorldSize()/mesh.GetNE();
   }
   // for (int lev = 0; lev < ref_levels; lev++)
   // {
   //    // mesh.UniformRefinement();
   //    pmesh.RandomRefinement(0.6, false, 1, 4);
   // }
   pmesh.ExchangeFaceNbrData();
   DG_FECollection fec(order, dim, BasisType::GaussLobatto);
   // H1_FECollection fec(order, dim);
   FiniteElementSpace serial_fes(&mesh, &fec);
   ParFiniteElementSpace fes(&pmesh, &fec);

   cout << "Number of unknowns: " << fes.GetVSize() << endl;

   Vector velocity_vector(dim);
   for (int i = 0; i < dim; ++i)
   {
      // velocity_vector[i] = 0.0;
      velocity_vector[i] = -1;
   }
   // velocity_vector[0] = 1.0;
   VectorConstantCoefficient velocity(velocity_vector);
   // VectorFunctionCoefficient velocity(dim, velocity_function);
   FunctionCoefficient inflow(inflow_function);
   FunctionCoefficient u0(u0_function);

   ParBilinearForm k_test(&fes), k_ref(&fes);
   // k_ref.SetAssemblyLevel(AssemblyLevel::PARTIAL);
   // k_test.SetAssemblyLevel(AssemblyLevel::FULL);
   k_test.SetAssemblyLevel(AssemblyLevel::PARTIAL);
   // k_test.SetAssemblyLevel(AssemblyLevel::ELEMENT);

   AddDGIntegrators(k_test, velocity);
   AddDGIntegrators(k_ref, velocity);

   tic_toc.Clear();
   tic_toc.Start();
   k_test.Assemble();
   k_test.Finalize();
   tic_toc.Stop();
   cout << "test assembly time: " << tic_toc.RealTime() << " sec." << endl;
   tic_toc.Clear();
   tic_toc.Start();
   k_ref.Assemble();
   k_ref.Finalize();
   tic_toc.Stop();
   cout << "ref assembly time: " << tic_toc.RealTime() << " sec." << endl;

   // std::cout << "FA matrix:" << std::endl;
   // std::cout  << k_ref.SpMat() << std::endl;
   ParFiniteElementSpace& test_fes = fes;

   ParGridFunction fx(&test_fes), fy(&test_fes);
   FunctionCoefficient coeffx( [&](const Vector &x) { return x[0]; } );
   FunctionCoefficient coeffy( [&](const Vector &x) { return x[1]; } );

   fx.ProjectCoefficient(coeffx);
   fy.ProjectCoefficient(coeffy);

   auto face_restr = test_fes.GetFaceRestriction(
                        ElementDofOrdering::LEXICOGRAPHIC,
                        FaceType::Interior,
                        L2FaceValues::DoubleValued);

   Vector resx(face_restr->Height());
   Vector resy(face_restr->Height());
   const int face_dofs = dim==2? order + 1 : (order+1)*(order+1);
   const int nf = test_fes.GetNFbyType(FaceType::Interior);

   face_restr->Mult(fx,resx);
   face_restr->Mult(fy,resy);

   const auto tx = Reshape(resx.HostReadWrite(),face_dofs,2,nf);
   const auto ty = Reshape(resy.HostReadWrite(),face_dofs,2,nf);

   if (mpi.WorldRank()!=0)
   {
      int val;
      MPI_Recv(&val,1,MPI_INT,mpi.WorldRank()-1,0,MPI_COMM_WORLD,NULL);
   }
   std::cout << "My rank is: " << mpi.WorldRank() << std::endl;
   double errorx = 0.0;
   double errory = 0.0;
   const double tol = 1e-12;
   for (int face = 0; face < nf; face++)
   {
      double face_error_x = 0.0;
      double face_error_y = 0.0;
      for (int dof = 0; dof < face_dofs; dof++)
      {
         const double loc_error_x = std::abs(tx(dof,0,face)-tx(dof,1,face));
         face_error_x = max(face_error_x, loc_error_x);
         const double loc_error_y = std::abs(ty(dof,0,face)-ty(dof,1,face));
         face_error_y = max(face_error_y, loc_error_y);
      }
      errorx = max(errorx, face_error_x);
      errory = max(errory, face_error_y);
      if (face_error_x > tol || face_error_y > tol)
         std::cout << "Error on face " << face << std::endl
                   << "face_error_x = " << face_error_x << " , face_error_y = " << face_error_y <<
                   std::endl
                   << pmesh.GetFaceInformation(face) << std::endl;
   }
   std::cout << "errorx = " << errorx << " , errory = " << errory << std::endl;
   // for (int i = 0; i < resx.Size(); i++)
   // {
   //    std::cout << resx[i] << ", " << resy[i] << std::endl;
   // }
   // std::cout << std::endl;
   if (mpi.WorldRank()!=mpi.WorldSize()-1)
   {
      int val = 0;
      MPI_Send(&val,1,MPI_INT,mpi.WorldRank()+1,0,MPI_COMM_WORLD);
   }

   GridFunction serial_u(&serial_fes);
   serial_u.Randomize(1);

   // ParGridFunction u(&fes), r_test(&fes), r_ref(&fes), diff(&fes);
   ParGridFunction u(&pmesh, &serial_u, partitioning), r_test(&fes), r_ref(&fes),
                   diff(&fes);
   // delete[] partitioning;
   // u.ProjectCoefficient(u0);
   // u.Randomize(1);
   // u = 1.0;
   Array<int> bdofs;
   OperatorHandle A_ref;
   k_ref.FormSystemMatrix(bdofs,A_ref);

   k_test.Mult(u, r_test);
   A_ref->Mult(u, r_ref);

   diff = r_test;
   diff -= r_ref;

   std::cout << "Ref norm on rank " << mpi.WorldRank() << " : " << r_ref.Norml2()
             << '\n';
   std::cout << "Test norm on rank " << mpi.WorldRank() << " : " << r_test.Norml2()
             << '\n';
   std::cout << "Difference norm on rank " << mpi.WorldRank() << " : " <<
             diff.Norml2() << '\n';

   double norm_ref = sqrt(InnerProduct(MPI_COMM_WORLD,r_ref,r_ref));
   double norm_test = sqrt(InnerProduct(MPI_COMM_WORLD,r_test,r_test));
   std::cout << "Ref l2 norm: " << norm_ref << '\n';
   std::cout << "Test l2 norm: " << norm_test << '\n';


   ConstantCoefficient coeff(0);
   std::cout << "Ref L2 norm: " << r_ref.ComputeL2Error(coeff) << '\n';
   std::cout << "Test L2 norm: " << r_test.ComputeL2Error(coeff) << '\n';

   {
      pmesh.Save("ex9ppa");
      u.Save("init");
      diff.Save("error");
      r_ref.Save("ref");
      r_test.Save("test");
      ofstream omesh("ex9pa.mesh");
      omesh.precision(precision);
      mesh.Print(omesh);
      // ofstream osol("ex9-init.gf");
      // osol.precision(precision);
      // u.Save(osol);

      // SaveSolution("resid_error.gf", diff);
      // SaveSolution("resid_ref.gf", r_ref);
      // SaveSolution("resid_test.gf", r_test);
   }

   return 0;
}

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      double center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   switch (problem)
   {
      case 0:
      {
         // Translations in 1D, 2D, and 3D
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = sqrt(2./3.); v(1) = sqrt(1./3.); break;
            case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.);
               break;
         }
         break;
      }
      case 1:
      case 2:
      {
         // Clockwise rotation in 2D around the origin
         const double w = M_PI/2;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = w*X(1); v(1) = -w*X(0); break;
            case 3: v(0) = w*X(1); v(1) = -w*X(0); v(2) = 0.0; break;
         }
         break;
      }
      case 3:
      {
         // Clockwise twisting rotation in 2D around the origin
         const double w = M_PI/2;
         double d = max((X(0)+1.)*(1.-X(0)),0.) * max((X(1)+1.)*(1.-X(1)),0.);
         d = d*d;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = d*w*X(1); v(1) = -d*w*X(0); break;
            case 3: v(0) = d*w*X(1); v(1) = -d*w*X(0); v(2) = 0.0; break;
         }
         break;
      }
   }
}

// Initial condition
double u0_function(const Vector &x)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      double center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   switch (problem)
   {
      case 0:
      case 1:
      {
         switch (dim)
         {
            case 1:
               return exp(-40.*pow(X(0)-0.5,2));
            case 2:
            case 3:
            {
               double rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
               if (dim == 3)
               {
                  const double s = (1. + 0.25*cos(2*M_PI*X(2)));
                  rx *= s;
                  ry *= s;
               }
               return ( erfc(w*(X(0)-cx-rx))*erfc(-w*(X(0)-cx+rx)) *
                        erfc(w*(X(1)-cy-ry))*erfc(-w*(X(1)-cy+ry)) )/16;
            }
         }
      }
      case 2:
      {
         double x_ = X(0), y_ = X(1), rho, phi;
         rho = hypot(x_, y_);
         phi = atan2(y_, x_);
         return pow(sin(M_PI*rho),2)*sin(3*phi);
      }
      case 3:
      {
         const double f = M_PI;
         return sin(f*X(0))*sin(f*X(1));
      }
   }
   return 0.0;
}

// Inflow boundary condition (zero for the problems considered in this example)
double inflow_function(const Vector &x)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3: return 0.0;
   }
   return 0.0;
}
