// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "../tmop.hpp"
#include "tmop_pa.hpp"
#include "../tmop_tools.hpp"
#include "../../general/forall.hpp"
#include "../../linalg/kernels.hpp"

namespace mfem
{

MFEM_REGISTER_TMOP_KERNELS(double, MinDetJpr_Kernel_2D,
                           const int NE,
                           const Array<double> &b_,
                           const Array<double> &g_,
                           const Vector &x_,
                           Vector &DetJ,
                           const int d1d,
                           const int q1d)
{
   constexpr int DIM = 2;
   constexpr int NBZ = 1;

   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;

   const auto b = Reshape(b_.Read(), Q1D, D1D);
   const auto g = Reshape(g_.Read(), Q1D, D1D);
   const auto X = Reshape(x_.Read(), D1D, D1D, DIM, NE);

   auto E = Reshape(DetJ.Write(), Q1D, Q1D, NE);

   MFEM_FORALL_2D(e, NE, Q1D, Q1D, NBZ,
   {
      const int D1D = T_D1D ? T_D1D : d1d;
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      constexpr int NBZ = 1;
      constexpr int MQ1 = T_Q1D ? T_Q1D : T_MAX;
      constexpr int MD1 = T_D1D ? T_D1D : T_MAX;

      MFEM_SHARED double BG[2][MQ1*MD1];
      MFEM_SHARED double XY[2][NBZ][MD1*MD1];
      MFEM_SHARED double DQ[4][NBZ][MD1*MQ1];
      MFEM_SHARED double QQ[4][NBZ][MQ1*MQ1];

      kernels::internal::LoadX<MD1,NBZ>(e,D1D,X,XY);
      kernels::internal::LoadBG<MD1,MQ1>(D1D,Q1D,b,g,BG);

      kernels::internal::GradX<MD1,MQ1,NBZ>(D1D,Q1D,BG,XY,DQ);
      kernels::internal::GradY<MD1,MQ1,NBZ>(D1D,Q1D,BG,DQ,QQ);

      MFEM_FOREACH_THREAD(qy,y,Q1D)
      {
         MFEM_FOREACH_THREAD(qx,x,Q1D)
         {
            double J[4];
            kernels::internal::PullGrad<MQ1,NBZ>(qx,qy,QQ,J);
            E(qx,qy,e) = kernels::Det<2>(J);
         }
      }
   });

   return DetJ.Min();
}

double TMOPNewtonSolver::MinDetJpr_2D(const FiniteElementSpace *fes,
                                      const Vector &X) const
{
   const ElementDofOrdering ordering = ElementDofOrdering::LEXICOGRAPHIC;
   const Operator *R = fes->GetElementRestriction(ordering);
   Vector XE(R->Height(), Device::GetDeviceMemoryType());
   XE.UseDevice(true);
   R->Mult(X, XE);

   const DofToQuad &maps = fes->GetFE(0)->GetDofToQuad(ir, DofToQuad::TENSOR);
   const int NE = fes->GetMesh()->GetNE();
   const int NQ = ir.GetNPoints();
   const int D1D = maps.ndof;
   const int Q1D = maps.nqpt;
   const int id = (D1D << 4 ) | Q1D;
   const Array<double> &B = maps.B;
   const Array<double> &G = maps.G;

   Vector E(NE*NQ);
   E.UseDevice(true);

   MFEM_LAUNCH_TMOP_KERNEL(MinDetJpr_Kernel_2D,id,NE,B,G,XE,E);
}

} // namespace mfem
