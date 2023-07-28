// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "../bilininteg.hpp"
#include "../gridfunc.hpp"
#include "../ceed/integrators/diffusionmass/diffusionmass.hpp"

using namespace std;

namespace mfem
{

void DiffusionMassIntegrator::AssemblePA(const FiniteElementSpace &fes)
{
   Mesh *mesh = fes.GetMesh();
   if (mesh->GetNE() == 0) { return; }
   if (DeviceCanUseCeed())
   {
      delete ceedOp;
      if (MQd) { ceedOp = new ceed::PADiffusionMassIntegrator(*this, fes, MQd, Qm); }
      else if (VQd) { ceedOp = new ceed::PADiffusionMassIntegrator(*this, fes, VQd, Qm); }
      else { ceedOp = new ceed::PADiffusionMassIntegrator(*this, fes, Qd, Qm); }
      return;
   }

   // Assumes tensor-product elements
   // const FiniteElement &el = *fes.GetFE(0);
   // ElementTransformation &T = *mesh->GetElementTransformation(0);
   // const IntegrationRule *ir = IntRule ? IntRule : &GetRule(*el, T);
   MFEM_ABORT("Error: DiffusionMassIntegrator::AssemblePA only implemented with"
              " libCEED");
}

void DiffusionMassIntegrator::AssemblePABoundary(const FiniteElementSpace &fes)
{
   Mesh *mesh = fes.GetMesh();
   if (mesh->GetNBE() == 0) { return; }
   if (DeviceCanUseCeed())
   {
      delete ceedOp;
      if (MQd) { ceedOp = new ceed::PADiffusionMassIntegrator(*this, fes, MQd, Qm, true); }
      else if (VQd) { ceedOp = new ceed::PADiffusionMassIntegrator(*this, fes, VQd, Qm, true); }
      else { ceedOp = new ceed::PADiffusionMassIntegrator(*this, fes, Qd, Qm, true); }
      return;
   }

   // Assumes tensor-product elements
   // const FiniteElement &el = *fes.GetBE(0);
   // ElementTransformation &T = *mesh->GetBdrElementTransformation(0);
   // const IntegrationRule *ir = IntRule ? IntRule : &GetRule(*el, T);
   MFEM_ABORT("Error: DiffusionMassIntegrator::AssemblePABoundary only implemented with"
              " libCEED");
}

void DiffusionMassIntegrator::AssembleDiagonalPA(Vector &diag)
{
   if (DeviceCanUseCeed())
   {
      if (ceedOp) { ceedOp->GetDiagonal(diag); }
   }
   else
   {
      MFEM_ABORT("Error: DiffusionMassIntegrator::AssembleDiagonalPA only"
                 " implemented with libCEED");
   }
}

void DiffusionMassIntegrator::AddMultPA(const Vector &x, Vector &y) const
{
   if (DeviceCanUseCeed())
   {
      if (ceedOp) { ceedOp->AddMult(x, y); }
   }
   else
   {
      MFEM_ABORT("Error: DiffusionMassIntegrator::AddMultPA only implemented with"
                 " libCEED");
   }
}

}