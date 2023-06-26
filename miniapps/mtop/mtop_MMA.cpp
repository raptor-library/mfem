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

#include "mtop_MMA.hpp"

namespace mma
{

   MMA::MMA(int n, int m, double * x, int sx)
   {

   }

   MMA::MMA(MPI_Comm Comm,int n, int m, double * x, int sx)
   {

   }

   void MMA::Update(double* xval, double* dfdx, double* gx, double* dgdx, double* xmin, double xmax)
   {

   }

   void MMA::KKTresidual(double* xval, double* dfdx, double* gx, double* dgdx, double* xmin, double* xmax, double* norm2,
                               double* normInf)
   {

   }


    void MMA::Restart(double* xo1, double* xo2, double* U, double* L)
    {

    }

    void MMA::SetAsymptotes(double init, double decrease, double increase)
    {

    }


} // end mma namespace
