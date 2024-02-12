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

#ifndef MFEM_RAPTOR
#define MFEM_RAPTOR

#include "../config/config.hpp"
#include "linalg/vector.hpp"

#ifdef MFEM_USE_RAPTOR
#ifdef MFEM_USE_MPI

#include <memory>

#include <raptor/core/matrix.hpp>
#include <raptor/core/par_matrix.hpp>
#include <raptor/core/par_vector.hpp>
#include <raptor/ruge_stuben/par_ruge_stuben_solver.hpp>

#include "handle.hpp"
#include "hypre.hpp"

namespace mfem {

class RaptorParVector : public Vector {
public:
	RaptorParVector() : owns_vec(false), vec(NULL) {}

	RaptorParVector(raptor::index_t gsize, int lsize);
	RaptorParVector(const RaptorParVector & other);
	RaptorParVector(RaptorParVector && other);

	void WrapMemoryRead(const Memory<double> & mem);
	void WrapMemoryReadWrite(Memory<double> & mem);

	operator raptor::ParVector*() const { return vec; }

private:
	void _SetDataAndSize_();

	bool owns_vec;
	raptor::ParVector *vec;
};


class RaptorParMatrix : public Operator {
public:
	RaptorParMatrix();
	// block-diagonal square parallel matrix
	RaptorParMatrix(MPI_Comm comm, HYPRE_BigInt glob_size,
	                HYPRE_BigInt *row_starts, SparseMatrix *diag,
	                Operator::Type tid);
	explicit RaptorParMatrix(raptor::ParMatrix * m, bool owner = true);
	explicit RaptorParMatrix(const HypreParMatrix *ha,
	                         Operator::Type tid = Operator::RAPTOR_ParCSR);
	~RaptorParMatrix();

	operator raptor::ParMatrix*() const { return mat; }
	operator raptor::ParCSRMatrix*() const;
	operator raptor::ParBSRMatrix*() const;

	void Mult(double a, const Vector & x, double b, Vector & y) const;
	void Mult(const Vector & x, Vector & y) const override {
		Mult(1., x, 0., y);
	};
	Operator::Type GetType() const;

	void MakeRef(const RaptorParMatrix & m);

	RaptorParMatrix * EliminateRowsCols(const Array<int> & rows_cols);

	int GetGlobalNumRows() const {
		return mat->global_num_rows;
	}

	int GetGlobalNumCols() const {
		return mat->global_num_cols;
	}

private:
	void ConstructBlockDiagCSR(MPI_Comm comm, HYPRE_BigInt glob_size,
	                           HYPRE_BigInt *row_starts, SparseMatrix *diag);
	void ConstructBlocKDiagBSR(MPI_Comm comm, HYPRE_BigInt glob_size,
	                           HYPRE_BigInt *row_starts, SparseMatrix *diag);
	void CopyCSR(raptor::CSRMatrix & csr, const SparseMatrix & diag);
	raptor::ParMatrix * Convert(const HypreParMatrix & ha,
	                            Operator::Type tid);
	void Destroy();

	raptor::ParMatrix *mat;
	mutable RaptorParVector *X, *Y;
	bool owns_mat;
};


class RaptorSolver : public Solver
{
public:
	RaptorSolver(raptor::ParMultilevel *);
	RaptorSolver(raptor::ParMultilevel *, const RaptorParMatrix *);

	virtual void Setup() const;

	virtual void Mult(const RaptorParVector &, RaptorParVector&) const;
	virtual void Mult(const Vector &, Vector &) const;
	using Operator::Mult;
	void SetOperator(const Operator &) override
	{
		mfem_error("RaptorSolvers do not support SetOperator!");
	}

	virtual ~RaptorSolver();

	void WrapVectors(const Vector &, Vector &) const;


protected:
	const RaptorParMatrix *A;
	mutable RaptorParVector *B, *X;
	mutable bool setup_called;
	mutable raptor::ParMultilevel *target;
};


class RaptorRugeStuben : public RaptorSolver
{
public:
	RaptorRugeStuben();
	RaptorRugeStuben(const RaptorParMatrix &);
};

// Return the matrix P^t * A * P
RaptorParMatrix *RAP(RaptorParMatrix *A, RaptorParMatrix *P);

void EliminateBC(RaptorParMatrix & A, RaptorParMatrix & Ae,
                 const Array<int> & ess_dof_list,
                 const Vector & X, Vector & B);
}

#endif // MFEM_USE_MPI
#endif // MFEM_USE_RAPTOR

#endif
