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

#include "../config/config.hpp"

#ifdef MFEM_USE_MPI
#ifdef MFEM_USE_RAPTOR

#include <fstream>

#include <raptor/core/par_matrix.hpp>
#include <raptor/external/hypre_wrapper.hpp>
#include <raptor/core/matrix.hpp>
#include <raptor/core/types.hpp>
#include <raptor/krylov/par_cg.hpp>
#include <raptor/ruge_stuben/ruge_stuben_solver.hpp>

#include "hypre.hpp"
#include "linalg.hpp"
#include "general/error.hpp"


namespace {
template <class T>
void host_copy(T * dst, const mfem::Memory<T> & src, std::size_t n)
{
	std::memcpy(dst, src.Read(mfem::MemoryClass::HOST, n), n * sizeof(T));
}
}


namespace mfem {

void RaptorParVector::_SetDataAndSize_()
{
	SetDataAndSize(vec->local.data(), vec->local_n);
}

RaptorParVector::RaptorParVector(raptor::index_t gsize, int lsize) :
	owns_vec(1), vec(new raptor::ParVector(gsize, lsize))
{
	_SetDataAndSize_();
}


void RaptorParVector::WrapMemoryRead(const Memory<double> & mem)
{
	MFEM_ASSERT(mem.Capacity() >= size, "");

	data.Delete();
	vec->local.set_base(const_cast<double*>(mem.Read(MemoryClass::HOST, size)));
	data.MakeAlias(mem, 0, size);
}


void RaptorParVector::WrapMemoryReadWrite(Memory<double> & mem)
{
	MFEM_ASSERT(mem.Capacity() >= size, "");

	data.Delete();
	vec->local.set_base(mem.ReadWrite(MemoryClass::HOST, size));
	data.MakeAlias(mem, 0, size);
}

RaptorParMatrix::RaptorParMatrix()
	: Operator(0, 0), mat(NULL), X(NULL), Y(NULL), owns_mat(false) {}


RaptorParMatrix::RaptorParMatrix(raptor::ParMatrix *m, bool owner)
	: Operator(m->on_proc->n_rows, m->on_proc->n_cols), mat(m), X(NULL), Y(NULL),
	  owns_mat(owner) {}



RaptorParMatrix::RaptorParMatrix(MPI_Comm comm, HYPRE_BigInt glob_size,
                                 HYPRE_BigInt *row_starts, SparseMatrix *diag,
                                 Operator::Type tid)
	: Operator(diag->Height(), diag->Width()),
	  X(NULL), Y(NULL),
	  owns_mat(true)
{
	if (tid == Operator::Type::RAPTOR_ParCSR) {
		ConstructBlockDiagCSR(comm, glob_size, row_starts, diag);
	} else {
		mfem_error("Raptor BSR not yet implemented");
	}
}


RaptorParMatrix::RaptorParMatrix(const HypreParMatrix *ha, Operator::Type tid)
	: Operator(ha->Height(), ha->Width()),
	  X(NULL), Y(NULL),
	  owns_mat(true)
{
	mat = Convert(*ha, tid);
}



void RaptorParMatrix::Destroy()
{
	if (owns_mat) {
		if (mat) delete mat;
	}
	if (X != NULL) delete X;
	if (Y != NULL) delete Y;
}


RaptorParMatrix::~RaptorParMatrix()
{
	Destroy();
}


void RaptorParMatrix::MakeRef(const RaptorParMatrix & other)
{
	Destroy();

	owns_mat = false;

	height = other.height;
	width = other.width;

	mat = other.mat;
}

raptor::ParMatrix * RaptorParMatrix::Convert(const HypreParMatrix & ha,
                                             Operator::Type tid)
{
	return convert(static_cast<hypre_ParCSRMatrix*>(ha), ha.GetComm());
}


void RaptorParMatrix::CopyCSR(raptor::CSRMatrix & dst,
                              const SparseMatrix & src)
{
	auto & mem_I = src.GetMemoryI();
	auto & mem_J = src.GetMemoryJ();
	auto & mem_data = src.GetMemoryData();

	const int num_rows = src.Height();
	const int nnz = src.NumNonZeroElems();

	host_copy(dst.idx1.data(), mem_I, num_rows + 1);

	dst.idx2.resize(nnz);
	dst.vals.resize(nnz);

	host_copy(dst.idx2.data(), mem_J, nnz);
	host_copy(dst.vals.data(), mem_data, nnz);
}


void RaptorParMatrix::ConstructBlockDiagCSR(MPI_Comm comm,
                                            HYPRE_BigInt glob_size,
                                            HYPRE_BigInt *row_starts,
                                            SparseMatrix *diag)
{
	auto *part = new raptor::Partition(glob_size, glob_size,
	                                   row_starts[1] - row_starts[0],
	                                   row_starts[1] - row_starts[0],
	                                   row_starts[0], row_starts[0]);
	auto *on_proc = new raptor::CSRMatrix(part->local_num_rows, part->local_num_cols, diag->NumNonZeroElems());
	CopyCSR(*on_proc, *diag);
	auto *off_proc = new raptor::CSRMatrix(0, 0, 0);

	mat = new raptor::ParCSRMatrix(part, on_proc, off_proc);
	mat->off_proc->idx1.resize(mat->off_proc->n_rows + 1);
	std::fill(mat->off_proc->idx1.begin(), mat->off_proc->idx1.end(), 0.);
	// mat->off_proc->resize(0, 0);
}


void RaptorParMatrix::Mult(double a, const Vector & x, double b, Vector & y) const
{
   MFEM_ASSERT(x.Size() == Width(), "invalid x.Size() = " << x.Size()
               << ", expected size = " << Width());
   MFEM_ASSERT(y.Size() == Height(), "invalid y.Size() = " << y.Size()
               << ", expected size = " << Height());

   if (X == NULL)
   {
	   X = new RaptorParVector(GetGlobalNumCols(),
	                           Width());
	   Y = new RaptorParVector(GetGlobalNumRows(),
	                           Height());
   }

   X->WrapMemoryRead(x.GetMemory());
   Y->WrapMemoryReadWrite(y.GetMemory());

   // todo: avoid tmp with raptor support for this
   raptor::ParVector tmp(GetGlobalNumRows(), Height());
   raptor::ParCSRMatrix &A = dynamic_cast<raptor::ParCSRMatrix&>(*mat);
   raptor::ParVector *xr = *X;
   raptor::ParVector *yr = *Y;

   A.mult(*xr, tmp);
   if (b != 0.)
	   yr->scale(b);
   else
	   yr->set_const_value(0.);
   yr->axpy(tmp, a);
}


namespace {
Array<int> get_sorted(const Array<int> & rows_cols)
{
	rows_cols.HostRead();
	Array<int> ret;
	ret.SetSize(rows_cols.Size());
	bool sorted = true;
	for (int i = 0; i < rows_cols.Size(); ++i) {
		ret[i] = rows_cols[i];
		if (i && rows_cols[i-1] > rows_cols[i]) sorted = false;
	}
	if (!sorted) ret.Sort();
	return ret;
}


void csr_elim_create(raptor::CSRMatrix &A, raptor::CSRMatrix &Ae, int nrows,
                     int *rows, int ncols, int *cols,
                     int *col_mark)
{
	auto & Ae_i = Ae.idx1;
	auto & A_i = A.idx1;
	auto & A_j = A.idx2;
	int nnz = 0;

	for (int i = 0; i < A.n_rows; ++i) {
		Ae_i[i] = nnz;

		int A_beg = A_i[i];
		int A_end = A_i[i+1];

		if (hypre_BinarySearch(rows, i, nrows) >= 0) {
			// full row
			nnz += A_end - A_beg;

			if (col_mark)
			{
				for (int j = A_beg; j < A_end; ++j) {
					col_mark[A_j[j]] = 1;
				}
			}
		} else {
			// count columns
			for (int j = A_beg; j < A_end; ++j) {
				int col = A_j[j];
				if (hypre_BinarySearch(cols, col, ncols) >= 0) {
					++nnz;
					if (col_mark) {
						col_mark[col] = 1;
					}
				}
			}
		}
	}

	Ae_i[A.n_rows] = nnz;
	Ae.idx2.resize(nnz);
	Ae.vals.resize(nnz);
	Ae.nnz = nnz;
}


/*
  Eliminate rows and columns of A, store eliminated values in Ae.
  If 'diag' is nonzero, the eliminated diagonal of A is set to identity.
  If 'col_remap' is not NULL it specifies renumbering of columns of Ae.
*/
void csr_elim_rowscols(raptor::CSRMatrix & A, raptor::CSRMatrix & Ae,
                       int nrows, int *rows,
                       int ncols, int *cols,
                       bool diag, int *col_remap)
{
	auto & A_i = A.idx1;
	auto & A_j = A.idx2;
	auto & A_data = A.vals;

	auto & Ae_i = Ae.idx1;
	auto & Ae_j = Ae.idx2;
	auto & Ae_data = Ae.vals;

	for (int i = 0; i < A.n_rows; ++i) {
		int A_beg = A_i[i];
		int A_end = A_i[i+1];
		int Ae_beg = Ae_i[i];

		if (hypre_BinarySearch(rows, i, nrows) >= 0) {
			// eliminate row
			for (int j = A_beg, k = Ae_beg; j < A_end; ++j, ++k) {
				int col = A_j[j];
				Ae_j[k] = col_remap ? col_remap[col] : col;
				int a = (diag && col == i) ? 1. : 0.;
				Ae_data[k] = A_data[j] - a;
				A_data[j] = a;
			}
		} else {
			// eliminate columns
			for (int j = A_beg, k = Ae_beg; j < A_end; ++j) {
				int col = A_j[j];
				if (hypre_BinarySearch(cols, col, ncols) >= 0) {
					Ae_j[k] = col_remap ? col_remap[col] : col;
					Ae_data[k] = A_data[j];
					A_data[j] = 0.;
					++k;
				}
			}
		}
	}
}


raptor::ParMatrix * parcsr_eliminateAAe(raptor::ParCSRMatrix & mat,
                                        Array<int> & rowscols_elim)
{
	// todo: not sure why this doesn't work
	// auto *part = mat.partition;
	// ++part->num_shared;
	// auto *Ae = new raptor::ParCSRMatrix(part);

	if (!mat.comm) {
		mat.comm = new raptor::ParComm(mat.partition, mat.off_proc_column_map,
		                               mat.on_proc_column_map);
	}
	auto & comm = dynamic_cast<raptor::CommPkg&>(*mat.comm);

	std::vector<int> elim_diag_col(mat.on_proc_num_cols, 0);

	auto & A_diag = dynamic_cast<raptor::CSRMatrix&>(*mat.on_proc);
	// auto & Ae_diag = dynamic_cast<raptor::CSRMatrix&>(*Ae->on_proc);
	auto Ae_diag_ptr = new raptor::CSRMatrix(mat.local_num_rows,
	                                        mat.on_proc_num_cols);
	auto & Ae_diag = *Ae_diag_ptr;

	for (int i = 0; i < rowscols_elim.Size(); ++i) {
		elim_diag_col[rowscols_elim[i]] = 1;
	}

	comm.init_comm(elim_diag_col);

	csr_elim_create(A_diag, Ae_diag,
	                rowscols_elim.Size(), rowscols_elim.GetData(),
	                rowscols_elim.Size(), rowscols_elim.GetData(),
	                NULL);
	csr_elim_rowscols(A_diag, Ae_diag,
	                  rowscols_elim.Size(), rowscols_elim.GetData(),
	                  rowscols_elim.Size(), rowscols_elim.GetData(),
	                  true, NULL);

	Ae_diag.move_diag();

	auto & elim_offd_col = comm.complete_comm<int>();

	// enumerate offd column indices to eliminate.
	int num_offd_cols_to_elim = 0;
	for (int i = 0; i < mat.off_proc_num_cols; ++i) {
		if (elim_offd_col[i]) ++num_offd_cols_to_elim;
	}

	std::vector<int> offd_cols_to_elim(num_offd_cols_to_elim);
	num_offd_cols_to_elim = 0;
	for (int i = 0; i < mat.off_proc_num_cols; ++i) {
		if (elim_offd_col[i])
			offd_cols_to_elim[num_offd_cols_to_elim++] = i;
	}

	std::vector<int> col_mark(mat.off_proc_num_cols);
	std::vector<int> col_remap(mat.off_proc_num_cols);

	auto & A_offd = dynamic_cast<raptor::CSRMatrix&>(*mat.off_proc);
	// auto & Ae_offd = dynamic_cast<raptor::CSRMatrix&>(*Ae->on_proc);
	auto Ae_offd_ptr = new raptor::CSRMatrix(mat.local_num_rows,
	                                        mat.off_proc_num_cols);
	auto & Ae_offd = *Ae_offd_ptr;

	csr_elim_create(A_offd, Ae_offd,
	                rowscols_elim.Size(), rowscols_elim.GetData(),
	                num_offd_cols_to_elim, offd_cols_to_elim.data(),
	                col_mark.data());
	for (int i = 0, k = 0; i < mat.off_proc_num_cols; ++i) {
		if (col_mark[i]) col_remap[i] = k++;
	}

	csr_elim_rowscols(A_offd, Ae_offd,
	                  rowscols_elim.Size(), rowscols_elim.GetData(),
	                  num_offd_cols_to_elim, offd_cols_to_elim.data(),
	                  false, col_remap.data());

	int Ae_offd_ncols = 0;
	for (int i = 0; i < mat.off_proc_num_cols; ++i) {
		if (col_mark[i]) ++Ae_offd_ncols;
	}

	auto *part = mat.partition;
	// TODO: can't use this ctor since it calls ParComm ctor without fixing indices
	// auto Ae = new raptor::ParCSRMatrix(part, Ae_diag_ptr, Ae_offd_ptr);
	auto Ae = new raptor::ParCSRMatrix(part, mat.global_num_rows, mat.global_num_cols,
	                                   mat.local_num_rows, Ae_diag.n_cols, Ae_offd_ncols,
	                                   0, false);
	Ae->on_proc = Ae_diag_ptr;
	Ae->off_proc = Ae_offd_ptr;

	auto & Ae_colmap_offd = Ae->get_off_proc_column_map();
	auto & A_colmap_offd = mat.get_off_proc_column_map();
	// Ae_colmap_offd.resize(Ae_offd_ncols);
	Ae->off_proc_num_cols = Ae_offd_ncols;

	Ae->finalize(false);

	Ae_offd_ncols = 0;
	for (int i = 0; i < mat.off_proc_num_cols; ++i) {
		if (col_mark[i])
			Ae_colmap_offd[Ae_offd_ncols++] = A_colmap_offd[i];
	}

	Ae->comm = new raptor::ParComm(Ae->partition, Ae->off_proc_column_map,
	                               Ae->on_proc_column_map);

	return Ae;
}
}


RaptorParMatrix * RaptorParMatrix::EliminateRowsCols(const Array<int> & rows_cols)
{
	auto rc_sorted = get_sorted(rows_cols);

	MFEM_VERIFY(mat, "no associated Raptor matrix object");
	auto *csrmat = dynamic_cast<raptor::ParCSRMatrix*>(mat);
	MFEM_VERIFY(csrmat, "EliminateRowsCols only supported for ParCSR");
	auto rmat = parcsr_eliminateAAe(*csrmat, rc_sorted);

	return new RaptorParMatrix(rmat, true);
}


Operator::Type RaptorParMatrix::GetType() const
{
	MFEM_VERIFY(mat, "no associated Raptor matrix object");

	if (dynamic_cast<raptor::ParBSRMatrix*>(mat))
		return Operator::Type::RAPTOR_ParBSR;
	else {
		MFEM_VERIFY(dynamic_cast<raptor::ParCSRMatrix*>(mat),
		            "Raptor GetType() invalid matrix type");
		return Operator::Type::RAPTOR_ParCSR;
	}
}

RaptorParMatrix::operator raptor::ParCSRMatrix * () const
{
	MFEM_VERIFY(mat, "no asspociated Raptor matrix object");
	MFEM_VERIFY(GetType() == Operator::Type::RAPTOR_ParCSR,
	            "Invalid ParCSRMatrix type for conversion to raptor::ParCSRMatrix*");
	return dynamic_cast<raptor::ParCSRMatrix*>(mat);
}


RaptorParMatrix::operator raptor::ParBSRMatrix * () const
{
	MFEM_VERIFY(mat, "no asspociated Raptor matrix object");
	MFEM_VERIFY(GetType() == Operator::Type::RAPTOR_ParBSR,
	            "Invalid ParCSRMatrix type for conversion to raptor::ParCSRMatrix*");
	return dynamic_cast<raptor::ParBSRMatrix*>(mat);
}


void RaptorParMatrix::Print(const char *fname) const
{
	auto & diag = *mat->on_proc;
	auto & offd = *mat->off_proc;

	std::ofstream ofile(fname);
	ofile.precision(14);
	for (std::size_t i = 0; i < diag.n_rows; ++i) {
		auto write_row = [&](const raptor::Matrix & m, const std::vector<int> & colmap) {
			const auto & rowptr = m.idx1;
			const auto & colind = m.idx2;
			const auto & values = m.vals;
			for (std::size_t off = rowptr[i]; off < rowptr[i+1]; ++off) {
				ofile << mat->get_local_row_map()[i] << " "
				      << colmap[colind[off]] << " "
				      << std::scientific << values[off] << '\n';
			}
		};

		write_row(diag, mat->get_on_proc_column_map());
		write_row(offd, mat->get_off_proc_column_map());
	}
}


RaptorSolver::RaptorSolver() : A(NULL), B(NULL), X(NULL) {}
RaptorSolver::RaptorSolver(const RaptorParMatrix * a) :
	Solver(a->Height(), a->Width()), A(a), B(NULL), X(NULL) {}

void RaptorSolver::Mult(const Vector & b, Vector & x) const
{
	WrapVectors(b, x);
	Mult(*B, *X);
}


void RaptorSolver::WrapVectors(const Vector &b, Vector &x) const
{
	if (B == NULL) {
		B = new RaptorParVector(A->GetGlobalNumRows(),
		                        A->Height());
		X = new RaptorParVector(A->GetGlobalNumCols(),
		                        A->Width());
	}

	B->WrapMemoryRead(b.GetMemory());
	X->WrapMemoryReadWrite(x.GetMemory());
}


void RaptorSolver::Mult(const RaptorParVector & b, RaptorParVector & x) const
{
	if (A == NULL) {
		mfem_error("RaptorSolver::Mult(...) : RaptorParMatrix A is missing");
		return;
	}

	raptor::ParVector *xr = x;
	const raptor::ParVector *br = b;

	if (!iterative_mode) {
		xr->set_const_value(0);
	}

	Setup();
	Solve(*br, *xr);
}


RaptorSolver::~RaptorSolver()
{
	if (B) delete B;
	if (X) delete X;
}


RaptorMultilevel::RaptorMultilevel(raptor::ParMultilevel * t) : RaptorSolver(),
                                                                setup_called(false),
                                                                target(t), iters(0) {}


RaptorMultilevel::RaptorMultilevel(raptor::ParMultilevel *t,
                                   const RaptorParMatrix *a)
	: RaptorSolver(a),
      setup_called(false), target(t), iters(0) {}


void RaptorMultilevel::Setup() const
{
	if (setup_called) return;

	raptor::ParCSRMatrix *Ar = *A;
	target->setup(Ar);

	setup_called = true;
}


void RaptorMultilevel::Solve(const raptor::ParVector & b, raptor::ParVector & x) const
{
	iters = target->solve(x, const_cast<raptor::ParVector&>(b));
}


void RaptorMultilevel::PrintResiduals() const
{
	target->print_residuals(iters);
}


void RaptorMultilevel::SetOperator(const Operator &op)
{
	const RaptorParMatrix *new_A = dynamic_cast<const RaptorParMatrix*>(&op);
	MFEM_VERIFY(new_A, "new Operator must be a RaptorParMatrix");

	height = new_A->Height();
	width  = new_A->Width();
	A = const_cast<RaptorParMatrix*>(new_A);
	setup_called = 0;
	delete X;
	delete B;
	X = B = NULL;
}

RaptorMultilevel::~RaptorMultilevel()
{
	delete target;
}


RaptorRugeStuben::RaptorRugeStuben()
	: RaptorMultilevel(new raptor::ParRugeStubenSolver()) {
	SetDefaultOptions();
}

RaptorRugeStuben::RaptorRugeStuben(const RaptorParMatrix &A)
	: RaptorMultilevel(new raptor::ParRugeStubenSolver(), &A) {
	SetDefaultOptions();
}

void RaptorRugeStuben::SetDefaultOptions()
{
	target->max_iterations = 1; // default for use as preconditioner
	target->relax_type = raptor::relax_t::Jacobi;
	target->strength_type = raptor::strength_t::Symmetric;
}


RaptorPCG::RaptorPCG(MPI_Comm) : maxiter_(1000), rtol_(1e-8), precond_(NULL) {
	iterative_mode = true;
}


void RaptorPCG::SetOperator(const Operator & op)
{
	const RaptorParMatrix *new_A = dynamic_cast<const RaptorParMatrix*>(&op);
	MFEM_VERIFY(new_A, "new Operator must be a RaptorParMatrix!");

   height = new_A->Height();
   width  = new_A->Width();

   A = const_cast<RaptorParMatrix*>(new_A);
   if (precond_) {
	   precond_->SetOperator(*A);
	   SetPreconditioner(*precond_);
   }

   delete X;
   delete B;
   B = X = NULL;
}


void RaptorPCG::SetPreconditioner(RaptorMultilevel & ml)
{
	precond_ = &ml;
}


void RaptorPCG::Solve(const raptor::ParVector & b, raptor::ParVector & x) const
{
	raptor::ParCSRMatrix *parcsr_A = *A;
	raptor::ParMultilevel *ml = *precond_;

	precond_->Setup();

	std::vector<double> res;
	raptor::PCG(parcsr_A, ml, x, const_cast<raptor::ParVector&>(b), res, rtol_, maxiter_);
}



RaptorParMatrix * RAP(RaptorParMatrix * A, RaptorParMatrix * P)
{
	MFEM_VERIFY(A->Width() == P->Height(),
	            "Raptor RAP: Number of local cols of A " << A->Width() <<
	            " differs from the number of local rows of P " << P->Height());
	MFEM_VERIFY(A->GetType() == P->GetType(), "Raptor RAP: mixing operator types is not supported");

	if (A->GetType() == Operator::Type::RAPTOR_ParCSR) {
		raptor::ParCSRMatrix *Ar = *A;
		raptor::ParCSRMatrix *Pr = *P;
		raptor::ParCSRMatrix *AP = Ar->mult(Pr);
		raptor::ParCSRMatrix *PtAP = AP->mult_T(Pr);

		delete AP;

		return new RaptorParMatrix(PtAP);
	} else {
		mfem_error("Raptor RAP: BSR not implemented");
		return nullptr;
	}
}


void EliminateBC(RaptorParMatrix & A, RaptorParMatrix & Ae,
                 const Array<int> & ess_dof_list,
                 const Vector &x, Vector & b)
{
	// b -= Ae * x
	Ae.Mult(-1., x, 1., b);

	if (ess_dof_list.Size() == 0) return;

	x.HostRead();
	b.HostReadWrite();

	auto & A_diag = *((raptor::ParCSRMatrix*)A)->on_proc;
	auto & data = A_diag.vals;
	auto & I = A_diag.idx1;

	for (int i = 0; i < ess_dof_list.Size(); ++i) {
		int r = ess_dof_list[i];
		b(r) = data[I[r]] * x(r);
	}
}

}

#endif // MFEM_USE_MPI
#endif // MFEM_USE_RAPTOR
