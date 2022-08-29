#include "mfem.hpp"
#include "sys_operator.hpp"
using namespace mfem;
using namespace std;


double GetMaxError(LinearForm &res) {
  double *resv = res.GetData();
  int size = res.FESpace()->GetTrueVSize();

  double max_val = - numeric_limits<double>::infinity();
  for (int i = 0; i < size; ++i) {
    if (abs(resv[i]) > max_val) {
      max_val = abs(resv[i]);
    }
  }
  return max_val;
}

void SysOperator::Mult(const Vector &psi, Vector &y) const {
  // diff_operator * psi - plasma_term(psi) * psi - coil_term

  GridFunction x(fespace);
  x = psi;
  int ind_min, ind_max;
  double min_val, max_val;
  int iprint = 1;
  compute_plasma_points(x, *mesh, vertex_map, ind_min, ind_max, min_val, max_val, iprint);
  NonlinearGridCoefficient nlgcoeff1(model, 1, &x, min_val, max_val, attr_lim);
  LinearForm plasma_term(fespace);
  plasma_term.AddDomainIntegrator(new DomainLFIntegrator(nlgcoeff1));
  plasma_term.Assemble();
  
  y = *coil_term;
  add(y, -1.0, plasma_term, y);
  diff_operator->AddMult(psi, y);
}

Operator &SysOperator::GetGradient(const Vector &psi) const {
  // diff_operator - sum_{i=1}^3 diff_plasma_term_i(psi)

  delete Mat;
  GridFunction x(fespace);
  x = psi;

  int ind_min, ind_max;
  double min_val, max_val;
  int iprint = 0;
  compute_plasma_points(x, *mesh, vertex_map, ind_min, ind_max, min_val, max_val, iprint);
 
  NonlinearGridCoefficient nlgcoeff_2(model, 2, &x, max_val, min_val, attr_lim);
  BilinearForm diff_plasma_term_2(fespace);
  diff_plasma_term_2.AddDomainIntegrator(new MassIntegrator(nlgcoeff_2));
  diff_plasma_term_2.Assemble();

  NonlinearGridCoefficient nlgcoeff_3(model, 3, &x, max_val, min_val, attr_lim);
  LinearForm diff_plasma_term_3(fespace);
  diff_plasma_term_3.AddDomainIntegrator(new DomainLFIntegrator(nlgcoeff_3));
  diff_plasma_term_3.Assemble();

  NonlinearGridCoefficient nlgcoeff_4(model, 4, &x, max_val, min_val, attr_lim);
  LinearForm diff_plasma_term_4(fespace);
  diff_plasma_term_4.AddDomainIntegrator(new DomainLFIntegrator(nlgcoeff_4));
  diff_plasma_term_4.Assemble();

  const SparseMatrix M1 = diff_operator->SpMat();
  SparseMatrix M2 = diff_plasma_term_2.SpMat();
  M2.Finalize();
  
  int m = fespace->GetTrueVSize();
  Mat = new SparseMatrix(m, m);
  for (int k = 0; k < m; ++k) {
    Mat->Add(k, ind_max, diff_plasma_term_3[k]);
    Mat->Add(k, ind_min, diff_plasma_term_4[k]);
  }

  // diff operator
  int height;
  const auto II1 = M1.ReadI();
  const auto JJ1 = M1.ReadJ();
  const auto AA1 = M1.ReadData();
  height = M1.Height();
  for (int i = 0; i < height; ++i) {
    const int begin1 = II1[i];
    const int end1 = II1[i+1];
    
    int j;
    for (j = begin1; j < end1; j++) {
      Mat->Add(i, JJ1[j], AA1[j]);
    }
  }
  const auto II2 = M2.ReadI();
  const auto JJ2 = M2.ReadJ();
  const auto AA2 = M2.ReadData();
  
  height = M2.Height();
  for (int i = 0; i < height; ++i) {
    const int begin2 = II2[i];
    const int end2 = II2[i+1];

    int j;
    for (j = begin2; j < end2; j++) {
      Mat->Add(i, JJ2[j], AA2[j]);
    }
  }  
  Mat->Finalize();
  
  return *Mat;
}
