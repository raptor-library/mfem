#ifndef EXACT
#define EXACT

#include "mfem.hpp"
#include "plasma_model.hpp"
using namespace mfem;
using namespace std;

class ExactCoefficient : public Coefficient
{
private:
  double r0;
  double z0;
  double k;
  bool use_manufactured;
public:
  ExactCoefficient(double r0_, double z0_, double k_, bool use_manufactured_) : r0(r0_), z0(z0_), k(k_), use_manufactured(use_manufactured_) { }
  ExactCoefficient() : r0(0.0), z0(0.0), k(1.0) { }
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  virtual ~ExactCoefficient() { }
};

class ExactForcingCoefficient : public Coefficient
{
private:
  double r0;
  double z0;
  double k;
  PlasmaModelBase model;
  Vector *coil_current;
  bool use_manufactured;
public:
  ExactForcingCoefficient(double r0_, double z0_, double k_, PlasmaModelBase model_, bool use_manufactured_) : r0(r0_), z0(z0_), k(k_), model(model_), use_manufactured(use_manufactured_) { }
  void set_coil_current(Vector *coil_current_) {
    coil_current = coil_current_;
  }
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  virtual ~ExactForcingCoefficient() { }
};

double psi_exact(double r, double z, double r0, double z0, double k);
#endif