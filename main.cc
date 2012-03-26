#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

int main(int argc, char **argv)
{
  hamiltonian ham
    (
     new exp_wz(GaN),
     new exp_coulomb(),
     new exp_overlap()
    );
  printf("GaN: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ga,Be));
  printf("Be_Ga: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ga,Mg));
  printf("Mg_Ga: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ga,Zn));
  printf("Zn_Ga: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ga,Ca));
  printf("Ca_Ga: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(N,C));
  printf("C_N: %g\n", ham.get_eval(0));
  ham.set_crystal(new exp_gwz(ZnGeN2));
  ham.set_impurity(new exp_coulomb());
  printf("ZnGeN2: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Zn,Li));
  printf("Li_Zn: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ge,B));
  printf("B_Ge: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Zn,Na));
  printf("Na_Zn: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ge,Al));
  printf("Al_Ge: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Zn,Cu));
  printf("Cu_Zn: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ge,Ga));
  printf("Ga_Ge: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Zn,K));
  printf("K_Zn: %g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(N,C));
  printf("C_N: %g\n", ham.get_eval(0));
  return 0;
}
