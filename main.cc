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
  printf("GaN: Energy vs crystal field splitting\n");
  for(double i = ham.get_crystal_parameter(_d1); i > -37;
      i = ham.get_crystal_parameter(_d1))
    {
      printf("%.1f: %f\n", i, ham.get_eval(0));
      i = ham.get_crystal_parameter(_d1);
      ham.set_crystal_parameter(_d1, i-1);
    }

  ham.set_crystal(new exp_wz(InN));
  ham.set_crystal_parameter(_d1, ham.get_crystal_parameter(_d1)+20);
  printf("InN: Energy vs crystal field splitting\n");
  for(double i = ham.get_crystal_parameter(_d1); i > -60;
      i = ham.get_crystal_parameter(_d1))
    {
      printf("%.1f: %f\n", i, ham.get_eval(0));
      i = ham.get_crystal_parameter(_d1);
      ham.set_crystal_parameter(_d1, i-1);
    }

  ham.set_crystal(new exp_wz(AlN));
  printf("AlN BE: %g\n", ham.get_eval(0));
  ham.set_crystal_parameter(_A6,1.05);
  ham.set_crystal_parameter(_d1,245);
  ham.set_crystal_parameter(_d2,-6.3);
  ham.set_crystal_parameter(_d3,-6.3);
  printf("AlN quasicubic BE: %g\n", ham.get_eval(0));

  ham.set_crystal(new exp_wz(GaN));
  printf("GaN BE: %g\n", ham.get_eval(0));
  ham.set_crystal_parameter(_A6,3.31);
  ham.set_crystal_parameter(_d1,-18.2);
  ham.set_crystal_parameter(_d2,-1.8);
  ham.set_crystal_parameter(_d3,-1.8);
  printf("GaN quasicubic BE: %g\n", ham.get_eval(0));

  ham.set_crystal(new exp_wz(InN));
  printf("InN BE: %g\n", ham.get_eval(0));
  ham.set_crystal_parameter(_A6,9.45);
  ham.set_crystal_parameter(_d1,-43.4);
  ham.set_crystal_parameter(_d2,3.07);
  ham.set_crystal_parameter(_d3,3.07);
  printf("InN quasicubic BE: %g\n", ham.get_eval(0));

  /*
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
  printf("C_N: %g\n", ham.get_eval(0));*/
  return 0;
}
