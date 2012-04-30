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
  printf("\nWang+Chen potentials\n");

  double rb = 0.98, r1 = 0.68;
  ham.set_impurity(new exp_wang(1000, 1.06, rb, r1));
  printf("Be_Ga: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, Be));
  printf("Be_Ga: %f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, 1000);
  ham.set_impurity_parameter(_ra, 1.06);
  ham.set_impurity_parameter(_rb, rb);
  ham.set_impurity_parameter(_r1, r1);
  printf("Be_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(1300, 1.40, rb, r1));
  printf("Mg_Ga: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, Mg));
  printf("Mg_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(1450, 1.74, rb, r1));
  printf("Ca_Ga: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, Ca));
  printf("Ca_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(3400, 1.31, rb, r1));
  printf("Zn_Ga: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, Zn));
  printf("Zn_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(4050, 1.48, rb, r1));
  printf("Cd_Ga: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, Cd));
  printf("Cd_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(-8050, 0.77, rb, r1));
  printf("C_N: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, C));
  printf("C_N: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(2050, 1.17, rb, r1));
  printf("Si_N: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, Si));
  printf("Si_N: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(2950, 1.22, rb, r1));
  printf("Ge_N: %f\n", ham.get_eval(0));
  ham.set_impurity(new exp_wang(GaN, Ge));
  printf("Ge_N: %f\n", ham.get_eval(0));

  /*  printf("AlN: Energy vs crystal field splitting\n");
  ham.set_crystal_parameter(_d1,105);
  for(double i = ham.get_crystal_parameter(_d1); i < 386;
      i = ham.get_crystal_parameter(_d1))
    {
      printf("%.1f: %f\n", i, ham.get_eval(0));
      i = ham.get_crystal_parameter(_d1);
      ham.set_crystal_parameter(_d1, i+5);
    }

  ham.set_crystal(new exp_wz(AlN));
  ham.set_impurity(new exp_LCZ(Al, Mg));
  printf("Mg_Al in AlN: %f\n", ham.get_eval(0));

  ham.set_crystal(new exp_wz(GaN));
  ham.set_impurity(new exp_LCZ(Ga, Mg));
  printf("Mg_Ga in GaN: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga, Zn));
  printf("Zn_Ga in GaN (3d as core): %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga_val, Zn));
  printf("Zn_Ga in GaN (3d as valance): %f\n", ham.get_eval(0));

  ham.set_crystal(new exp_wz(InN));
  ham.set_impurity(new exp_LCZ(Ga, Zn));
  printf("Zn_Ga in InN (3d as core): %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga_val, Zn));
  printf("Zn_Ga in InN (3d as valance): %f\n", ham.get_eval(0));*/

  return 0;
}
