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
  /*
  printf("Mireles potentials\n");
  printf("Coulomb_GaN: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga, Be));
  printf("Be_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity_parameter(_dielectric_ratio, 0.5);
  printf("half dielectric: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga, Mg));
  printf("Mg_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity_parameter(_dielectric_ratio, 0.5);
  printf("half dielectric: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga, Ca));
  printf("Ca_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity_parameter(_dielectric_ratio, 0.5);
  printf("half dielectric: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga, Zn));
  printf("Zn_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity_parameter(_dielectric_ratio, 0.5);
  printf("half dielectric: %f\n", ham.get_eval(0));

  ham.set_crystal(new exp_wz(AlN));
  ham.set_impurity(new exp_coulomb());
  printf("Coulomb_AlN: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Al, Be));
  printf("Be_Al: %f\n", ham.get_eval(0));

  ham.set_impurity_parameter(_dielectric_ratio, 0.5);
  printf("half dielectric: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Al, Mg));
  printf("Mg_Al: %f\n", ham.get_eval(0));

  ham.set_impurity_parameter(_dielectric_ratio, 0.5);
  printf("half dielectric: %f\n", ham.get_eval(0));
  */
#ifdef off
  printf("\nWang+Chen potentials\n");
  printf("new params\n");
  printf("Coulomb_GaN: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(1000, 1.06, 0.98, 0.68));
  printf("Be_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(1300, 1.40, 0.98, 0.68));
  printf("Mg_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(1450, 1.74, 0.98, 0.68));
  printf("Ca_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(3400, 1.31, 0.98, 0.68));
  printf("Zn_Ga: %f\n", ham.get_eval(0));
#endif
  printf("old params\n");
  ham.set_crystal(new exp_zb(GaN));
  double A = 8.12, B = 2.24, C = 9.96, eta = 1.12;
  ham.set_crystal_parameter(_g1, (A+2*B)/3*eta);
  ham.set_crystal_parameter(_g2, (A-B)/6*eta);
  ham.set_crystal_parameter(_g3, C/6*eta);
  ham.set_crystal_parameter(_d0, 80);
  ham.set_crystal_parameter(_dielectric, 11.02);
  printf("Coulomb_GaN: %f\n", ham.get_eval(0));

  double rb = 1.18, r1 = 0.82;
  ham.set_impurity(new exp_wang(1000, 1.06, rb, r1));
  printf("Be_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(1300, 1.40, rb, r1));
  printf("Mg_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(1450, 1.74, rb, r1));
  printf("Ca_Ga: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(3400, 1.31, rb, r1));
  printf("Zn_Ga: %f\n", ham.get_eval(0));

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
