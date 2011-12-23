#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

int main(int argc, char **argv)
{
  gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
  printf("%d\n", old_handler);
  hamiltonian ham
    (
     new exp_zb(1.0, 0.0, 0.0, 0.0, 1.0),
     new exp_coulomb(),
     new exp_overlap()
    );
  for(int i = 0; i < 5; i++)
    printf("%g\n", ham.get_eval(i));
  return 0;
}
