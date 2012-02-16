#include <math.h>
#include <stdio.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

#include "matrix_term.hh"
#include "defs.hh"

// defined in matrix_term.cc
extern int matrix_term_block_size;

// ############################## exp_zb #############################
exp_zb::exp_zb
(double _g1, double _g2, double _g3, double _d0, double _dielectric)
{
  g1 = _g1;
  g2 = _g2;
  g3 = _g3;
  d0 = _d0 * RYD_PER_MEV;
  dielectric = _dielectric;
  inv_radius = 1.0 / _g1 / _dielectric;
}

exp_zb::exp_zb(crystals_t c)
{
  switch(c)
    {
    case GaN:
      g1 = 2.46;
      g2 = 0.65;
      g3 = 0.98;
      d0 = 19 * RYD_PER_MEV;
      dielectric = 10.13;
      break;
    case AlN:
      g1 = 1.40;
      g2 = 0.35;
      g3 = 0.59;
      d0 = 19 * RYD_PER_MEV;
      dielectric = 9.67;
      break;
    default:
      printf("Error: exp_zb instantiated with unknown crystal");
      exit(-1);
    }
  inv_radius = 1.0 / g1 / dielectric;

}

gsl_matrix_complex *exp_zb::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_zb_def.hh"
  return output;
}


// ############################## exp_wz #############################
exp_wz::exp_wz(double _A1, double _A2, double _A3, double _A4, double _A5,
	       double _A6, double _d1, double _d2, double _d3,
	       double _dielectric)
{
  A1 = _A1;
  A2 = _A2;
  A3 = _A3;
  A4 = _A4;
  A5 = _A5;
  A6 = _A6;
  d1 = _d1 * RYD_PER_MEV;
  d2 = _d2 * RYD_PER_MEV;
  d3 = _d3 * RYD_PER_MEV;
  dielectric = _dielectric;
  inv_radius = 1.0 / (_A2 + _A4) / _dielectric;
}

exp_wz::exp_wz(crystals_t c)
{
  switch(c)
    {
    case GaN:
      A1 = 5.98;
      A2 = 0.58;
      A3 = -5.44;
      A4 = 2.46;
      A5 = 2.53;
      A6 = 1.55;
      d1 = -12 * RYD_PER_MEV;
      d2 = -3.9 * RYD_PER_MEV;
      d3 = -5.4 * RYD_PER_MEV;
      dielectric = 10.13;
      break;
    case AlN:
      A1 = 4.05;
      A2 = 0.28;
      A3 = -3.71;
      A4 = 1.71;
      A5 = 1.90;
      A6 = 1.05;
      d1 = 245 * RYD_PER_MEV;
      d2 = -6.2 * RYD_PER_MEV;
      d3 = -7.5 * RYD_PER_MEV;
      dielectric = 9.67;
      break;
    case InN:
      A1 = 15.72;
      A2 = 0.7;
      A3 = -15.23;
      A4 = 5.79;
      A5 = 5.9;
      A6 = 4.58;
      d1 = 14.8 * RYD_PER_MEV;
      d2 = 3.4 * RYD_PER_MEV;
      d3 = -2.2 * RYD_PER_MEV;
      dielectric = 14;
      break;
    default:
      printf("Error: exp_wz instantiated with unknown crystal");
      exit(-1);
    }
  inv_radius = 1.0 / (A2 + A4) / dielectric;
}


gsl_matrix_complex *exp_wz::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_wz_def.hh"
  return output;
}


// ############################# exp_gwz #############################
exp_gwz::exp_gwz(double _A1, double _A2, double _A3, double _B1, double _B2,
		 double _B3, double _C1, double _C2, double _C3, double _D1,
		 double _D2, double _D3, double _d1c, double _d2c,
		 double _d1so, double _d2so, double _d3so, double _dielectric)
{
  A1 = _A1;
  A2 = _A2;
  A3 = _A3;
  B1 = _B1;
  B2 = _B2;
  B3 = _B3;
  C1 = _C1;
  C2 = _C2;
  C3 = _C3;
  D1 = _D1;
  D2 = _D2;
  D3 = _D3;
  d1c = _d1c * RYD_PER_MEV;
  d2c = _d2c * RYD_PER_MEV;
  d1so = _d1so * RYD_PER_MEV;
  d2so = _d2so * RYD_PER_MEV;
  d3so = _d3so * RYD_PER_MEV;
  dielectric = _dielectric;
  inv_radius = 1.0 / (_B1 + _B2) / _dielectric;
}

exp_gwz::exp_gwz(crystals_t c)
{
  switch(c)
    {
    case ZnGeN2:
      A1 = 6.82;
      A2 = -6.44;
      A3 = -0.01;
      B1 = 0.51;
      B2 = 2.19;
      B3 = -0.09;
      C1 = 0.02;
      C2 = 0.05;
      C3 = -2.30;
      D1 = -4.6;
      D2 = -2.76;
      D3 = -2.76;
      d1c = -115 * RYD_PER_MEV;
      d2c = 14 * RYD_PER_MEV;
      d1so = 0 * RYD_PER_MEV;
      d2so = 0 * RYD_PER_MEV;
      d3so = 0 * RYD_PER_MEV;
      dielectric = 9.7;
      break;
    case ZnSnN2:
      A1 = 8.57;
      A2 = -8.10;
      A3 = -0.02;
      B1 = 0.53;
      B2 = 3.12;
      B3 = -0.11;
      C1 = 0.03;
      C2 = 0.05;
      C3 = -3.18;
      D1 = -6.36;
      D2 = -4.62;
      D3 = -4.62;
      d1c = -82 * RYD_PER_MEV;
      d2c = 94 * RYD_PER_MEV;
      d1so = 0 * RYD_PER_MEV;
      d2so = 0 * RYD_PER_MEV;
      d3so = 0 * RYD_PER_MEV;
      dielectric = 12.71;
      break;
    default:
      printf("Error: exp_gwz instantiated with unknown crystal");
      exit(-1);
    }
  inv_radius = 1.0 / (B1 + B2) / dielectric;
}

gsl_matrix_complex *exp_gwz::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_gwz_def.hh"
  return output;
}


// ############################ exp_coulomb ##########################
exp_coulomb::exp_coulomb()
{
  dielectric = inv_radius = 1.0;
}

gsl_matrix_complex *exp_coulomb::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_coulomb_def.hh"
  return output;
}


// ############################# exp_wang ############################
exp_wang::exp_wang(double _V, double _ra, double _rb, double _r1)
{
  V = _V;
  ra = _ra;
  rb = _rb;
  r1 = _r1;
  dielectric = inv_radius = 1.0;
}

gsl_matrix_complex *exp_wang::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_well_def.hh"
  gsl_matrix_complex *out2 = output;
  output = gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_dielectric_def.hh"
  gsl_matrix_complex_add(output, out2);
  gsl_matrix_complex_free(out2);
  return output;
}


// ############################# exp_LCZ #############################
exp_LCZ::exp_LCZ(elements_t _host, elements_t _impurity)
{
  host = new exp_LCZ_atom(_host);
  impurity = new exp_LCZ_atom(_impurity);
  coulomb = new exp_coulomb();
}

gsl_matrix_complex *exp_LCZ::matrix(double min, double max, size_t num)
{
  gsl_matrix_complex *h = host->matrix(min, max, num),
    *i = impurity->matrix(min, max, num),
    *c = coulomb->matrix(min, max, num);
  gsl_matrix_complex_sub(h,i);
  gsl_matrix_complex_add(h,c);
  gsl_matrix_complex_free(i);
  gsl_matrix_complex_free(c);
  return h;
}

// no-op because blocks are calculated in the lower-level terms
gsl_matrix_complex *exp_LCZ::matrix_block(double a1, double a2)
{
  return NULL;
}

void exp_LCZ::on_set_inv_radius(double r)
{
  host->set_inv_radius(r);
  impurity->set_inv_radius(r);
  coulomb->set_inv_radius(r);
}

void exp_LCZ::on_set_dielectric_constant(double k)
{
  host->set_dielectric_constant(k);
  impurity->set_dielectric_constant(k);
  coulomb->set_dielectric_constant(k);
}

void exp_LCZ::on_delete()
{
  delete host;
  delete impurity;
  delete coulomb;
}


// ########################## exp_LCZ_atom ###########################
exp_LCZ_atom::exp_LCZ_atom(elements_t atom)
{
  set_d_core(atom);
}

void exp_LCZ_atom::set_d_core(elements_t atom)
{
  switch(atom)
    {
    case Li:
      Zc = 2;
      C1 = 2.9990;
      C2 = 1.0472;
      C3 = 1.9686;
      break;
    case Be:
      Zc = 2;
      C1 = 3.0169;
      C2 = 1.4229;
      C3 = 2.9456;
      break;
    case B:
      Zc = 2;
      C1 = 3.0315;
      C2 = 1.8014;
      C3 = 3.9668;
      break;
    case C:
      Zc = 2;
      C1 = 3.0731;
      C2 = 2.1922;
      C3 = 5.0062;
      break;
    case N:
      Zc = 2;
      C1 = 3.0487;
      C2 = 2.5638;
      C3 = 6.0734;
      break;
    case O:
      Zc = 2;
      C1 = 3.0500;
      C2 = 2.9600;
      C3 = 6.8700;
      break;
    case F:
      Zc = 2;
      C1 = 3.0744;
      C2 = 3.3609;
      C3 = 7.7131;
      break;
    case Ne:
      Zc = 2;
      C1 = 3.0765;
      C2 =  3.7456;
      C3 = 8.7366;
      break;
    case Na:
      Zc = 10;
      C1 = 8.6691;
      C2 = 1.354;
      C3 = 2.3326;
      break;
    case Mg:
      Zc = 10;
      C1 = 8.0022;
      C2 = 1.4094;
      C3 = 2.6122;
      break;
    case Al:
      Zc = 10;
      C1 = 8.0223;
      C2 = 1.5274;
      C3 = 3.0200;
      break;
    case Si:
      Zc = 10;
      C1 = 8.0929;
      C2 = 1.6490;
      C3 = 3.4389;
      break;
    case P:
      Zc = 10;
      C1 = 8.3029;
      C2 = 1.8115;
      C3 = 3.8615;
      break;
    case S:
      Zc = 10;
      C1 = 8.4635;
      C2 = 1.9634;
      C3 = 4.3209;
      break;
    case Cl:
      Zc = 10;
      C1 = 8.6365;
      C2 = 2.1205;
      C3 = 4.8235;
      break;
    case Ar:
      Zc = 10;
      C1 = 8.8225;
      C2 = 2.2825;
      C3 = 5.4118;
      break;
    case K:
      Zc = 18;
      C1 = 10.615;
      C2 = 0.8133;
      C3 = 1.5698;
      break;
    case Ca:
      Zc = 18;
      C1 = 11.154;
      C2 = 0.8898;
      C3 = 1.8822;
      break;
    case Sc:
      Zc = 18;
      C1 = 11.706;
      C2 = 0.9809;
      C3 = 2.1414;
      break;
    case Ti:
      Zc = 18;
      C1 = 12.210;
      C2 = 1.0726;
      C3 = 2.3685;
      break;
    case V:
      Zc = 18;
      C1 = 12.655;
      C2 = 1.1619;
      C3 = 2.5793;
      break;
    case Cr:
      Zc = 18;
      C1 = 12.781;
      C2 = 1.2449;
      C3 = 2.7634;
      break;
    case Mn:
      Zc = 18;
      C1 = 13.257;
      C2 = 1.3321;
      C3 = 2.9744;
      break;
    case Fe:
      Zc = 18;
      C1 = 13.661;
      C2 = 1.4247;
      C3 = 3.1611;
      break;
    case Co:
      Zc = 18;
      C1 = 13.761;
      C2 = 1.4997;
      C3 = 3.3435;
      break;
    case Ni:
      Zc = 18;
      C1 = 14.134;
      C2 = 1.5925;
      C3 = 3.5240;
      break;
    case Cu:
      Zc = 18;
      C1 = 14.230;
      C2 = 1.6794;
      C3 = 3.6701;
      break;
    case Zn:
      Zc = 18;
      C1 = 14.557;
      C2 = 1.7601;
      C3 = 3.8772;
      break;
    case Ga:
      Zc = 28;
      C1 = 14.934;
      C2 = 1.8452;
      C3 = 3.3476;
      break;
    case Ge:
      Zc = 28;
      C1 = 14.888;
      C2 = 1.8923;
      C3 = 3.5107;
      break;
    case As:
      Zc = 28;
      C1 = 14.954;
      C2 = 1.9544;
      C3 = 3.6166;
      break;
    case Se:
      Zc = 28;
      C1 = 15.171;
      C2 = 2.0278;
      C3 = 3.8322;
      break;
    case Br:
      Zc = 28;
      C1 = 15.224;
      C2 = 2.0910;
      C3 = 4.0707;
      break;
    default:
      printf("Error: LCZ pseudopotential not available for %s\n",
	     elements_to_string(atom));
    }
}

gsl_matrix_complex *exp_LCZ_atom::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_lcz_def.hh"
  return output;
}


// ########################### exp_overlap ###########################
gsl_matrix_complex *exp_overlap::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "exp_overlap_def.hh"
  return output;
}
