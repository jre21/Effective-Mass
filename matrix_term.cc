#include <math.h>
#include <stdio.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

#include "matrix_term.hh"

// size of matrix blocks.  This is the number of angular momentum
// states (four spin-3/2 states plus two spin-1/2 states) multiplied
// by the number of spherical harmonics used (currently one l=0 plus
// five l=2).  This should only be modified if the mathematica-
// generated matrix blocks are changed to include more spherical
// harmonics.
#define BLOCK_SIZE 36

// ########################### matrix_term ###########################
matrix_term::matrix_term()
{
  this->inv_radius = this->dielectric = 1.0;
}

matrix_term::~matrix_term()
{
  this->on_delete();
}

gsl_matrix_complex *matrix_term::matrix(double min, double max, size_t num)
{
  gsl_matrix_complex *work, *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE * num, BLOCK_SIZE * num);
  gsl_matrix_complex_view view;

  // rescale min and max
  min *= this->inv_radius;
  max *= this->inv_radius;

  for(size_t i = 0; i < num; i++)
    for(size_t j = 0; j < num; j++)
      {
	// copy contents of each block into output
	work = this->matrix_block
	  (
	   min * pow(max/min,(double)i/(num-1)),
	   min * pow(max/min,(double)j/(num-1))
	  );
	view = gsl_matrix_complex_submatrix
	  (output, BLOCK_SIZE * i, BLOCK_SIZE * j, BLOCK_SIZE, BLOCK_SIZE);
	gsl_matrix_complex_memcpy(&view.matrix, work);
	gsl_matrix_complex_free(work);
      }
  return output;
}

double matrix_term::get_inv_radius()
{
  return this->inv_radius;
}

double matrix_term::set_inv_radius(double r)
{
  this->on_set_inv_radius(r);
  return this->inv_radius = r;
}

double matrix_term::get_dielectric_constant()
{
  return this->dielectric;
}

double matrix_term::set_dielectric_constant(double k)
{
  this->on_set_dielectric_constant(k);
  return this->dielectric = k;
}

// hooks provided for use by child classes
void matrix_term::on_set_inv_radius(double r) {}
void matrix_term::on_set_dielectric_constant(double k) {}
void matrix_term::on_delete() {}


// ############################## exp_zb #############################
exp_zb::exp_zb
(double _g1, double _g2, double _g3, double _d0, double _dielectric)
{
  this->g1 = _g1;
  this->g2 = _g2;
  this->g3 = _g3;
  this->d0 = _d0;
  this->dielectric = _dielectric;
  this->inv_radius = 1.0 / _g1 / _dielectric;
}

gsl_matrix_complex *exp_zb::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_zb_def.hh"
  return output;
}


// ############################## exp_wz #############################
exp_wz::exp_wz(double _A1, double _A2, double _A3, double _A4, double _A5,
	       double _A6, double _d1, double _d2, double _d3,
	       double _dielectric)
{
  this->A1 = _A1;
  this->A2 = _A2;
  this->A3 = _A3;
  this->A4 = _A4;
  this->A5 = _A5;
  this->A6 = _A6;
  this->d1 = _d1;
  this->d2 = _d2;
  this->d3 = _d3;
  this->dielectric = _dielectric;
  this->inv_radius = 1.0 / (_A2 + _A4) / _dielectric;
}

gsl_matrix_complex *exp_wz::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_wz_def.hh"
  return output;
}


// ############################# exp_gwz #############################
exp_gwz::exp_gwz(double _A1, double _A2, double _A3, double _B1, double _B2,
		 double _B3, double _C1, double _C2, double _C3, double _D1,
		 double _D2, double _D3, double _d1c, double _d2c,
		 double _d1so, double _d2so, double _d3so, double _dielectric)
{
  this->A1 = _A1;
  this->A2 = _A2;
  this->A3 = _A3;
  this->B1 = _B1;
  this->B2 = _B2;
  this->B3 = _B3;
  this->C1 = _C1;
  this->C2 = _C2;
  this->C3 = _C3;
  this->D1 = _D1;
  this->D2 = _D2;
  this->D3 = _D3;
  this->d1c = _d1c;
  this->d2c = _d2c;
  this->d1so = _d1so;
  this->d2so = _d2so;
  this->d3so = _d3so;
  this->dielectric = _dielectric;
  this->inv_radius = 1.0 / (_B1 + _B2) / _dielectric;
}

gsl_matrix_complex *exp_gwz::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_gwz_def.hh"
  return output;
}


// ############################ exp_coulomb ##########################
exp_coulomb::exp_coulomb()
{
  this->dielectric = this->inv_radius = 1.0;
}

gsl_matrix_complex *exp_coulomb::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_coulomb_def.hh"
  return output;
}


// ############################# exp_wang ############################
exp_wang::exp_wang(double _V, double _ra, double _rb, double _r1)
{
  this->V = _V;
  this->ra = _ra;
  this->rb = _rb;
  this->r1 = _r1;
  this->dielectric = this->inv_radius = 1.0;
}

gsl_matrix_complex *exp_wang::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_well_def.hh"
  gsl_matrix_complex *out2 = output;
  output = gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_dielectric_def.hh"
  gsl_matrix_complex_add(output, out2);
  gsl_matrix_complex_free(out2);
  return output;
}


// ############################# exp_LCZ #############################
exp_LCZ::exp_LCZ
(elements_t _host, elements_t _impurity)
{
  this->host = new exp_LCZ_atom(_host);
  this->impurity = new exp_LCZ_atom(_impurity);
  this->coulomb = new exp_coulomb();
}

gsl_matrix_complex *exp_LCZ::matrix(double min, double max, size_t num)
{
  gsl_matrix_complex *h = this->host->matrix(min, max, num),
    *i = this->impurity->matrix(min, max, num),
    *c = this->coulomb->matrix(min, max, num);
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
  this->host->set_inv_radius(r);
  this->impurity->set_inv_radius(r);
  this->coulomb->set_inv_radius(r);
}

void exp_LCZ::on_set_dielectric_constant(double k)
{
  this->host->set_dielectric_constant(k);
  this->impurity->set_dielectric_constant(k);
  this->coulomb->set_dielectric_constant(k);
}

void exp_LCZ::on_delete()
{
  delete this->host;
  delete this->impurity;
  delete this->coulomb;
}


// ########################## exp_LCZ_atom ###########################
exp_LCZ_atom::exp_LCZ_atom(elements_t atom)
{
  this->set_d_core(atom);
}

void exp_LCZ_atom::set_d_core(elements_t atom)
{
  switch(atom)
    {
    case Li:
      this->Zc = 2;
      this->C1 = 2.9990;
      this->C2 = 1.0472;
      this->C3 = 1.9686;
      break;
    case Be:
      this->Zc = 2;
      this->C1 = 3.0169;
      this->C2 = 1.4229;
      this->C3 = 2.9456;
      break;
    case B:
      this->Zc = 2;
      this->C1 = 3.0315;
      this->C2 = 1.8014;
      this->C3 = 3.9668;
      break;
    case C:
      this->Zc = 2;
      this->C1 = 3.0731;
      this->C2 = 2.1922;
      this->C3 = 5.0062;
      break;
    case N:
      this->Zc = 2;
      this->C1 = 3.0487;
      this->C2 = 2.5638;
      this->C3 = 6.0734;
      break;
    case O:
      this->Zc = 2;
      this->C1 = 3.0500;
      this->C2 = 2.9600;
      this->C3 = 6.8700;
      break;
    case F:
      this->Zc = 2;
      this->C1 = 3.0744;
      this->C2 = 3.3609;
      this->C3 = 7.7131;
      break;
    case Ne:
      this->Zc = 2;
      this->C1 = 3.0765;
      this->C2 =  3.7456;
      this->C3 = 8.7366;
      break;
    case Na:
      this->Zc = 10;
      this->C1 = 8.6691;
      this->C2 = 1.354;
      this->C3 = 2.3326;
      break;
    case Mg:
      this->Zc = 10;
      this->C1 = 8.0022;
      this->C2 = 1.4094;
      this->C3 = 2.6122;
      break;
    case Al:
      this->Zc = 10;
      this->C1 = 8.0223;
      this->C2 = 1.5274;
      this->C3 = 3.0200;
      break;
    case Si:
      this->Zc = 10;
      this->C1 = 8.0929;
      this->C2 = 1.6490;
      this->C3 = 3.4389;
      break;
    case P:
      this->Zc = 10;
      this->C1 = 8.3029;
      this->C2 = 1.8115;
      this->C3 = 3.8615;
      break;
    case S:
      this->Zc = 10;
      this->C1 = 8.4635;
      this->C2 = 1.9634;
      this->C3 = 4.3209;
      break;
    case Cl:
      this->Zc = 10;
      this->C1 = 8.6365;
      this->C2 = 2.1205;
      this->C3 = 4.8235;
      break;
    case Ar:
      this->Zc = 10;
      this->C1 = 8.8225;
      this->C2 = 2.2825;
      this->C3 = 5.4118;
      break;
    case K:
      this->Zc = 18;
      this->C1 = 10.615;
      this->C2 = 0.8133;
      this->C3 = 1.5698;
      break;
    case Ca:
      this->Zc = 18;
      this->C1 = 11.154;
      this->C2 = 0.8898;
      this->C3 = 1.8822;
      break;
    case Sc:
      this->Zc = 18;
      this->C1 = 11.706;
      this->C2 = 0.9809;
      this->C3 = 2.1414;
      break;
    case Ti:
      this->Zc = 18;
      this->C1 = 12.210;
      this->C2 = 1.0726;
      this->C3 = 2.3685;
      break;
    case V:
      this->Zc = 18;
      this->C1 = 12.655;
      this->C2 = 1.1619;
      this->C3 = 2.5793;
      break;
    case Cr:
      this->Zc = 18;
      this->C1 = 12.781;
      this->C2 = 1.2449;
      this->C3 = 2.7634;
      break;
    case Mn:
      this->Zc = 18;
      this->C1 = 13.257;
      this->C2 = 1.3321;
      this->C3 = 2.9744;
      break;
    case Fe:
      this->Zc = 18;
      this->C1 = 13.661;
      this->C2 = 1.4247;
      this->C3 = 3.1611;
      break;
    case Co:
      this->Zc = 18;
      this->C1 = 13.761;
      this->C2 = 1.4997;
      this->C3 = 3.3435;
      break;
    case Ni:
      this->Zc = 18;
      this->C1 = 14.134;
      this->C2 = 1.5925;
      this->C3 = 3.5240;
      break;
    case Cu:
      this->Zc = 18;
      this->C1 = 14.230;
      this->C2 = 1.6794;
      this->C3 = 3.6701;
      break;
    case Zn:
      this->Zc = 18;
      this->C1 = 14.557;
      this->C2 = 1.7601;
      this->C3 = 3.8772;
      break;
    case Ga:
      this->Zc = 28;
      this->C1 = 14.934;
      this->C2 = 1.8452;
      this->C3 = 3.3476;
      break;
    case Ge:
      this->Zc = 28;
      this->C1 = 14.888;
      this->C2 = 1.8923;
      this->C3 = 3.5107;
      break;
    case As:
      this->Zc = 28;
      this->C1 = 14.954;
      this->C2 = 1.9544;
      this->C3 = 3.6166;
      break;
    case Se:
      this->Zc = 28;
      this->C1 = 15.171;
      this->C2 = 2.0278;
      this->C3 = 3.8322;
      break;
    case Br:
      this->Zc = 28;
      this->C1 = 15.224;
      this->C2 = 2.0910;
      this->C3 = 4.0707;
      break;
    default:
      printf("Error: LCZ pseudopotential not available for %s\n",
	     elements_to_string(atom));
    }
}

gsl_matrix_complex *exp_LCZ_atom::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_lcz_def.hh"
  return output;
}


// ########################### exp_overlap ###########################
gsl_matrix_complex *exp_overlap::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(BLOCK_SIZE, BLOCK_SIZE);
  #include "exp_overlap_def.hh"
  return output;
}
