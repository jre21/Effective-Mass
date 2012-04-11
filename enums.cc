#include "enums.hh"

const char *elements_to_string(elements_t e)
{
  switch(e)
    {
    case H:  return "H";
    case He: return "He";
    case Li: return "Li";
    case Be: return "Be";
    case B:  return "B";
    case C:  return "C";
    case N:  return "N";
    case O:  return "O";
    case F:  return "F";
    case Ne: return "Ne";
    case Na: return "Na";
    case Mg: return "Mg";
    case Al: return "Al";
    case Si: return "Si";
    case P:  return "P";
    case S:  return "S";
    case Cl: return "Cl";
    case Ar: return "Ar";
    case K:  return "K";
    case Ca: return "Ca";
    case Sc: return "Sc";
    case Ti: return "Ti";
    case V:  return "V";
    case Cr: return "Cr";
    case Mn: return "Mn";
    case Fe: return "Fe";
    case Co: return "Co";
    case Ni: return "Ni";
    case Cu: return "Cu";
    case Zn: return "Zn";
    case Ga: return "Ga";
    case Ge: return "Ge";
    case As: return "As";
    case Se: return "Se";
    case Br: return "Br";
    case Kr: return "Kr";
    case Rb: return "Rb";
    case Sr: return "Sr";
    case Y:  return "Y";
    case Zr: return "Zr";
    case Nb: return "Nb";
    case Mo: return "Mo";
    case Tc: return "Tc";
    case Ru: return "Ru";
    case Rh: return "Rh";
    case Pd: return "Pd";
    case Ag: return "Ag";
    case Cd: return "Cd";
    case In: return "In";
    case Sn: return "Sn";
    case Sb: return "Sb";
    case Te: return "Te";
    case I:  return "I";
    case Xe: return "Xe";
    case Cs: return "Cs";
    case Ba: return "Ba";
    case La: return "La";
    case Ce: return "Ce";
    case Pr: return "Pr";
    case Nd: return "Nd";
    case Pm: return "Pm";
    case Sm: return "Sm";
    case Eu: return "Eu";
    case Gd: return "Gd";
    case Tb: return "Tb";
    case Dy: return "Dy";
    case Ho: return "Ho";
    case Er: return "Er";
    case Tm: return "Tm";
    case Yb: return "Yb";
    case Lu: return "Lu";
    case Hf: return "Hf";
    case Ta: return "Ta";
    case W:  return "W";
    case Re: return "Re";
    case Os: return "Os";
    case Ir: return "Ir";
    case Pt: return "Pt";
    case Au: return "Au";
    case Hg: return "Hg";
    case Tl: return "Tl";
    case Pb: return "Pb";
    case Bi: return "Bi";
    case Po: return "Po";
    case At: return "At";
    case Rn: return "Rn";
    case Fr: return "Fr";
    case Ra: return "Ra";
    case Ac: return "Ac";
    case Th: return "Th";
    case Pa: return "Pa";
    case U:  return "U";
    case Np: return "Np";
    case Pu: return "Pu";
    case Am: return "Am";
    case Cm: return "Cm";
    case Bk: return "Bk";
    case Cf: return "Cf";
    case Es: return "Es";
    case Fm: return "Fm";
    case Md: return "Md";
    case No: return "No";
    case Lr: return "Lr";
    case Rf: return "Rf";
    case Db: return "Db";
    case Sg: return "Sg";
    case Bh: return "Bh";
    case Hs: return "Hs";
    case Mt: return "Mt";
    case Ds: return "Ds";
    case Rg: return "Rg";
    case Cn: return "Cn";
    case Ga_val: return "Ga: valance";
    case Ge_val: return "Ge: valance";
    case As_val: return "As: valance";
    case Se_val: return "Se: valance";
    case Br_val: return "Br: valance";
    default: return "unknown";
    }
}

const char *crystals_to_string(crystals_t c)
{
  switch(c)
    {
    case GaN: return "GaN";
    case AlN: return "AlN";
    case InN: return "InN";
    case ZnGeN2: return "ZnGeN2";
    case ZnSnN2: return "ZnSnN2";
    default: return "unknown";
    }
}

const char *crystal_parameters_to_string(crystal_parameters_t c)
{
  switch(c)
    {
    case _g1: return "g1";
    case _g2: return "g2";
    case _g3: return "g3";
    case _d0: return "d0";
    case _A1: return "A1";
    case _A2: return "A2";
    case _A3: return "A3";
    case _A4: return "A4";
    case _A5: return "A5";
    case _A6: return "A6";
    case _d1: return "d1";
    case _d2: return "d2";
    case _d3: return "d3";
    case _B1: return "B1";
    case _B2: return "B2";
    case _B3: return "B3";
    case _C1: return "C1";
    case _C2: return "C2";
    case _C3: return "C3";
    case _D1: return "D1";
    case _D2: return "D2";
    case _D3: return "D3";
    case _d1c: return "d1c";
    case _d2c: return "d2c";
    case _d1so: return "d1so";
    case _d2so: return "d2so";
    case _d3so: return "d3so";
    case _dielectric: return "dielectric";
    default: return "unknown";
    }
}

const char *impurity_parameters_to_string(impurity_parameters_t c)
{
  switch(c)
    {
    case _dielectric_ratio: return "dielectric_ratio";
    default: return "unknown";
    }
}
