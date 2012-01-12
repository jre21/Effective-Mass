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
    default: return "unknown";
    }
}

const char *crystals_to_string(crystals_t c)
{
  switch(c)
    {
    case GaN: return "GaN";
    case AlN: return "AlN";
    case InN1: return "InN1";
    case InN2: return "InN2";
    case InN3: return "InN3";
    case ZnGeN2: return "ZnGeN2";
    case ZnSnN2: return "ZnSnN2";
    default: return "unknown";
    }
}
