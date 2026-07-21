#define TMB_LIB_INIT R_init_accmeta
#include <TMB.hpp>
#include "tlmm.h"
#include "tglmm.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(MODEL);
  if (MODEL == "tlmm") return tlmm(this);
  if (MODEL == "tglmm") return tglmm(this);
  Rf_error("unknown MODEL '%s'", MODEL.c_str());
  return Type(0);
}
