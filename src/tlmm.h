#ifndef ACCMETA_TLMM_H
#define ACCMETA_TLMM_H

#undef  TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type tlmm(objective_function<Type>* obj)
{
  // Data
  DATA_MATRIX(EST);          // n x 3: eta, xi, gamma
  DATA_MATRIX(WVAR);         // n x 3: var_eta, var_xi, var_gamma
  DATA_SCALAR(DEGREES);      // Wishart degrees of freedom
  DATA_SCALAR(SCALE);        // Wishart scale
  DATA_STRING(WORKPAR);      // "PinheiroBates" or "Joe"

  // params
  PARAMETER_VECTOR(MU);      // mu_eta, mu_xi, mu_gamma
  PARAMETER_VECTOR(ALPHA);   // PinheiroBates: log(L11), L21, log(L22), L31, L32, log(L33)
                             // Joe: log(s1), log(s2), log(s3), atanh(r12), atanh(r13), atanh(r23|1)

  using namespace density;

  // log-Cholesky of Sigma
  vector<Type> LOGCHOL = ALPHA;
  if (WORKPAR == "Joe") {
    LOGCHOL(1) = exp(ALPHA(1)) * tanh(ALPHA(3));
    LOGCHOL(2) = ALPHA(1) - log(cosh(ALPHA(3)));
    LOGCHOL(3) = exp(ALPHA(2)) * tanh(ALPHA(4));
    LOGCHOL(4) = exp(ALPHA(2)) * tanh(ALPHA(5)) / cosh(ALPHA(4));
    LOGCHOL(5) = ALPHA(2) - log(cosh(ALPHA(4))) - log(cosh(ALPHA(5)));
  }
  REPORT(LOGCHOL);

  // Sigma mat
  matrix<Type> L(3, 3);
  L.setZero();

  int ind = 0;
  for(int i = 0; i < 3; i++)
    for(int j = 0; j <= i; j++) {
      L(i, j) = (i == j) ? exp(LOGCHOL(ind)) : LOGCHOL(ind);
      ind += 1;
    }
  matrix<Type> Sigma = L * L.transpose();

  // Negative log-likelihood
  Type nll = 0.0;
  for(int i = 0; i < EST.rows(); i++){
    matrix<Type> Sigmai = Sigma;
    for(int j = 0; j < 3; j++) Sigmai(j, j) += WVAR(i, j);
    vector<Type> r = vector<Type>(EST.row(i)) - MU;
    nll += MVNORM<Type>(Sigmai)(r);
  }

  // Wishart(DEGREES, A*I_3) prior. 
  nll -= (DEGREES - Type(4.0)) * (LOGCHOL(0) + LOGCHOL(2) + LOGCHOL(5));
  nll += Sigma.trace() / (Type(2.0) * SCALE);        

  ADREPORT(Sigma);

  // Sensitivity
  Type SE = invlogit(MU(0));
  ADREPORT(SE);

  // Specificity
  Type SP = Type(1.0) - invlogit(MU(1));
  ADREPORT(SP);

  // Prevalence
  Type PREV = invlogit(MU(2));
  ADREPORT(PREV);

  return nll;
}

#undef  TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
