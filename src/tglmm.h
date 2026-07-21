#ifndef ACCMETA_TGLMM_H
#define ACCMETA_TGLMM_H

// The integrand of one study on the log scale
template<class Type>
Type tglmm_h(const vector<Type>& u, const vector<Type>& MU, const vector<Type>& A,
             const vector<int>& yi, const vector<int>& di) {
  Type a1 = MU(0) + A(0) * u(0);
  Type a2 = MU(1) + A(1) * u(0) + A(2) * u(1);
  Type a3 = MU(2) + A(3) * u(0) + A(4) * u(1) + A(5) * u(2);
  return yi(0) * a1 - di(0) * logspace_add(Type(0), a1)
       + yi(1) * a2 - di(1) * logspace_add(Type(0), a2)
       + yi(2) * a3 - di(2) * logspace_add(Type(0), a3)
       - Type(0.5) * (u * u).sum();
}

#undef  TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type tglmm(objective_function<Type>* obj)
{
  // data
  DATA_IVECTOR(y);             //numerators in y1, y2, y3 order
  DATA_IVECTOR(den);           //denominators
  DATA_FACTOR(f);             //group indicator  1,..,n
  DATA_INTEGER(niter);
  DATA_VECTOR(ws);            // these are w * exp(z^2)
  DATA_VECTOR(z);
  DATA_SCALAR(DEGREES);      // Wishart degrees of freedom
  DATA_SCALAR(SCALE);        // Wishart scale

  // params 
  PARAMETER_VECTOR(MU);      // mu_eta, mu_xi, mu_gamma
  PARAMETER_VECTOR(ALPHA);   // log-Cholesky of Sigma: log(L11), L21, log(L22), L31, L32, log(L33)

  using namespace density;

  // Cholesky of Sigma on the natural scale
  vector<Type> A(6);
  A(0) = exp(ALPHA(0));
  A(1) = ALPHA(1);
  A(2) = exp(ALPHA(2));
  A(3) = ALPHA(3);
  A(4) = ALPHA(4);
  A(5) = exp(ALPHA(5));

  Type ll = 0.0;

  int n = y.size() / 3;

  matrix<Type> L(3,3);
  L.setZero();

  int ind = 0;
  for(int i=0;i<3;i++)
    for(int j=0;j<=i;j++) {
      L(i,j) = A(ind);
      ind += 1;
    }
  matrix<Type> Sigma = L * L.transpose();

  // Compute modes
  Type INNER_GRAD = 0.0;   // running max of the squared inner gradient norm
  vector<vector<int> > yspl = split(y, f);
  vector<vector<int> > dspl = split(den, f);
  for(int i = 0; i < n; i++){
    vector<Type> ui(3);
    ui(0) = 0.0;
    ui(1) = 0.0;
    ui(2) = 0.0;
    Type logLhat = 0.0; 
    Type l11 = 0.0; //these l are the Cholesky elements of H^-1
    Type l21 = 0.0;
    Type l22 = 0.0;
    Type l31 = 0.0;
    Type l32 = 0.0;
    Type l33 = 0.0;
    vector<int> yi = yspl(i);
    vector<int> di = dspl(i);
    
    for(int k = 0; k <= niter; k++){
      vector<Type>  g = -1.0 * ui;
      matrix<Type> H(3,3);
      H(0,0) = 1.0;
      H(0,1) = 0.0;
      H(0,2) = 0.0;
      H(1,0) = 0.0;
      H(1,1) = 1.0;
      H(1,2) = 0.0;
      H(2,0) = 0.0;
      H(2,1) = 0.0;
      H(2,2) = 1.0;
      Type arg1 = MU(0) + A(0) * ui(0);
      Type p1 = invlogit(arg1);
      Type arg2 = MU(1) + A(1) * ui(0) + A(2) * ui(1);
      Type p2 = invlogit(arg2);
      Type arg3 = MU(2) + A(3) * ui(0) + A(4) * ui(1) + A(5) * ui(2);
      Type p3 = invlogit(arg3);
      g(0) += A(0) * (yi(0) - p1 * di(0)) +  A(1) * (yi(1) - p2 * di(1)) + A(3) * (yi(2) - p3 * di(2));
      g(1) += A(2) * (yi(1) - p2 * di(1)) + A(4) * (yi(2) - p3 * di(2));
      g(2) += A(5) * (yi(2) - p3 * di(2));
      H(0,0) += di(0) * pow(A(0), 2.0) * p1 * (1 - p1) + di(1) * pow(A(1), 2.0) * p2 * (1 - p2) + di(2) * pow(A(3), 2.0) * p3 * (1 - p3);
      H(1,1) += di(1) * pow(A(2), 2.0) * p2 * (1 - p2) + di(2) * pow(A(4), 2.0) * p3 * (1 - p3);
      H(2,2) += di(2) * pow(A(5), 2.0) * p3 * (1 - p3);
      H(0,1) += di(1) * A(1) * A(2) * p2 * (1 - p2) + di(2) * A(3) * A(4) * p3 * (1 - p3);
      H(0,2) += di(2) * A(3) * A(5) * p3 * (1 - p3);
      H(1,2) += di(2) * A(4) * A(5) * p3 * (1 - p3);
      H(1,0) = H(0,1);
      H(2,0) = H(0,2);
      H(2,1) = H(1,2);
      matrix<Type> H1 = H;
      Type detH = H(0,0) * H(1,1) * H(2,2)
        + Type(2.0) * H(0,1) * H(0,2) * H(1,2)
        - H(0,0) * pow(H(1,2), Type(2.0))
        - H(1,1) * pow(H(0,2), Type(2.0))
        - H(2,2) * pow(H(0,1), Type(2.0));
      detH = CppAD::CondExpLt(detH, Type(1e-12), Type(1e-12), detH);
      H1(0,0) = (H(1,1) * H(2,2) - pow(H(1,2), 2.0))  / detH;
      H1(1,1) = (H(0,0) * H(2,2) - pow(H(0,2), 2.0))  / detH;
      H1(2,2) = (H(0,0) * H(1,1) - pow(H(0,1), 2.0))  / detH;
      H1(0,1) = (H(0,2) * H(1,2) - H(0,1) * H(2,2)) / detH;
      H1(0,2) = (H(0,1) * H(1,2) - H(0,2) * H(1,1)) / detH;
      H1(1,2) = (H(0,1) * H(0,2) - H(0,0) * H(1,2)) / detH;
      H1(1,0) = H1(0,1);
      H1(2,0) = H1(0,2);
      H1(2,1) = H1(1,2);
      
      //backtracking
      if (k < niter) {
        vector<Type> step = H1 * g;
        Type h0 = tglmm_h(ui, MU, A, yi, di);
        Type t = Type(1.0);
        for (int j = 0; j < 20; j++) {
          vector<Type> cand = ui + t * step;
          Type hc = tglmm_h(cand, MU, A, yi, di);
          t = CppAD::CondExpGt(hc, h0, t, Type(0.5) * t);
        }
        ui += t * step;
      }
      // worst inner gradient over studies
      if (k == niter) {
        Type gnorm2 = (g * g).sum();
        INNER_GRAD = CppAD::CondExpGt(gnorm2, INNER_GRAD, gnorm2, INNER_GRAD);
      }
      Type tmp11 = H1(0,0);
      tmp11 = CppAD::CondExpLt(tmp11, Type(1e-12), Type(1e-12), tmp11);
      l11 = sqrt(tmp11);
      l21 = H1(0,1) / l11;
      l31 = H1(0,2) / l11;

      Type tmp22 = H1(1,1) - pow(l21, Type(2.0));
      tmp22 = CppAD::CondExpLt(tmp22, Type(1e-12), Type(1e-12), tmp22);
      l22 = sqrt(tmp22);
      l32 = (H1(1,2) - l31 * l21) / l22;

      Type tmp33 = H1(2,2) - pow(l31, Type(2.0)) - pow(l32, Type(2.0));
      tmp33 = CppAD::CondExpLt(tmp33, Type(1e-12), Type(1e-12), tmp33);
      l33 = sqrt(tmp33);
      logLhat =  log(l11) + log(l22) + log(l33);
    }

    //at this point we got the mode ui and the log(L)


    // negative log-likelihood
    Type logsum = -INFINITY;
    int nq = ws.size();
    for (int j1 = 0; j1 < nq; j1++) {
      for (int j2 = 0; j2 < nq; j2++) {
        for (int j3 = 0; j3 < nq; j3++) {
          Type addint = 0.0;
          Type zj1 = sqrt(Type(2.0)) * l11 * z(j1) + ui(0);
          Type zj2 = sqrt(Type(2.0)) * (l21 * z(j1) +  l22 * z(j2)) + ui(1);
          Type zj3 = sqrt(Type(2.0)) * (l31 * z(j1) +  l32 * z(j2) + l33 * z(j3)) + ui(2);
          Type arg1 = MU(0) + A(0) * zj1;
          Type arg2 = MU(1) + A(1) * zj1 + A(2) * zj2;
          Type arg3 = MU(2) + A(3) * zj1 + A(4) * zj2 + A(5) * zj3;
          addint += yi(0) * arg1 - di(0) * logspace_add(Type(0), arg1);
          addint += yi(1) * arg2 - di(1) * logspace_add(Type(0), arg2);
          addint += yi(2) * arg3 - di(2) * logspace_add(Type(0), arg3);
          addint += dnorm(zj1, Type(0), Type(1), true) + dnorm(zj2, Type(0), Type(1), true) +  dnorm(zj3, Type(0), Type(1), true);
          addint += log(ws[j1]) + log(ws[j2]) + log(ws[j3]);
          logsum = logspace_add(logsum, addint);
        }
      }
    }

    // Once per study, not once per node: the substitution is
    // u = u_hat + sqrt(2) L z, whose Jacobian is 2^(d/2) |L| for d = 3, and
    // exp(logLhat) supplies |L| alone.
    ll += logsum + logLhat + Type(1.5) * log(Type(2.0));
  }
  Type nll = ll * (-1.0);

  Type INNER_GRAD_MAX = sqrt(INNER_GRAD);
  REPORT(INNER_GRAD_MAX);

  // Wishart(DEGREES, A*I_3) prior. 
  nll -= (DEGREES - Type(4.0)) * (ALPHA(0) + ALPHA(2) + ALPHA(5));
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
