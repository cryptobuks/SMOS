#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//////////////////////////////////////////////////////////////////////
// Component functions
//////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
mat at(mat m, mat G) {
  return((G*m.t()).t());
}

// [[Rcpp::export]]
mat Rt(mat G, mat C, mat W) {
  return(G*C*G.t() + W);
}

// [[Rcpp::export]]
mat ft(mat F, mat a) {
  return((F*a.t()).t());
}

// [[Rcpp::export]]
mat Qt(mat F, mat R, mat V) {
  return(F*R*F.t() + V);
}

// [[Rcpp::export]]
mat mt(double Y, mat F, mat a, mat R, mat f, mat Q) {
  return((a.t() + R*F.t()*inv(Q)*(Y - f)).t());
}

// [[Rcpp::export]]
mat Ct(mat F, mat R, mat Q) {
  return(R - R*F.t()*inv(Q)*F*R);
}

// [[Rcpp::export]]
mat ht(mat G, mat C, mat R, mat theta, mat a, mat m) {
  return((m.t() + C*G.t()*inv(R)*(theta - a).t()).t());
}

// [[Rcpp::export]]
mat Ht(mat C, mat G, mat R) {
  return(C - C*G.t()*inv(R)*G*C);
}

// [[Rcpp::export]]
mat ht2(mat G, mat C, mat W, mat theta, mat H, mat m) {
  return((H*(inv(C)*m.t() + G.t()*inv(W)*theta.t())).t());
}

// [[Rcpp::export]]
mat Ht2(mat C, mat G, mat W) {
  return(inv(inv(C) + G.t()*inv(W)*G));
}

// [[Rcpp::export]]
int nomiss(mat y) {
  int T = y.n_rows;
  int sum=0;
  Rcout << sum << std::endl;
  for (int i=0; i<T; i++){
    sum += !R_IsNA(y[i]);
  }
  return sum;
}

//////////////////////////////////////////////////////////////////////
// Drawing the Multivariate Normal
//////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
mat mvrnormChol(int n, mat mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n) + Y * chol(sigma);
}

// [[Rcpp::export]]
mat mvrnormSVD(int n, mat mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  mat U;
  vec d;
  mat V;
  svd(U, d, V, sigma);
  return repmat(mu, 1, n) + Y * U*diagmat(sqrt(d))*V;
}

// [[Rcpp::export]]
double normden(double x, double mu, double sd) {
  return (1/(pow(2*M_PI,0.5)*sd))*exp(-0.5*pow((x-mu)/sd, 2));
}

// [[Rcpp::export]]
mat rwishart(int df, mat S){
  // Dimension of returned wishart
  int m = S.n_rows;
  
  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  mat Z(m,m);
  
  // Fill the diagonal
  for(int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  
  // Fill the lower matrix with random guesses
  for(int j = 0; j < m; j++){  
    for(int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  // Lower triangle * chol decomp
  mat C = trimatl(Z).t() * chol(S);
  
  // Return random wishart
  return C.t()*C;
}

//////////////////////////////////////////////////////////////////////
// Kalman Filter
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List KalmanFilter(vec Y, mat V, mat W, List model){
  // Set initial values
  mat F = as<mat>(model["FF"]);
  mat G = as<mat>(model["GG"]);
  //  mat V = as<mat>(model["V"]);
  //  mat W = as<mat>(model["W"]);
  vec m0 = as<vec>(model["m0"]);
  mat C0 = as<mat>(model["C0"]);
  
  
  // Try SVD on variances for speed
  int T = Y.n_rows;
  int p = G.n_rows;
  
  // storing output
  mat mt_keep(T+1,p);
  mat at_keep(T,p);
  mat ft_keep(T,1);
  List Rt_keep(T);
  List Qt_keep(T);
  List Ct_keep(T+1);
  
  // reused during kalman filter
  mat R;
  mat C;
  mat Q;
  
  mt_keep.row(0) = m0.t();
  Ct_keep[0] = C0;
  C = C0;
  
  for (int i=0; i<T; i++){
    at_keep.row(i) = at(mt_keep.row(i), G);
    R = Rt(G, C, W);
    
    ft_keep.row(i) = ft(F, at_keep.row(i));
    Q = Qt(F, R, V);
    
    if (R_IsNA(Y[i])) {
      mt_keep.row(i+1) = at_keep.row(i);
      C = R;
    } else {
      mt_keep.row(i+1) = mt(Y[i], F, at_keep.row(i), R, ft_keep.row(i), Q);
      C = Ct(F, R, Q);
    }
    Rt_keep[i] = R;
    Qt_keep[i] = Q;
    Ct_keep[i+1] = C;
  }
  return List::create(
    Named("at") = at_keep,
    Named("ft") = ft_keep,
    Named("mt") = mt_keep,
    Named("Rt") = Rt_keep,
    Named("Qt") = Qt_keep,
    Named("Ct") = Ct_keep);
}

//////////////////////////////////////////////////////////////////////
// Fast Forward Backward Sampling (FFBS)
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List FFBS(mat G, mat W, List Kalman, bool SVD = false){
  mat at = as<mat>(Kalman["at"]);
  mat mt = as<mat>(Kalman["mt"]);
  List Rt = as<List>(Kalman["Rt"]);
  List Ct = as<List>(Kalman["Ct"]);
  
  int nrow = at.n_rows;
  int T = nrow - 1;
  int p = G.n_rows;
  
  // storing output
  mat ht_keep(nrow,p);
  mat theta_keep(nrow+1,p);
  List Ht_keep(nrow);
  
  // Draw theta_T from filtering distribution
  if(SVD) {
    theta_keep.row(T+1) = mvrnormSVD(1, mt.row(T+1), Ct[T+1]);
  } else {
    theta_keep.row(T+1) = mvrnormChol(1, mt.row(T+1), Ct[T+1]);
  }
  
  for (int i=1; i<=nrow; i++){
    Ht_keep[nrow-i] = Ht2(Ct[nrow-i], G, W);
    ht_keep.row(nrow-i) = ht2(G, Ct[nrow-i], W, theta_keep.row(nrow-i+1), Ht_keep[nrow-i], mt.row(nrow-i));
    
    // Draw theta from posterior
    if (SVD) {
      theta_keep.row(nrow-i) = mvrnormSVD(1, ht_keep.row(nrow-i), Ht_keep[nrow-i]);
    } else {
      theta_keep.row(nrow-i) = mvrnormChol(1, ht_keep.row(nrow-i), Ht_keep[nrow-i]);
    }
  }
  return List::create(
    Named("theta") = theta_keep,
    Named("ht") = ht_keep,
    Named("Ht") = Ht_keep);
}

//////////////////////////////////////////////////////////////////////
// Calculating Sums of Squares
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
mat calc_SSE_error(mat y, mat theta, mat F) {
  int T = y.n_rows;
  mat SSE(1,1);
  SSE.zeros();
  for (int i=0; i<T; i++) {
    if (R_IsNA(y[i])) {
      SSE += 0;
    } else {
      SSE += (y[i]-F*theta.row(i+1).t())*(y[i]-F*theta.row(i+1).t());
    }
  }
  return SSE;
}

// [[Rcpp::export]]
mat calc_SSE_theta(mat theta, mat G) {
  int T = theta.n_rows;
  mat SSE(1,1);
  SSE.zeros();
  for (int i=1; i<T; i++) 
    SSE += ((theta.row(i).t()-G*theta.row(i-1).t()).t()*(theta.row(i).t()-G*theta.row(i-1).t()));
  return SSE;
}

// [[Rcpp::export]]
mat calc_SSE_theta_wishart(mat theta, mat G) {
  int T = theta.n_rows;
  int p = theta.n_cols;
  mat SSE(p,p);
  SSE.zeros();
  for (int i=1; i<T; i++) 
    SSE += (theta.row(i).t()-G*theta.row(i-1).t())*(theta.row(i).t()-G*theta.row(i-1).t()).t();
  return SSE;
}

//////////////////////////////////////////////////////////////////////
// Sampling Functions
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double cauchy_MH(double par_current, double shape, double rate, double c) {
  double par_proposed = 1/R::rgamma(shape, 1/rate);
  double log_rho = log(1+par_current/(c*c)) - log(1+par_proposed/(c*c));
  return ifelse(log(runif(1)) < log_rho, par_proposed, par_current)[0];
}

// [[Rcpp::export]]
double log_RW(double par_current, double shape, double SSE, double propsd, double c) {
  double lcurrent = log(par_current);
  double lproposed = lcurrent + propsd*R::rnorm(0,1);
  double par_proposed = exp(lproposed);
  
  double log_rho = -log(1 + pow(par_proposed/c,2)) - shape*log(par_proposed) - 0.5*pow(par_proposed,-2)*SSE + 
    log(1 + pow(par_current/c,2)) + shape*log(par_current) + 0.5*pow(par_current,-2)*SSE + 
    log(par_proposed) - log(par_current);
  return ifelse(log(runif(1)) < log_rho, par_proposed, par_current)[0];
} 

// [[Rcpp::export]]
double trunc_MHc(double par_current, double shape, double SSE, double propsd, double c) {
  double par_proposed = -1.0;
  while(par_proposed <= 0.0){
    par_proposed = R::rcauchy(par_current, propsd);
  }
  double log_rho = -log(1 + pow(par_proposed,2)/(c*c)) - shape*log(par_proposed) - 1/(2*pow(par_proposed,2))*SSE + 
    log(1 + pow(par_current,2)/(c*c)) + shape*log(par_current) + 1/(2*pow(par_current,2))*SSE + 
    log(R::pcauchy(par_current, 0.0, propsd, 1, 0)) -
    log(R::pcauchy(par_proposed, 0.0, propsd, 1, 0)); 
  return ifelse(log(runif(1)) < log_rho, par_proposed, par_current)[0];
}

//////////////////////////////////////////////////////////////////////
// MCMC
//////////////////////////////////////////////////////////////////////

// Cauchy priors, MH - IG proposal
// [[Rcpp::export]]
List mcmc_seas_IG(
    int n_reps, 
    mat dat, 
    List initial_values,
    NumericVector c,
    List model,
    bool svd = false,
    bool save_theta = true){
  
  // Set initial values
  double sig2_e = as<double>(initial_values["sig2_e"]);
  double sig2_w = as<double>(initial_values["sig2_w"]);
  
  // Model things
  mat F = as<mat>(model["F"]);
  mat G = as<mat>(model["G"]);
  mat W = as<mat>(model["W"]);
  mat V = as<mat>(model["V"]);
  
  
  List keep_theta(n_reps);
  NumericVector keep_sigma_e(n_reps);
  NumericVector keep_sigma_w(n_reps);
  
  int T = dat.n_rows;
  int p = G.n_cols;
  
  //reused during sampling
  List kalman;
  mat rate1, rate2;
  double sh1 = (T-1)/2, sh2 = (T*p-1)/2;
  
  for (int i=0; i<n_reps; i++) {
    kalman = KalmanFilter(dat, V, W, model); // Kalman Filter
    mat theta = FFBS(G, W, kalman, svd)[0]; // draw theta_0:T
    
    rate1 = calc_SSE_error(dat, theta, F)/2;
    double rate1b = rate1[0]; 
    sig2_e = cauchy_MH(sig2_e, sh1, rate1b, c[0]); 
    
    rate2 = calc_SSE_theta(theta, G)/2;
    double rate2b = rate2[0]; 
    sig2_w = cauchy_MH(sig2_w, sh2, rate2b, c[1]); 
    
    // update V    
    V = sig2_e;
    
    // update W
    for (int j=0; j < p; j++){
      W(j,j) = sig2_w;
    }
    
    // Update storage
    if(save_theta){
      keep_theta[i] = theta;
    }
    keep_sigma_e[i] = sig2_e;
    keep_sigma_w[i] = sig2_w;  
  }
  if(save_theta){
    return List::create(
      Named("sigma2_e") = keep_sigma_e,
      Named("sigma2_w") = keep_sigma_w,
      Named("theta") = keep_theta);
  }else{
    return List::create(
      Named("sigma2_e") = keep_sigma_e,
      Named("sigma2_w") = keep_sigma_w);
  }
}

// Conjugate IG priors
// [[Rcpp::export]]
List mcmc_local_seas(
    int n_reps,
    int thin,
    mat dat, 
    List initial_values,
    List prior,
    List model,
    bool svd = false){
  
  // Set initial values
  double psi1 = as<double>(initial_values["psi1"]);
  double psi2 = as<double>(initial_values["psi2"]);
  
  // Set prior values
  double a1 = as<double>(prior["a1"]); 
  double a2 = as<double>(prior["a2"]);
  double b1 = as<double>(prior["b1"]);
  double b2 = as<double>(prior["b2"]); 
  
  // Model things
  mat F = as<mat>(model["FF"]);
  mat G = as<mat>(model["GG"]);
  mat W = as<mat>(model["W"]);
  mat V = as<mat>(model["V"]);
  
  List keep_theta(n_reps);
  NumericVector keep_sigma_e(n_reps);
  NumericVector keep_sigma_w(n_reps);
  
  int T = dat.n_rows;
  int T2=nomiss(dat);
  int p = G.n_cols;
  int every = thin+1;
  int mc = n_reps*every;
  int save = 0;
  //reused during sampling
  List kalman;
  mat rate1, rate2;
  
  double sh1 = a1 + T2/2, sh2 = a2 + T*p/2;
  
  for (int i=0; i<mc; i++) {
    kalman = KalmanFilter(dat, V, W, model); // Kalman Filter
    mat theta = FFBS(G, W, kalman, svd)[0]; // draw theta_0:T
    
    rate1 = b1 + calc_SSE_error(dat, theta, F)/2;
    double rate1b = rate1[0]; 
    psi1 = rgamma(1, sh1, 1/rate1b)[0];
    
    rate2 = b2 + calc_SSE_theta(theta, G)/2;
    double rate2b = rate2[0]; 
    psi2 = rgamma(1, sh2, 1/rate2b)[0];
    
    // update V    
    V = 1/psi1;
    
    // update W
    for (int j=0; j < p; j++){
      W(j,j) = 1/psi2;
    }
    
    // Update storage
    if ((i%every==0)) {
    keep_theta(save) = theta;
    keep_sigma_e(save) = 1/psi1;
    keep_sigma_w(save) = 1/psi2;  
    save=save+1;
    }
  }
  
  return List::create(
    Named("theta") = keep_theta,
    Named("sigma2_e") = keep_sigma_e,
    Named("sigma2_w") = keep_sigma_w);
}


/*** R
#Put R code Here
*/
