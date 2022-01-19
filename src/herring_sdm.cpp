////////////////////////////////////////////////////
// Description:
// Model for the spatial CPUE by age i.e CPUE(age) at a given location
// Each age class follow a Tweedie distribution (or we can decide per age group) 
// AR1 process for between age correlation (distance is calculated as age difference)
// Each age group has its own spatial field

// Author: 
// Kotaro Ono
//
// Version: 1
// Detail: 
//
// To do: 
// 1. for the fixed effect design matrix, think about including a spline structure (if needed) - if so, not forget about penalization
// 2. If I want the annual barrier effect, then don't use the average spatial field


#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Function to important barrier-SPDE code
template<class Type>
struct spde_barrier_t{
  vector<Type> C0;
  vector<Type> C1;
  Eigen::SparseMatrix<Type> D0;
  Eigen::SparseMatrix<Type> D1;
  Eigen::SparseMatrix<Type> I;
  spde_barrier_t(SEXP x){           // x = List passed from R 
    C0 = asVector<Type>(getListElement(x,"C0"));
    C1 = asVector<Type>(getListElement(x,"C1"));
    D0 = tmbutils::asSparseMatrix<Type>(getListElement(x,"D0"));
    D1 = tmbutils::asSparseMatrix<Type>(getListElement(x,"D1"));
    I = tmbutils::asSparseMatrix<Type>(getListElement(x,"I"));
  }
};

// Function to calculate Q (precision) matrix using barrier-SPDE
template<class Type>
Eigen::SparseMatrix<Type> Q_spde(spde_barrier_t<Type> spde_barrier, Type kappa, vector<Type> c){
  //using namespace Eigen;
  vector <Type> range(2);
  range(0) = sqrt(8.0)/kappa*c(0);
  range(1) = range(0)*c(1);
  Type pi = 3.141592;
  
  int dimLatent = spde_barrier.D0.row(0).size();
  vector<Type> Cdiag(dimLatent);
  Eigen::SparseMatrix<Type > Cinv(dimLatent,dimLatent);
  
  Cdiag = spde_barrier.C0*pow(range(0),2.0) + spde_barrier.C1*pow(range(1),2.0);
  for(int i =0; i<dimLatent; ++i){
    Cinv.coeffRef(i,i) = 1/Cdiag(i);
  }
  
  Eigen::SparseMatrix<Type>A = spde_barrier.I;
  A = A + (pow(range(0),2.0)/8.0) * spde_barrier.D0 + (pow(range(1),2.0)/8.0) * spde_barrier.D1;
  
  Eigen::SparseMatrix<Type> Q = A.transpose() * Cinv * A/pi *2 * 3;
  
  return Q;
}

// Repeat vector  
template <class Type>
vector<Type> RepeatVector(vector<Type> x, int times)
{
  int n = x.size() * times;
  vector<Type> res(n);
  int k = 0;
  for (int i = 0; i < times; i++) {
    for (int j = 0; j < x.size(); j++) {
      res[k] = x(j);
      k++;
    }
  }
  return res;
}

// Parameter transform for the autocorrelation coefficient
// approach 1
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(0.5) * x)) - Type(1);}

// approach 2
template <class Type>
Type zerofive_to_one(Type x)
{
  return Type(1) - Type(0.5) * invlogit(x);
}


// some likelihood functions
template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0)
{
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type dskewednorm(Type x, Type mean, Type sigma, Type df, int give_log = 0)
{
  // from 
  Type logres = dsn((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0)
{
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}


// some specifications of available options
enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2,
  inverse_link  = 3
};

template <class Type>
Type InverseLink(Type eta, int link)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = exp(eta);
    break;
  case logit_link:
    out = invlogit(eta);
    break;
  case inverse_link:
    out = Type(1.0) / eta;
    break;
  default:
    error("Link not implemented.");
  }
  return out;
}

enum valid_family {
  tweedie_family  = 0,
  // poisson_family    = 2,
  gamma_family      = 1,
  // nb_family         = 4,
  student_family    = 2,
  lognormal_family  = 3,
  gaussian_family   = 4,
  skewednormal_family = 5
};



template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  
  // Data section
  DATA_MATRIX(X);                             //Design matrix for the fixed effects (same for each age groups)
  DATA_VECTOR(yobs);                          //the observed age samples at each station  
  DATA_INTEGER(Nobs);                         //Number of observations
  // Barrier effect
  DATA_STRUCT(spde_barrier,spde_barrier_t);   		//Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier1994,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier1995,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier1996,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier1997,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier1998,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier1999,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2000,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2001,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2002,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2003,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2004,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2005,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2006,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2007,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2008,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_STRUCT(spde_barrier2009,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2009
  DATA_STRUCT(spde_barrier2010,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2010
  DATA_STRUCT(spde_barrier2011,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2011
  DATA_STRUCT(spde_barrier2012,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2012
  DATA_STRUCT(spde_barrier2013,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2013
  DATA_STRUCT(spde_barrier2014,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2014
  DATA_STRUCT(spde_barrier2015,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2015
  DATA_STRUCT(spde_barrier2016,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2016
  DATA_STRUCT(spde_barrier2017,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2017
  DATA_STRUCT(spde_barrier2018,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2018
  DATA_STRUCT(spde_barrier2019,spde_barrier_t);   //Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2019
  DATA_VECTOR(barrier_scaling);               //scaling of range
  
  // INLA feature (how to move from mesh to actual data point)      
  DATA_SPARSE_MATRIX(Aobs);                   //Matrix for interpolating points within triangles (for the observation) - used for the mean spatial field
  DATA_SPARSE_MATRIX(Ast);                    //Same but now divided for each year (each observation in a year) - used for the spatio-temporal field
  DATA_IVECTOR(A_spatial_index);              //Vector of stations to match up A_st output
  DATA_FACTOR(year_i);                        //Year index for the stations 
  DATA_INTEGER(Nyear);                        //Number of years (needed if including the spatio-temporal effect)
  DATA_INTEGER(Nmesh);                        //Number of mesh
  DATA_INTEGER(Npred);                        //Number of prediction points
  DATA_INTEGER(do_predict);                   //Doing prediction or not (save time if just exploring) 0 = not, 1 = yes
  DATA_INTEGER(calc_se);                      //Calculating SD around P(age) (save time if just exploring) 0 = not, 1 = yes
  DATA_MATRIX(X_proj);                        //The design matrix for the prediction
  DATA_SPARSE_MATRIX(A_proj);                 //The A matrix for the prediction
  // Distribution
  DATA_INTEGER(family);
  DATA_INTEGER(link);
  DATA_INTEGER(AR1_spatial);                  //Include AR1 structure to the spatial spatio-temporal field? 0 = not, 1 = yes
  DATA_SCALAR(df);                            //the df for the student t distribution
  DATA_INTEGER(exclude_RE_pred);              //whether to exclude the RE in the prediction
  
  // Parameters to estimate       
  PARAMETER_VECTOR(beta);                     //coefficient associated with the predictors: Npredictor
  PARAMETER_VECTOR(omega);                 //The mean spatial effects (by age class): Nmesh (not used for annual barrier model)
  PARAMETER_ARRAY(epsilon_st);                //The spatio-temporal effect: Nmesh x Ntime
  PARAMETER(transf_rho);                      //The autocorrelation value between year
  
  PARAMETER(logKappa);                 //Spatial scale parameter in Matern covariance structures
  PARAMETER(logTauO);                  //Precision parameter for the average spatial field
  PARAMETER(logTauE);                  //Precision parameter for the spatio-temporal variation 

  PARAMETER(thetaf);                   // tweedie only
  PARAMETER(ln_phi);                   // sigma / dispersion / etc. i.e parameter related to observation error
  
  // derived parameters
  Type rho=f(transf_rho);
  Type range = sqrt(Type(8.0)) / exp(logKappa);
  Type kappa = exp(logKappa);
  Type sigma_E = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * logTauE) * exp(Type(2.0) * logKappa));
  Type sigma_O = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * logTauO) * exp(Type(2.0) * logKappa));
  
  
  // ======================== Calculate the linear predictors (link space) then back-transform ========================
  
  vector<Type> eta(Nobs);     // this is at the link scale 
  vector<Type> mu(Nobs);      // this is at the real scale 
  
  // Step 1: Add the contribution of the fixed effects
    eta = X * beta;

  // Step 2: Add the contribution of the spatial random field (on the intercept)
  // Here we are "projecting" the spatiotemporal and spatial random effects to the
  // locations of the data using the INLA 'A' matrices, Aobs. Because everything is calculated at the mesh vertice level in INLA
    vector<Type> omega_A(Nobs);
    omega_A = Aobs* omega;
    eta += omega_A;

  // Step 3: Add the contribution of the spatio-temporal effect (nrow = unique observation locations of ALL time)
  // Same using the A matrices but a little more complicated because interporlation needs to be done for each year's observation separately
    // Begin with calculating the effect at each observation level across year and age
      matrix<Type> epsilon_st_A(Ast.rows(), Nyear);
      for (int i = 0; i < Nyear; i++){
        epsilon_st_A.col(i) = Ast * vector<Type>(epsilon_st.col(i));
      }
  
    // now finding the YEAR the observation takes place and assign that effect value from epsilon_st_A
      for (int i = 0; i < Nobs; i++){
        eta(i) += epsilon_st_A(A_spatial_index(i), year_i(i));
      }

  // Step 4: back transform to the real scale
    for (int i = 0; i < Nobs; i++){
      mu(i) = InverseLink(eta(i), link);
    }
    
  
  // ======================== The likelihood components ========================
  
  // Defining the NLL
  Type NLL = 0;
  
  // The spatial random effect
    Eigen::SparseMatrix<Type> Q = Q_spde(spde_barrier, kappa, barrier_scaling);
    NLL += SCALE(GMRF(Q,false), 1.0/exp(logTauO))(omega);

  // The spatio-temporal random effect (same kappa = same spatial range)
    Eigen::SparseMatrix<Type> Q1= Q_spde(spde_barrier1994, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q2= Q_spde(spde_barrier1995, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q3= Q_spde(spde_barrier1996, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q4= Q_spde(spde_barrier1997, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q5= Q_spde(spde_barrier1998, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q6= Q_spde(spde_barrier1999, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q7= Q_spde(spde_barrier2000, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q8= Q_spde(spde_barrier2001, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q9= Q_spde(spde_barrier2002, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q10= Q_spde(spde_barrier2003, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q11= Q_spde(spde_barrier2004, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q12= Q_spde(spde_barrier2005, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q13= Q_spde(spde_barrier2006, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q14= Q_spde(spde_barrier2007, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q15= Q_spde(spde_barrier2008, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q16= Q_spde(spde_barrier2009, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q17= Q_spde(spde_barrier2010, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q18= Q_spde(spde_barrier2011, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q19= Q_spde(spde_barrier2012, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q20= Q_spde(spde_barrier2013, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q21= Q_spde(spde_barrier2014, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q22= Q_spde(spde_barrier2015, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q23 = Q_spde(spde_barrier2016, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q24 = Q_spde(spde_barrier2017, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q25 = Q_spde(spde_barrier2018, kappa, barrier_scaling);
    Eigen::SparseMatrix<Type> Q26 = Q_spde(spde_barrier2019, kappa, barrier_scaling);

    if (AR1_spatial== 1) NLL += SCALE(SEPARABLE(AR1(rho), GMRF(Q,false)), 1.0/exp(logTauE))(epsilon_st);
    
    if (AR1_spatial== 0) {
      for (int t = 0; t<Nyear; t++){
        if (t==0) {
          NLL += SCALE(GMRF(Q1,false), 1.0/exp(logTauE))(epsilon_st.col(0));
        }
        if (t==1) {
          NLL += SCALE(GMRF(Q2,false), 1.0/exp(logTauE))(epsilon_st.col(1));
        }
        if (t==2) {
          NLL += SCALE(GMRF(Q3,false), 1.0/exp(logTauE))(epsilon_st.col(2));
        }
        if (t==3) {
          NLL += SCALE(GMRF(Q4,false), 1.0/exp(logTauE))(epsilon_st.col(3));
        }
        if (t==4) {
          NLL += SCALE(GMRF(Q5,false), 1.0/exp(logTauE))(epsilon_st.col(4));
        }
        if (t==5) {
          NLL += SCALE(GMRF(Q6,false), 1.0/exp(logTauE))(epsilon_st.col(5));
        }
        if (t==6) {
          NLL += SCALE(GMRF(Q7,false), 1.0/exp(logTauE))(epsilon_st.col(6));
        }
        if (t==7) {
          NLL += SCALE(GMRF(Q8,false), 1.0/exp(logTauE))(epsilon_st.col(7));
        }
        if (t==8) {
          NLL += SCALE(GMRF(Q9,false), 1.0/exp(logTauE))(epsilon_st.col(8));
        }
        if (t==9) {
          NLL += SCALE(GMRF(Q10,false), 1.0/exp(logTauE))(epsilon_st.col(9));
        }
        if (t==10) {
          NLL += SCALE(GMRF(Q11,false), 1.0/exp(logTauE))(epsilon_st.col(10));
        }
        if (t==11) {
          NLL += SCALE(GMRF(Q12,false), 1.0/exp(logTauE))(epsilon_st.col(11));
        }
        if (t==12) {
          NLL += SCALE(GMRF(Q13,false), 1.0/exp(logTauE))(epsilon_st.col(12));
        }
        if (t==13) {
          NLL += SCALE(GMRF(Q14,false), 1.0/exp(logTauE))(epsilon_st.col(13));
        }
        if (t==14) {
          NLL += SCALE(GMRF(Q15,false), 1.0/exp(logTauE))(epsilon_st.col(14));
        }
        if (t==15) {
          NLL += SCALE(GMRF(Q16,false), 1.0/exp(logTauE))(epsilon_st.col(15));
        }
        if (t==16) {
          NLL += SCALE(GMRF(Q17,false), 1.0/exp(logTauE))(epsilon_st.col(16));
        }
        if (t==17) {
          NLL += SCALE(GMRF(Q18,false), 1.0/exp(logTauE))(epsilon_st.col(17));
        }
        if (t==18) {
          NLL += SCALE(GMRF(Q19,false), 1.0/exp(logTauE))(epsilon_st.col(18));
        }
        if (t==19) {
          NLL += SCALE(GMRF(Q20,false), 1.0/exp(logTauE))(epsilon_st.col(19));
        }
        if (t==20) {
          NLL += SCALE(GMRF(Q21,false), 1.0/exp(logTauE))(epsilon_st.col(20));
        }
        if (t==21) {
          NLL += SCALE(GMRF(Q22,false), 1.0/exp(logTauE))(epsilon_st.col(21));
        }
        if (t==22) {
          NLL += SCALE(GMRF(Q23,false), 1.0/exp(logTauE))(epsilon_st.col(22));
        }
        if (t==23) {
          NLL += SCALE(GMRF(Q24,false), 1.0/exp(logTauE))(epsilon_st.col(23));
        }
        if (t==24) {
          NLL += SCALE(GMRF(Q25,false), 1.0/exp(logTauE))(epsilon_st.col(24));
        }
        if (t==25) {
          NLL += SCALE(GMRF(Q26,false), 1.0/exp(logTauE))(epsilon_st.col(25));
        }
      }
    } 


  // The observation likelihood: Only tweedie at the moment
    Type s1 =0;
    Type s2 =0; 
    Type phi = exp(ln_phi);
    for (int i=0; i<Nobs; i++) {
      if (!isNA(yobs(i))) {
        switch (family) {
          case tweedie_family:
            s1 = invlogit(thetaf) + Type(1.0);
            NLL -= dtweedie(yobs(i), mu(i), phi, s1, true);
            break;
         case student_family:
            NLL -= dstudent(yobs(i), mu(i), phi, df, true) ;
            break;
         case lognormal_family:
            NLL -= dlnorm(yobs(i), log(mu(i)) - pow(phi, Type(2)) / Type(2), phi, true);
            break;
         case gamma_family:
            s2 = mu(i) / phi;  // scale
            NLL -= dgamma(yobs(i), phi, s2, true);
            break;
         case gaussian_family:
            NLL -= dnorm(yobs(i), mu(i), phi, true);
            break;
         case skewednormal_family:
            NLL -= dskewednorm(yobs(i), mu(i), phi, df, true) ;
            break;
         default:
            error("Family not implemented.");
        }
      }
    }

  
  // ======================== Prediction on new data ========================
  
  if (do_predict==1) {
    
    matrix<Type> eta_proj(Npred, Nyear);   // Combined projection in link space
    matrix<Type> mu_proj(Npred, Nyear);    // combined projection in probability scale
    vector<Type> omega_A_proj(Npred);
    matrix<Type> epsilon_st_A_proj(Npred, Nyear);
    eta_proj.setZero();
    mu_proj.setZero();
    omega_A_proj.setZero();
    epsilon_st_A_proj.setZero();
    
    // the fixed effects
    matrix<Type> fixed_proj(Npred, Nyear);
    for (int i = 0; i < Nyear; i++){
      fixed_proj.col(i) = matrix<Type>(X_proj.block(i*Npred, 0, Npred, beta.size())) * beta;
    }

    if (exclude_RE_pred == 0) {
      // The random effects
      // the average spatial field
      omega_A_proj = A_proj * omega;
  
      // the spatio-temporal field
      for (int i = 0; i < Nyear; i++){
        epsilon_st_A_proj.col(i) = A_proj * vector<Type>(epsilon_st.col(i));
      }
    }
    
    // Now add all the effect together
    for (int i = 0; i < Nyear; i++){
      eta_proj.col(i) = vector<Type>(fixed_proj.col(i)) + vector<Type>(epsilon_st_A_proj.col(i)) + omega_A_proj;
    }

    // Now backtransform the variable
    for (int i = 0; i < Npred; i++){
      for (int k = 0; k < Nyear; k++){
        mu_proj(i,k) = InverseLink(eta_proj(i,k), link);
      }
    }
    
    // Report outputs  
    REPORT(fixed_proj);         // the predicted fixed effect
    REPORT(omega_A_proj);       // the predicted average spatial field
    REPORT(epsilon_st_A_proj);  // the predicted spatio-temporal changed in spatial field
    REPORT(eta_proj);           // the predicted combined effect in link space
    REPORT(mu_proj);            // the predicted map of P(age) over space and time
    
    if (calc_se == 1) {
      ADREPORT(mu_proj);
    }
    
  }
  
  // ======================== Reporting ========================
  
  // Model parameters
  REPORT(beta);
  REPORT(logTauE);
  REPORT(sigma_E);
  REPORT(range);
  
  // Derived parameters on species distribution
  REPORT(epsilon_st);
  REPORT(epsilon_st_A);

  // REPORT(Q);		
  // REPORT(Q1);
  // REPORT(Q3);		
  
  // Some outputs for the diagnostics
  REPORT(mu);
  // REPORT(nll_st);
  
  return NLL;
  
}
