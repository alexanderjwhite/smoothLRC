// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


arma::rowvec vecnorm_row(arma::mat x){
  arma::rowvec norm_x(x.n_cols,fill::zeros);
  for (int i = 0; i < x.n_cols; i++){
    norm_x(i) = norm(x.col(i), 2);
  }
  return(norm_x);
}


arma::mat vecnorm_diag(arma::mat x){
  arma::mat norm_x = zeros<arma::mat>(x.n_cols,x.n_cols);
  for (int i = 0; i < x.n_cols; i++){
    norm_x(i,i) = norm(x.col(i), 2);
  }
  return(norm_x);
}


List uv_norm(arma::mat u,
             arma::mat v){

  List result;

  arma::mat unorm = normalise(u);
  arma::mat vnorm = normalise(v);
  arma::rowvec unorm_vec = vecnorm_row(u);
  arma::rowvec vnorm_vec = vecnorm_row(v);

  for (int i = 0; i < u.n_cols; i++){
    u.col(i) = unorm.col(i) * sqrt(unorm_vec(i)*vnorm_vec(i));
    v.col(i) = vnorm.col(i) * sqrt(unorm_vec(i)*vnorm_vec(i));
  }
  result["u"] = u;
  result["v"] = v;
  return(result);
}



List fct_c_opt_adam(arma::mat gradients,
                    List state){
  if (state.size() == 0){
    state["beta1"] = 0.9;
    state["beta2"] = 0.999;
    state["epsilon"] = 1e-8;
    state["iteration"] = 1;
    state["m"] = 0*gradients;
    state["v"] = 0*gradients;
    state["alpha"] = 1e-2;
  }

  double beta1 = as<double>(state["beta1"]);
  double beta2 = as<double>(state["beta2"]);
  double epsilon = as<double>(state["epsilon"]);
  int iteration = as<int>(state["iteration"]);
  arma::mat m = as<arma::mat>(state["m"]);
  arma::mat v = as<arma::mat>(state["v"]);
  double alpha = as<double>(state["alpha"]);

  m = beta1 * m + (1.0 - beta1)*gradients;
  v = beta2 * v + (1.0 - beta2)*(gradients % gradients);
  arma::mat mhat = m / (1.0 - pow(beta1, iteration));
  arma::mat vhat = v / (1.0 - pow(beta2, iteration));
  state["updates"] = alpha * mhat / (sqrt(vhat) + epsilon);

  state["iteration"] = iteration + 1;

  state["m"] = m;
  state["v"] = v;

  return(state);
}


double pdist(arma::mat x, arma::sp_mat w, arma::mat index){
  double result = 0;
  arma::mat diff;
  index = index - 1;
  for (int i = 0; i < index.n_rows; i++){
    for (int j = 0; j < index.n_cols; j++){
      diff = x.col(i) - x.col(index(i,j));
      result = result + w(i,index(i,j)) * sqrt(accu(diff % diff));
    }
  }
  return(result);
}

// [[Rcpp::export]]
List fct_c_optimize(arma::sp_mat x,
                    arma::mat u,
                    arma::mat v,
                    arma::sp_mat w,
                    arma::mat index,
                    double lambda,
                    double epsilon,
                    int maxiter){

  bool converged = false;
  double objective = -std::numeric_limits<double>::max();
  double objective_prev;
  double lik;
  double diff;
  double osp;
  double u_diff;
  double v_diff;
  List results;
  List u_state;
  List v_state;
  arma::mat u_prev;
  arma::mat v_prev;
  arma::mat uv_t;
  arma::mat uv_exp;
  arma::mat u_gradient;
  arma::mat v_gradient;
  arma::mat lik_store(maxiter, 1, fill::zeros);
  arma::mat obj_store(maxiter, 1, fill::zeros);
  arma::mat udiff_store(maxiter, 1, fill::zeros);
  arma::mat vdiff_store(maxiter, 1, fill::zeros);
  arma::mat osp_store(maxiter, 1, fill::zeros);
  List u_grad_desc;
  List v_grad_desc;
  List norms;
  arma::mat u_update;
  arma::mat v_update;
  arma::mat j;
  j.ones(u.n_cols, x.n_cols);
  int i = 1;

  uv_t = u * v.t();
  uv_exp = exp(uv_t);

  while (!converged) {

    u_prev = u;
    v_prev = v;
    objective_prev = objective;

    if (i % 1 == 0){
      Rcpp::checkUserInterrupt();
    }

    u_gradient = (x - uv_exp) * v;
    v_gradient = ((x - uv_exp).t() * u)  - lambda*((2 * (w + w.t()) * j.t()) % v - (w + w.t()) * v);


    u_grad_desc = fct_c_opt_adam(u_gradient, u_state);
    v_grad_desc = fct_c_opt_adam(v_gradient, v_state);

    u_state = u_grad_desc;
    v_state = v_grad_desc;

    u_update = as<arma::mat>(u_grad_desc["updates"]);
    v_update = as<arma::mat>(v_grad_desc["updates"]);

    u = u_prev + u_update;
    v = v_prev + v_update;

    norms = uv_norm(u, v);
    u = as<arma::mat>(norms["u"]);
    v = as<arma::mat>(norms["v"]);

    uv_t = u * v.t();
    uv_exp = exp(uv_t);

    lik = accu(uv_t % x) - accu(uv_exp) ;

    osp = lambda*pdist(v.t(), w, index);
    objective = lik-osp;
    diff = abs(objective - objective_prev)/abs(objective_prev);
    u_diff = accu(abs(u_prev - u))/accu(abs(u_prev));
    v_diff = accu(abs(v_prev - v))/accu(abs(v_prev));

    lik_store[i] = lik;
    osp_store[i] = osp;
    obj_store[i] = objective;
    udiff_store[i] = u_diff;
    vdiff_store[i] = v_diff;

    if (i % 1 == 0){
      Rcout << "iteration: " << i << " | convergence:" << u_diff << " | " << v_diff << " | " << diff << "\n";
    }

    if((i >= maxiter) || ((diff < epsilon) && (u_diff < epsilon) && (v_diff < epsilon))){
      converged = true;
    }

    i = i + 1;

  }
  results["u"] = u;
  results["v"] = v;
  results["lik_penal"] = objective;
  results["likelihood"] = lik_store;
  results["osp"] = osp_store;
  results["objective"] = obj_store;
  results["udiff"] = udiff_store;
  results["vdiff"] = vdiff_store;
  return(results);

}

// [[Rcpp::export]]
arma::rowvec obs_log_like(arma::sp_mat test_x,
                          arma::mat u,
                          arma::mat v,
                          arma::mat test_nn,
                          arma::mat index_map){
  arma::rowvec cv_log_like(test_nn.n_rows, fill::zeros);
  double cv_log_like_i;
  double param_index;
  arma::uvec uv_index_loc;
  int uv_index;
  arma::mat uv_t;
  arma::mat uv_exp;
  test_nn = test_nn - 1;
  index_map = index_map - 1;
  uv_t = u * v.t();
  uv_exp = exp(uv_t);

  for (int i = 0; i < test_nn.n_rows; i++){
    for (int j = 0; j < test_nn.n_cols; j++){
      param_index = test_nn(i,j);
      uv_index_loc = find(index_map == param_index);
      uv_index = uv_index_loc(0,0);
      cv_log_like_i = accu((test_x.col(i) % uv_t.col(uv_index)) - uv_exp.col(uv_index));
      cv_log_like(i) = cv_log_like_i;
    }

  }
  return(cv_log_like);
}

