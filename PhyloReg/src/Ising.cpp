// [[Rcpp::depends(RcppArmadillo)]]
//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;



void Gibbsonce(arma::mat & state,
               const arma::vec & row_thr,
               const arma::vec & col_thr,
               const arma::mat & row_graph,
               const arma::mat & col_graph,
               int nrow,int ncol){
    arma::mat randmat(size(state), fill::randu);
    double pPlus = 0;
    for(int i = 0 ; i < nrow ; i++){
        for(int j = 0 ; j < ncol ; j++){
            pPlus = as_scalar( row_thr(i) + col_thr(j) + state.row(i) * col_graph.col(j) + 
                    row_graph.row(i) * state.col(j));
            state(i,j) = randmat(i,j) < exp(pPlus)/(exp(pPlus)+exp(-pPlus)) ? 1 : -1;
        }
    }
}

// [[Rcpp::export]]
arma::mat bootstrap_helper(const arma::vec &beta_host,
                      const arma::vec &beta_para,
                      const double &gamma_host,
                      const double &gamma_para,
                      arma::mat host_para_mat, // row for host and columns for parasites
                      const arma::mat &dist_host,
                      const arma::mat &dist_para,
                      const arma::mat &design_host,
                      const arma::mat &design_para,
                      int n_burn_in, int n_iter, int thin_by){



    arma::vec thr_host = design_host * beta_host;
    arma::vec thr_para = design_para * beta_para;
    arma::mat graph_host = dist_host * gamma_host;
    arma::mat graph_para = dist_para * gamma_para;
    int nrow = host_para_mat.n_rows;
    int ncol = host_para_mat.n_cols;
    int n_save = floor(n_iter/thin_by);
    arma::mat res(n_save, nrow*ncol, fill::zeros);
    int i_save = 0;

    for(int i = 0 ; i < (n_burn_in+n_iter) ; i++){
        Gibbsonce(host_para_mat,thr_host,thr_para,graph_host,graph_para,nrow,ncol);
        if( (i-n_burn_in)>=0 && (i+1-n_burn_in)%thin_by ==0 ){
            res.row(i_save) = trans(vectorise(host_para_mat));
            ++i_save;
        }
    }
    return(res);
}
