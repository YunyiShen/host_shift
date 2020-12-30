// [[Rcpp::depends(RcppArmadillo)]]
//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

using namespace roptim;


// calculate log pseudolikelihood and gradient
class logPL
{
public:
    arma::vec beta_host_grad, beta_para_grad;
    double Pl, gamma_host_grad, gamma_para_grad;
    void calculate(const arma::vec &beta_host,
                   const arma::vec &beta_para,
                   const double &gamma_host,
                   const double &gamma_para,
                   const arma::mat &host_para_mat,
                   const arma::mat &dist_host,
                   const arma::mat &dist_para,
                   const arma::mat &design_host,
                   const arma::mat &design_para);
};

void logPL::calculate(const arma::vec &beta_host,
                      const arma::vec &beta_para,
                      const double &gamma_host,
                      const double &gamma_para,
                      const arma::mat &host_para_mat, // row for host and columns for parasites
                      const arma::mat &dist_host,
                      const arma::mat &dist_para,
                      const arma::mat &design_host,
                      const arma::mat &design_para)
{
    int n_host = host_para_mat.n_rows;
    int n_para = host_para_mat.n_cols;
    double pPlus;
    beta_host_grad = 0 * beta_host;
    beta_para_grad = 0 * beta_para;
    gamma_host_grad = 0;
    gamma_para_grad = 0;
    double dPl_dpPlus = 0;
    Pl = 0;

    for (int i = 0; i < n_host; i++)
    {
        for (int j = 0; j < n_para; j++)
        {
            pPlus = as_scalar(design_host.row(i) * beta_host +                         // predictors for hosts
                              design_para.row(j) * beta_para +                         // predictors for parasites
                              gamma_host * (dist_host.row(i) * host_para_mat.col(j)) + // phylosignal for host
                              gamma_para * (host_para_mat.row(i) * dist_para.col(j))); // phylosignal for parasite

            dPl_dpPlus = (host_para_mat(i, j) - (exp(pPlus) - exp(-pPlus)) / (exp(pPlus) + exp(-pPlus)));
            beta_host_grad += dPl_dpPlus * trans(design_host.row(i));
            beta_para_grad += dPl_dpPlus * trans(design_para.row(j));
            gamma_host_grad += as_scalar(dPl_dpPlus * (dist_host.row(i) * host_para_mat.col(j)));
            gamma_para_grad += as_scalar(dPl_dpPlus * (host_para_mat.row(i) * dist_para.col(j)));
            Pl += pPlus * host_para_mat(i, j) - log(exp(pPlus) + exp(-pPlus));
        }
    }
    return;
}



class optlogPL : public Functor
{
public:
    const arma::mat &host_para_mat; // row for host and columns for parasites
    const arma::mat &dist_host;
    const arma::mat &dist_para;
    const arma::mat &design_host;
    const arma::mat &design_para;
    const int &p_para;
    const int &p_host;
    logPL Pl;
    // set paras
    optlogPL(const arma::mat &host_para_mat,
             const arma::mat &dist_host,
             const arma::mat &dist_para,
             const arma::mat &design_host,
             const arma::mat &design_para,
             const int &p_para,
             const int &p_host) : host_para_mat(host_para_mat),
                                  dist_host(dist_host),
                                  dist_para(dist_para),
                                  design_host(design_host),
                                  design_para(design_para),
                                  p_para(p_para),
                                  p_host(p_host) {}

    double operator()(const arma::vec &par) override
    {
        arma::vec beta_host = par.subvec(0, p_host - 1);
        arma::vec beta_para = par.subvec(p_host, p_host + p_para - 1);
        double gamma_host = par(p_host + p_para);
        double gamma_para = par(p_host + p_para + 1);
        Pl.calculate(beta_host, beta_para, gamma_host, gamma_para,
                     host_para_mat, dist_host, dist_para, design_host, design_para);
        return -Pl.Pl;
    }

    void Gradient(const arma::vec &par, arma::vec &gr) override
    {
        gr = 0 * par;
        gr.subvec(0, p_host - 1) = -Pl.beta_host_grad;
        gr.subvec(p_host, p_host + p_para - 1) = -Pl.beta_para_grad;
        gr(p_host + p_para) = -Pl.gamma_host_grad;
        gr(p_host + p_para + 1) = -Pl.gamma_para_grad;
        return;
    }
};

// [[Rcpp::export]]
List maxPL_cpp(const arma::mat &host_para_mat, // row for host and columns for parasites
                 const arma::mat &dist_host,
                 const arma::mat &dist_para,
                 const arma::mat &design_host,
                 const arma::mat &design_para,
                 const std::string &method)
{
    int p_host = design_host.n_cols;
    int p_para = design_para.n_cols;

    arma::vec par(p_host + p_para + 2, fill::randn);
    optlogPL logPL1(host_para_mat, dist_host, dist_para,
                    design_host, design_para, p_para, p_host);
    Roptim<optlogPL> opt(method);

    opt.control.trace = 0;
    opt.set_hessian(false);
    opt.minimize(logPL1, par);
    List res = Rcpp::List::create(
            Rcpp::Named("beta_host") = par.subvec(0, p_host - 1),
            Rcpp::Named("beta_para") = par.subvec(p_host, p_host + p_para - 1),
            Rcpp::Named("gamma_host") = par(p_host + p_para),
            Rcpp::Named("gamma_para") = par(p_host + p_para + 1));

    return res;
    //return(par);
}
