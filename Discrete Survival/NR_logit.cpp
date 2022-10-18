#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
//[[Rcpp::depends(RcppEigen)]]


using namespace Rcpp;
using Eigen::MatrixXd;
using namespace std;

// all the input must be sorted according to t, large->small.  
  
// [[Rcpp::export]]
List Update_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, int max_t, int c, int r){
  Eigen::MatrixXd score_t(max_t,1);
  score_t.setZero(max_t,1);
  Eigen::MatrixXd score_v(c,1);
  score_v.setZero(c,1);
  Eigen::MatrixXd info_t(max_t,1);
  info_t.setZero(max_t,1);
  Eigen::MatrixXd info_v(c,c);
  info_v.setZero(c,c);
  Eigen::MatrixXd info_tv(max_t,c);
  info_tv.setZero(max_t,c);

  for (int i = 0 ; i < r ; i++){ // i: which observation; r: number of observations
    for (int s = 1 ; s <= t(i,0) ; s++){ // s: time points of the observation i; maximum value is t(i,0); t is ordered decreasingly
        double lambda=1/(1+exp(-beta_t(s-1,0)-(X.row(i)*beta_v)(0,0)));   //lambda: logit(eta)
        score_t(s-1,0)=score_t(s-1,0)-lambda; //for ith observation, we add all time points together such that we can compute one row of score function
        score_v=score_v-lambda*X.row(i).transpose(); 
        info_t(s-1,0)=info_t(s-1,0)+lambda*(1-lambda);
        info_v=info_v+(lambda*(1-lambda))*(X.row(i).transpose()*X.row(i));
        info_tv.row(s-1)=info_tv.row(s-1)+(lambda*(1-lambda))*X.row(i);
        if (t(i,0) == s and ind(i,0) == 1){  //如果到达了当前这个人的最后一个时间，并且其在当前failure，则score要+1
          score_t(s-1,0) = score_t(s-1,0)+1;
          score_v = score_v + X.row(i).transpose();
        }    
    }
  }
  
  List result;
  result["score_t"]=score_t;
  result["score_v"]=score_v;
  result["info_t"]=info_t;
  result["info_tv"]=info_tv;
  result["info_v"]=info_v;
  
  return result;
}


// [[Rcpp::export]]
List NR_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, double tol, int max_iter){
  int r = X.rows();
  int c = X.cols();
  int max_t=t.maxCoeff();
  Eigen::MatrixXd beta2_t=beta_t;
  Eigen::MatrixXd beta2_v=beta_v;
  Eigen::MatrixXd beta_change(c,1);
  Eigen::MatrixXd A_inv(max_t,max_t);
  A_inv.setZero(max_t,max_t);
  Eigen::MatrixXd A(max_t,1);
  A.setZero(max_t,1);
  Eigen::MatrixXd B(max_t,c);
  B.setZero(max_t,c);
  Eigen::MatrixXd C(c,c);
  C.setZero(c,c);
  Eigen::MatrixXd schur(c,c);
  schur.setZero(c,c);
  Eigen::MatrixXd score_t(max_t,1);
  score_t.setZero(max_t,1);
  Eigen::MatrixXd score_v(c,1);
  score_v.setZero(c,1);

  List result;
  List update=Update_logit(t, X, ind, beta_t, beta_v, max_t, c, r);

  for (int i = 0 ; i <= max_iter; i++){
    double MaxChange_beta = 0; 
    A=update["info_t"];
    B=update["info_tv"];
    C=update["info_v"];
    score_t=update["score_t"];
    score_v=update["score_v"];
    A=A.array().inverse();
    A_inv=A.asDiagonal();
    schur=(C-(B.transpose())*(A_inv)*(B)).inverse();
    beta2_t=beta_t+(A_inv+A_inv*B*schur*(B.transpose())*A_inv)*score_t-A_inv*B*schur*score_v;
    beta2_v=beta_v+schur*score_v-schur*(B.transpose())*A_inv*score_t;
    beta_change = beta2_v - beta_v;
    update=Update_logit(t,X,ind,beta2_t,beta2_v,max_t,c,r);  
    result["beta_t"]=beta2_t;
    result["beta_v"]=beta2_v;
    result["iter"]=i;

    for (i = 0; i < c; i++){
      if (fabs(beta_change(i)) > MaxChange_beta) {
          MaxChange_beta = fabs(beta_change(i));
      }
    }
    if(MaxChange_beta < tol){ //改变converge criteria
      return result;
    }
    //cout << MaxChange_beta << endl;
    beta_t=beta2_t;
    beta_v=beta2_v;
  }

  return result;
}

