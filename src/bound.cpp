#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector optim_for_each_cpp(NumericVector x, Function f, int max_iteration, List paralist) {
  // List opt_results(max_iteration);
  // paralist[0] = para;
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];
  List path;
  for(int i = 1; i< max_iteration; ++i){
    // Rcpp::List opt_results = optim(Rcpp::_["par"]  = init_val,
    //                                Rcpp::_["fn"]     = Rcpp::InternalFunction(&verossimilhanca),
    //                                Rcpp::_["method"] = "BFGS", x)["par"];
    path = optim(Rcpp::_["par"]  = paralist[i-1],
                        Rcpp::_["fn"]     = f,
                        Rcpp::_["x"] = x);
    paralist[i] = path["par"];

  }
  List out = optim(Rcpp::_["par"]  = paralist[max_iteration-1],
                     Rcpp::_["fn"]     = f,
                     Rcpp::_["x"] = x);
  return out["value"];

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
