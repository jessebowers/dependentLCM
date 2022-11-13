// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "troubleshoot.h"

// only compile if troubleshooting...
#if TROUBLESHOOT > 0
#include <chrono>

/*****************************************************
 ****** Print Functions
 *****************************************************/
// Print progress when TROUBLESHOOT==1

//' @name trouble_start
//' @title trouble_start
//' @description When troubleshooting, run this at start of all functions. Tracks progress.
//' @param function_name name of the function in question
//' @returns An ID identifying this specific execution of this specific function
//' @keywords internal
unsigned long long int trouble_start_print(std::string function_name) {
  // Print progress (good for troubleshooting fatal errors when paired with R::sink)
  Rcpp::Rcout << function_name << " START"<< "\n";
  return 0;
}

//' @name trouble_end
//' @title trouble_end
//' @description When troubleshooting, run this at end of all functions. Concludes tracking of progress
//' @param trouble_id ID identifying the executing of the current function taken from trouble_start
//' @param function_name name of the current function
//' @keywords internal
void trouble_end_print(unsigned long long int trouble_id, std::string function_name) {
  // Print progress (good for troubleshooting fatal errors when paired with R::sink)
  Rcpp::Rcout << function_name << " END"<< "\n";
}

//' @name trouble_init
//' @title trouble_init
//' @description When troubleshooting, initializes globals
//' @keywords internal
void trouble_init_print() {
  Rcpp::Rcout << "TROUBLESHOOT==" << TROUBLESHOOT << ", printing... (recommend printing to file e.g. with 'sink(...)')" << "\n";
}

//' @name trouble_list
//' @title trouble_list
//' @description Used to output troubleshooting information
//' @keywords internal
Rcpp::List trouble_list_print() {
  return Rcpp::List::create();
}

/*****************************************************
 ****** Runtime Functions
 *****************************************************/
// Saves runtime when TROUBLESHOOT==2

typedef std::chrono::time_point<std::chrono::high_resolution_clock> ttime;

// global constants
std::vector<std::string> TROUBLE_FUNCTION_NAMES = {"colMax", "rDirichlet", "rCategorical", "count_unique", "lbeta", "which", "count_integers", "map_get", "minimum", "id2pattern", "insertSorted", "mmult", "equal_to_adjmat", "helper_compare_adjmat", "adjmat_to_equal", "product", "powl", "Hyperparameter::set_hparams #V1", "Hyperparameter::set_hparams #V2", "Hyperparameter::set_dataInfo", "Hyperparameter::print", "DomainCount::set_initial #V1", "DomainCount::set_pattern2id_map", "DomainCount::set_initial #V2", "DomainCount::reduce_items", "DomainCount::drop_item", "DomainCount::pattern2id", "DomainCount::get_ltheta", "DomainCount::id2pattern #V1", "DomainCount::id2pattern #V2", "DomainCount::countAdd", "DomainCount::list2domains", "DomainCount::copy", "DomainCount::itemsid_calc", "DomainCount::theta_alpha", "DomainCount::getloglik_marginal", "DomainCount::print", "BayesParameter::set_initial #V1", "BayesParameter::set_initial #V2", "BayesParameter::class_lprob #V1", "BayesParameter::class_lprob #V2", "BayesParameter::set_class_loglik", "BayesParameter::domain_resetCounts #V1", "BayesParameter::domain_resetCounts #V2", "BayesParameter::domain_addCount #V1", "BayesParameter::domain_addCount #V2", "BayesParameter::domain_addCounts #V1", "BayesParameter::domain_addCounts #V2", "BayesParameter::item2domainid_calc", "BayesParameter::domain_prior", "BayesParameter::get_superdomains", "BayesParameter::is_identifiable #V1", "BayesParameter::is_identifiable #V2", "BayesParameter::domain_id_new", "BayesParameter::class_pi_args", "BayesParameter::class_pi_next", "BayesParameter::classes_next", "BayesParameter::thetas_next", "BayesParameter::domain_proposal", "BayesParameter::domain_accept", "BayesParameter::domain_next", "BayesParameter::domains_next", "Archive::set_initial", "Archive::domains2mat", "Archive::add", "BayesContainer::set_initial", "BayesContainer::run", "BayesContainer::run_init", "dependentLCM_fit_cpp", "is_identifiable"};

// global variables (will receive updates!)
unsigned long long int _trouble_id = 0;
std::map<unsigned long long int, ttime> _trouble_start_times;
Rcpp::NumericVector _trouble_runtimes = Rcpp::NumericVector::create();
Rcpp::IntegerVector _trouble_runcounts = Rcpp::IntegerVector::create();

// see trouble_start above
unsigned long long int trouble_start_time(std::string function_name) {
  // Track runtime
  _trouble_runcounts(function_name) = _trouble_runcounts(function_name) + 1;
  _trouble_id += 1;
  _trouble_start_times[_trouble_id] = std::chrono::high_resolution_clock::now();
  return _trouble_id;
}

// see trouble_end above
void trouble_end_time(unsigned long long int trouble_id, std::string function_name) {
  // Track runtime
  ttime end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = end_time - _trouble_start_times[trouble_id];
  _trouble_runtimes(function_name) = _trouble_runtimes(function_name) + ms_double.count();
  _trouble_start_times.erase(trouble_id);
}

// see trouble_init above
void trouble_init_time() {
  Rcpp::Rcout << "TROUBLESHOOT==" << TROUBLESHOOT << ", saving runtimes" << "\n";
  
  // Reset times&counts&id
  _trouble_id = 0;
  int nfuns = int(TROUBLE_FUNCTION_NAMES.size());
  _trouble_runtimes = Rcpp::NumericVector(nfuns);
  _trouble_runcounts = Rcpp::IntegerVector(nfuns);
  _trouble_runtimes.names() = TROUBLE_FUNCTION_NAMES;
  _trouble_runcounts.names() = TROUBLE_FUNCTION_NAMES;
}

// see trouble_list above
Rcpp::List trouble_list_time() {
  return Rcpp::List::create(
    Rcpp::Named("runtimes") = _trouble_runtimes
  , Rcpp::Named("runcounts") = _trouble_runcounts
  );
}

#endif
