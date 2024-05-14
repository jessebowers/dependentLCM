// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "troubleshoot.h"

// only compile if troubleshooting...
#if TROUBLESHOOT > 0

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
void trouble_start_print(const std::string& function_name) {
  // Print progress (good for troubleshooting fatal errors when paired with R::sink)
  Rcpp::Rcout << function_name << " START"<< "\n";
}

//' @name trouble_end
//' @title trouble_end
//' @description When troubleshooting, run this at end of all functions. Concludes tracking of progress
//' @param trouble_id ID identifying the executing of the current function taken from trouble_start
//' @param function_name name of the current function
//' @keywords internal
void trouble_end_print(const std::string& function_name) {
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

// global variables (will receive updates!)
std::unordered_map<std::string,double> _trouble_runtimes;
std::unordered_map<std::string,unsigned int> _trouble_runcounts;

// see trouble_start above
ttime trouble_start_time(const std::string& function_name) {
  
  // Add Execution
  std::unordered_map<std::string,unsigned int>::iterator it = _trouble_runcounts.find(function_name);
  if( it != _trouble_runcounts.end() ) {
      it->second += 1;
  }
  else {
      _trouble_runcounts.insert(std::make_pair(function_name, 1));
  }
  
  return std::chrono::high_resolution_clock::now();
}

// see trouble_end above
void trouble_end_time(ttime start_time, const std::string& function_name) {
  // Track runtime
  ttime end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = end_time - start_time;
  
  // Add Runtime
  std::unordered_map<std::string, double>::iterator it = _trouble_runtimes.find(function_name);
  if( it != _trouble_runtimes.end() ) {
      it->second += ms_double.count();
  }
  else {
      _trouble_runtimes.insert(std::make_pair(function_name, ms_double.count()));
  }
  
}

// see trouble_init above
void trouble_init_time() {
  Rcpp::Rcout << "TROUBLESHOOT==" << TROUBLESHOOT << ", saving runtimes" << "\n";
  
  // Reset times&counts&id
  _trouble_runtimes.clear();
  _trouble_runcounts.clear();
}

// see trouble_list above
Rcpp::List trouble_list_time() {
  
  Rcpp::NumericVector runtimes;
  std::unordered_map<std::string,double>::const_iterator ditr;
  for (ditr = _trouble_runtimes.begin(); ditr != _trouble_runtimes.end(); ++ditr) {
    runtimes[ditr->first] = ditr->second;
  }
  
  Rcpp::IntegerVector runcounts;
  std::unordered_map<std::string,unsigned int>::const_iterator iitr;
  for (iitr = _trouble_runcounts.begin(); iitr != _trouble_runcounts.end(); ++iitr) {
    runcounts[iitr->first] = iitr->second;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("runtimes") = runtimes
  , Rcpp::Named("runcounts") = runcounts
  );
}

#endif
