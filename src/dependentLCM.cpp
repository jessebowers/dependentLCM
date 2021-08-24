#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

/*
 * OPEN ITEMS:
 * Enforce identifiability
 * Look at vector.push and note that this always creates a whole new vector which is slow. Though resizing upfront means only one copy
 * confirm proposal probability especially when both domains are of size 1 in BayesParameter::domain_proposal(.)
 */

/*
 * FUTURE ENHANCEMENTS:
 * Allow warmup without archiving
 * Incorporate optional checks/consistency checks eg domain size, nclasses, etc
 * Look at vector.push and note that this always creates a whole new vector which is slow. Though resizing upfront means only one copy
 * post-hoc handle label switching
 * Move counts from BayesParameter::domain_accept(.) to BayesParameter::domain_proposal(.) for simplicity / more clearly defined roles
 */


/*****************************************************
 ****** UTILITIES
 *****************************************************/

int TROUBLESHOOT = 0;

//' @name colMax
//' @title colMax
//' @description Get the max of each column column of matrix
//' @keywords internal
Rcpp::IntegerVector colMax(const Rcpp::IntegerMatrix& x) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "colMax" << "\n";}
  Rcpp::IntegerVector max = x(0, Rcpp::_);
  
  for (int irow=1; irow < x.nrow(); irow++) {
    for (int jcol=0; jcol < x.ncol(); jcol++) {
      if (x(irow, jcol) > max(jcol)) {
        max(jcol) = x(irow, jcol);
      }
    }
  }
  
  return(max);
}

//' @name rDirichlet
//' @title rDirichlet
//' @description Generate random values from dirichlet distribution
//' @param deltas vector of dirichlet concentration parameters
//' @keywords internal
Rcpp::NumericVector rDirichlet(const Rcpp::NumericVector& deltas) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "rDirichlet" << "\n";}
  int C = deltas.size();
  Rcpp::NumericVector Xgamma(C);
  
  // generating gamma(deltac,1)
  for (int c = 0; c < C; c++) {
    Xgamma(c) = R::rgamma(deltas(c), 1.0);
    //Xgamma(c) = Rcpp::rgamma(1, deltas(c), scale = 1.0);
  }
  return Xgamma / sum(Xgamma);
}

//' @name rCategorical
//' @title rCategorical
//' @description Generate random values from a polytomous categorical distribution
//' @param probs Vector of probabilities of each category from 0 to probs.size()-1. Should sum to 1.
//' @keywords internal
int rCategorical(const Rcpp::NumericVector& probs) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "rCategorical" << "\n";}
  int n = probs.size();
  float unif = R::runif(0, 1);
  
  float cutoff = 0;
  for (int i = 0; i < n; i++) {
    cutoff += probs(i);
    if (unif < cutoff) {
      return i;
    }
  }
  
  return probs.size()-1; // unif~1, or sum(probs)<<1
}

//' @name count_unique
//' @title count_unique
//' @description Count number of unique values in vector
//' @keywords internal
int count_unique(const Rcpp::IntegerVector& x) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "count_unique" << "\n";}
  std::unordered_set<int> xset(x.begin(), x.end());
  return xset.size();
}

//' @name lbeta
//' @title lbeta
//' @description Calculate the beta function on log scale
//' Log(Beta(alphas)) = Log([product Gamma(alpha_i)] / Gamma(sum(alphas)))
//' @keywords internal
float lbeta(const Rcpp::NumericVector& alpha) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "lbeta" << "\n";}
  float log_gamma_total = std::lgamma(Rcpp::sum(alpha));
  float log_gammas = Rcpp::sum(Rcpp::lgamma(alpha));
  
  return (log_gammas - log_gamma_total);
}

//' @name beta
//' @title beta
//' @description Calculate the beta function
//' Beta(alphas) = [product Gamma(alpha_i)] / Gamma(sum(alphas))
//' @keywords internal
float beta(const Rcpp::NumericVector& alpha) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "beta" << "\n";}
  return std::exp(lbeta(alpha));
}

//' @name which
//' @title which
//' @description Give the (integer) indices where vector is true
//' @keywords internal
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "which" << "\n";}
  int n = x.size();
  std::list<int> out; // linked list for fast append
  
  for(int i = 0; i < n; i++) {
    if (x[i]) { // If x is true
      out.push_back(i);
    }
  }
  
  return Rcpp::wrap(out);
}


//' @name count_integers
//' @title count_integers
//' @description For each unique value of x, count the number of times that value appears
//' @keywords internal
std::map<int,  int> count_integers(const Rcpp::IntegerVector& x) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "count_integers" << "\n";}
  std::map<int,  int> counts_map;
  int i;
  int nx = x.size();
  int ix;
  for (i=0; i < nx; i++) {
    ix = x[i];
    counts_map.insert({ix, 0});
    counts_map[ix] += 1;
  }
  
  // Rcpp::IntegerMatrix counts_mat = Rcpp::IntegerMatrix(counts_map.size(), 2);
  // std::map<int,  int>::iterator iter;
  // std::map<int,  int>::const_iterator iter_end = counts_map.end();
  // i = 0;
  // for (iter = counts_map.begin(); iter!=iter_end; ++iter) {
  //   counts_mat(i, 0) = iter->first;
  //   counts_mat(i, 1) = iter->second;
  //   i += 1;
  // }
  
  return counts_map;
}

//' @name map_get
//' @title map_get
//' @description Look up a value from a map. If no value found return default
//' @param map Map of key:value pairs
//' @param key Value you wish to look up
//' @param defaultvalue What we should return if key is missing from map
//' @keywords internal
template <typename K, typename V>
V map_get(const  std::map <K,V> & map, const K & key, const V & defaultvalue ) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "map_get" << "\n";}
  
  typename std::map<K,V>::const_iterator iter = map.find( key );
  if ( iter == map.end() ) {
    return defaultvalue;
  }
  else {
    return iter->second;
  }
}

//' @name minimum
//' @title minimum
//' @description Calculate the minimum of two values
//' @keywords internal
template <typename T>
T minimum(const T x1, const T x2) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "minimum" << "\n";}
  if (x1 < x2) {
    return x1;
  } else {
    return x2;
  }
}

//' @name id2pattern
//' @title id2pattern
//' @description Convert pattern id to pattern vector
//' See other instance of id2pattern(.) for details
//' @keywords internal
Rcpp::IntegerVector id2pattern(int xpattern, const Rcpp::IntegerVector& mapvec) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "id2pattern" << "\n";}
  int nmapvec = mapvec.size();
  Rcpp::IntegerVector unmapped_vec = Rcpp::IntegerVector(nmapvec);
  
  int idivisor;
  // int iquotient;
  for (int i = 0; i < nmapvec; i++) {
    idivisor = mapvec[i];
    if (idivisor < 2) {
      // Nothing to divide
      unmapped_vec[i] = -1;
      continue;
    }
    unmapped_vec[i] = xpattern % idivisor;
    xpattern = (int)xpattern / idivisor; // Compiler should know not to recalculate?
  }
  return unmapped_vec;
}


/*****************************************************
 ****** Hyperparameters
 *****************************************************/

// Fixed settings for our bayes modeling
class Hyperparameter {
public:
  // Hyperparameters [FIXED]
  int ndomains;
  int nclass;
  Rcpp::IntegerVector class2domain;
  Rcpp::NumericVector classPi_alpha;
  float domain_alpha;
  float domain_proposal_empty;
  float domain_proposal_swap;
  int domain_nproposals;
  int domain_maxitems;
  float theta_alpha;
  Rcpp::LogicalVector steps_active;
  // Data Info
  Rcpp::IntegerVector item_nlevels;
  int nobs;
  // Inferred. Saved for speed
  int nitem;
  int nclass2domain;
  float domain_alpha_one;
  
public:
  int nclass2domain_calc() {return count_unique(class2domain);};
  int nitem_calc() {return item_nlevels.size();};
  void set_hparams(int ndomains_in, int nclass_in, const Rcpp::IntegerVector& class2domain_in, const Rcpp::NumericVector& classPi_alpha_in, float domain_alpha_in, int domain_maxitems_in, float theta_alpha_in, float domain_proposal_empty_in, float domain_proposal_swap_in, int domain_nproposals_in, Rcpp::LogicalVector steps_active_in);
  void set_hparams(Rcpp::List hparams_in);
  void set_dataInfo(const Rcpp::IntegerMatrix& x);
  void print();
};

//' @name Hyperparameter::set_hparams
//' @title Hyperparameter::set_hparams
//' @description Set hyperparameter values
//' For details see dependentLCM.r >> getStart_hparams(.)
//' @keywords internal
void Hyperparameter::set_hparams(
    int ndomains_in
  , int nclass_in
  , const Rcpp::IntegerVector& class2domain_in
  , const Rcpp::NumericVector& classPi_alpha_in
  , float domain_alpha_in
  , int domain_maxitems_in
  , float theta_alpha_in
  , float domain_proposal_empty_in
  , float domain_proposal_swap_in
  , int domain_nproposals_in
  , Rcpp::LogicalVector steps_active_in) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "Hyperparameter::set_hparams" << "\n";}
  ndomains = ndomains_in;
  nclass = nclass_in;
  class2domain = class2domain_in;
  classPi_alpha = classPi_alpha_in;
  domain_alpha = domain_alpha_in;
  domain_maxitems = domain_maxitems_in;
  theta_alpha = theta_alpha_in;
  domain_proposal_empty = domain_proposal_empty_in;
  domain_proposal_swap = domain_proposal_swap_in;
  domain_nproposals = domain_nproposals_in;
  steps_active =  steps_active_in;
  
  // Inferred
  nclass2domain = nclass2domain_calc();
}

//' @name Hyperparameter::set_hparams
//' @title Hyperparameter::set_hparams
//' @description Set hyperparameter values
//' Same as other set_hparams but with list compatability for R
//' @param hparams_in list containing all arguments for Hyperparameter::set_hparams
//' @keywords internal
void Hyperparameter::set_hparams(Rcpp::List hparams_in) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "Hyperparameter::set_hparams" << "\n";}
  int ndomains = hparams_in("ndomains");
  int nclass = hparams_in("nclass");
  Rcpp::IntegerVector class2domain = hparams_in("class2domain");
  Rcpp::NumericVector classPi_alpha = hparams_in("classPi_alpha");
  float domain_alpha = hparams_in("domain_alpha");
  int domain_maxitems = hparams_in("domain_maxitems");
  float theta_alpha = hparams_in("theta_alpha");
  float domain_proposal_empty = hparams_in("domain_proposal_empty");
  float domain_proposal_swap = hparams_in("domain_proposal_swap");
  float domain_nproposals = hparams_in("domain_nproposals");
  Rcpp::LogicalVector steps_active = hparams_in("steps_active");
  
  set_hparams(ndomains, nclass, class2domain, classPi_alpha, domain_alpha, domain_maxitems, theta_alpha, domain_proposal_empty, domain_proposal_swap, domain_nproposals, steps_active);
}

//' @name Hyperparameter::set_dataInfo
//' @title Hyperparameter::set_dataInfo
//' @description Use the raw data we are modeling to set certain hparams settings (e.g. set number of items)
//' Assumptions: That there are no empty levels especially at end
//' @keywords internal
void Hyperparameter::set_dataInfo(const Rcpp::IntegerMatrix& x) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "Hyperparameter::set_dataInfo" << "\n";}
  item_nlevels = colMax(x) + 1; // length(0:n) = n+1
  nobs = x.nrow();
  nitem = nitem_calc();
}

//' @name Hyperparameter::print
//' @title Hyperparameter::print
//' @description Print Hyperparmeter (used mainly for troubleshooting)
//' @keywords internal
void Hyperparameter::print() {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "Hyperparameter::print" << "\n";}
  Rcpp::Rcout << "hparams.ndomains:" << ndomains << "\n";
  Rcpp::Rcout << "hparams.nclass:" << nclass << "\n";
  Rcpp::Rcout << "hparams.class2domain:" << class2domain << "\n";
  Rcpp::Rcout << "hparams.classPi_alpha:" << classPi_alpha << "\n";
  Rcpp::Rcout << "hparams.domain_alpha:" << domain_alpha << "\n";
  Rcpp::Rcout << "hparams.domain_maxitems:" << domain_maxitems << "\n";
  Rcpp::Rcout << "hparams.domain_proposal_empty:" << domain_proposal_empty << "\n";
  Rcpp::Rcout << "hparams.domain_proposal_swap:" << domain_proposal_swap << "\n";
  Rcpp::Rcout << "hparams.theta_alpha:" << theta_alpha << "\n";
  Rcpp::Rcout << "hparams.item_nlevels:" << item_nlevels << "\n";
  Rcpp::Rcout << "hparams.nobs:" << nobs << "\n";
  Rcpp::Rcout << "hparams.nclass2domain:" << nclass2domain << "\n";
  Rcpp::Rcout << "hparams.nitem:" << nitem << "\n";
  Rcpp::Rcout << "hparams.steps_active:" << steps_active << "\n";
}


/*****************************************************
 ****** Domains
 *****************************************************/

// Domains and corresponding (logged) probabilities
class DomainCount {
public:
  Rcpp::NumericVector lthetas; // probabilities of corresponding pattern on log scale
  Rcpp::IntegerVector items; // items in this domain
  Rcpp::IntegerVector pattern2id_map;
  
public:
  int npatterns = 0;
  Rcpp::IntegerVector counts; // number of times each pattern appears in the data
  
public:
  void set_initial(Rcpp::IntegerVector& items_in, Hyperparameter& hparams, const Rcpp::NumericVector& lthetas_in = Rcpp::NumericVector(0));
  void set_initial(Rcpp::List list_domain, Hyperparameter& hparams);
  static std::vector<std::map<int,  DomainCount> > list2domains(Rcpp::List list_list_domains, Hyperparameter& hparams);
  void set_pattern2id_map(Hyperparameter& hparams);
  int pattern2id(Rcpp::IntegerMatrix::ConstRow xobs); // maybe switch to template
  Rcpp::IntegerVector id2pattern(int id);
  double get_ltheta(Rcpp::IntegerMatrix::ConstRow xobs); // maybe switch to template
  
public:
  void countReset();
  void countAdd(Rcpp::IntegerMatrix::ConstRow xobs);
  
public:
  int ndomainitems_calc() {return items.size();}; // number of items in this domain
  int nitems_calc() {return pattern2id_map.size();}; // number of items in the data
  DomainCount copy();
  void print();
};

//' @name DomainCount::set_initial
//' @title DomainCount::set_initial
//' @description Set all (non-count) values for this domain
//' @param items_in Which items are in this domain?
//' @param hyparams Hyperparmeters
//' @param lthetas_in Log probablities of each response pattern of these items (optional)
//' @keywords internal
void DomainCount::set_initial(Rcpp::IntegerVector& items_in, Hyperparameter& hparams, const Rcpp::NumericVector& lthetas_in) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::set_initial" << "\n";}
  items = items_in;
  
  set_pattern2id_map(hparams);
  
  
  if (lthetas_in.size() > 0) {
    lthetas = lthetas_in;
  } else {
    lthetas = Rcpp::NumericVector(npatterns);
  }
}

//' @name DomainCount::set_pattern2id_map
//' @title DomainCount::set_pattern2id_map
//' @description Create 'conversion vector' for converting (vector) response pattern to representative ID
//' See id2pattern(.) for more details
//' @keywords internal
void DomainCount::set_pattern2id_map(Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::set_pattern2id_map" << "\n";}
  // Side effect. Sets npatterns
  pattern2id_map = Rcpp::NumericVector(hparams.nitem);
  pattern2id_map.fill(0);
  int cumprod_current = 1;
  int iitem;
  int ndomainitems = ndomainitems_calc();
  for (int i = 0; i < ndomainitems; i++) {
    iitem = items[i];
    pattern2id_map[iitem] = cumprod_current;
    cumprod_current *= hparams.item_nlevels[iitem];
  }
  
  npatterns = cumprod_current; // Should match npatterns = lthetas.size() and product(item_nlevels[theseItems])
}

//' @name DomainCount::set_initial
//' @title DomainCount::set_initial
//' @description As other set_initial() but with list compatibilities for R
//' @param list_domain List of all arguments for set_initial()
//' @param hparams hyperparmeters
//' @keywords internal
void DomainCount::set_initial(Rcpp::List list_domain, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::set_initial" << "\n";}
  Rcpp::NumericVector lthetas_in =  Rcpp::log(list_domain["domains"]);
  Rcpp::IntegerVector items_in = list_domain["items"];
  set_initial(items_in, hparams, lthetas_in);
}

//' @name DomainCount::pattern2id
//' @title DomainCount::pattern2id
//' @description Convert (vector) response pattern to representative ID
//' See id2pattern(.) for more details
//' @keywords internal
int DomainCount::pattern2id(Rcpp::IntegerMatrix::ConstRow xobs) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::pattern2id" << "\n";}
  return Rcpp::sum(xobs * pattern2id_map);
}

//' @name DomainCount::get_ltheta
//' @title DomainCount::get_ltheta
//' @description Look up the log-probability of seeing this response pattern in this domain
//' @param xobs vector of the FULL response pattern (not just the items in this domain)
//' @keywords internal
double DomainCount::get_ltheta(Rcpp::IntegerMatrix::ConstRow xobs) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::get_ltheta" << "\n";}
  return lthetas(pattern2id(xobs));
}

//' @name DomainCount::id2pattern
//' @title DomainCount::id2pattern
//' @description Convert pattern ID to original response pattern
//' See utilities ::id2pattern(.) for more details
//' Differs from ::id2pattern(.) in that pattern2id_map is a cumulative compared to mapvec
//' Assumes items are in same order as pattern2id_map;
//' @keywords internal
Rcpp::IntegerVector DomainCount::id2pattern(int id) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "" << "\n";}
  Rcpp::IntegerVector pattern(nitems_calc(), -1);
  
  int i_item;
  int i_divisor;
  int i_value;
  int ndomainItems = ndomainitems_calc();
  for (int i=ndomainItems-1; i > -1; i--) {
    i_item = items[i];
    i_divisor = pattern2id_map[i_item];
    if (i_divisor == 0) {
      Rcpp::warning("DomainCount::id2pattern:: Zero divisor.");
      break; // should never happen. Error out gracefully.
    }
    i_value = id / i_divisor;  // Silently truncates by design
    pattern[i_item] = i_value;
    id = id - i_value * i_divisor;
  }
  
  return pattern;
}

//' @name DomainCount::countReset
//' @title DomainCount::countReset
//' @description Set counts to zero. 
//' Counts measure the number of times each pattern appears in the data.
//' @keywords internal
void DomainCount::countReset() {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::id2pattern" << "\n";}
  counts = Rcpp::IntegerVector(npatterns);
  counts.fill(0);
}

//' @name DomainCount::countAdd
//' @title DomainCount::countAdd
//' @description Add one observation to this domain's counts.
//' Counts measure the number of times each pattern appears in the data.
//' @param xobs One FULL response pattern (not just the items in this domain)
//' @keywords internal
void DomainCount::countAdd(Rcpp::IntegerMatrix::ConstRow xobs) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::countAdd" << "\n";}
  counts[pattern2id(xobs)] += 1;
}

//' @name DomainCount::list2domains
//' @title DomainCount::list2domains
//' @description Convert of list of list of r::domains to vector of map of cpp:domains
//' Used to convert R domains to CPP domains
//' @param list_list_domains
//' The outermost list has one slot for each class
//' Each class slot is a list containing individual domains (names/ids are ignored).
//' The lowest level is a list containing information about that individual domain
//' @param hparams hyperparameters
//' @keywords internal
std::vector<std::map<int,  DomainCount> > DomainCount::list2domains(Rcpp::List list_list_domains, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::list2domains" << "\n";}
  std::vector<std::map<int,  DomainCount> > domains;
  domains.resize(list_list_domains.length());
  
  for (int iclass=0; iclass < list_list_domains.length(); iclass++) {
    Rcpp::List iclass_domains = list_list_domains[iclass]; // slow
    for (int jdomain=0; jdomain < iclass_domains.length(); jdomain++) {
      DomainCount jjdomain;
      jjdomain.set_initial(iclass_domains[jdomain], hparams);
      domains[iclass].insert({jdomain, jjdomain});
    }
  };
  
  return domains;
}

//' @name DomainCount::copy
//' @title DomainCount::copy
//' @description Creates a deep copy of this domain
//' @keywords internal
DomainCount DomainCount::copy() {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::copy" << "\n";}
  DomainCount newDomain;
  newDomain.lthetas = Rcpp::clone(lthetas);
  newDomain.items = Rcpp::clone(items);
  newDomain.pattern2id_map = Rcpp::clone(pattern2id_map);
  newDomain.npatterns = npatterns;
  newDomain.counts = Rcpp::clone(counts);
  return newDomain;
}

//' @name DomainCount::print
//' @title DomainCount::print
//' @description Prints this domain (used mainly for troubleshooting purposes)
//' @keywords internal
void DomainCount::print() {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "DomainCount::print" << "\n";}
  Rcpp::Rcout << "domain.lthetas:" << lthetas << "\n";
  Rcpp::Rcout << "domain.items:" << items << "\n";
  Rcpp::Rcout << "domain.pattern2id_map:" << pattern2id_map << "\n";
  Rcpp::Rcout << "domain.npatterns:" << npatterns << "\n";
  Rcpp::Rcout << "domain.counts:" << counts << "\n";
}

// Example blank domain. CONSTANT - DO NOT MODIFY!
DomainCount BLANK_DOMAIN;


/*****************************************************
 ****** Bayes Parameters
 *****************************************************/

// Output for BayesParameter::domainProposal(.)
struct domainProposalOut {
  int swap_type; // -100 undecided, 2 split/put item into empty domain, 0 swap/exchange items between 2 domains, 1 transfer/put item from one domain into another
  int domain_id1;
  int domain_id2;
  int item1_old;
  int item2_old;
  float forwardProb;
  float backwardProb;
  DomainCount* domain_old1;
  DomainCount* domain_old2;
  DomainCount domain_new1;
  DomainCount domain_new2;
  Rcpp::IntegerVector domain_classes;
};

// Output for BayesParameter::domainAccept(.)
struct domainAcceptOut {
  std::vector<std::map<int,  DomainCount> > domains_new;
  float loglik_old;
  float loglik_new;
  float unif;
  float log_cutoff;
  bool accept;
};

// Output for BayesParameter::domainPropProb(.)
struct domainPropProbOut {
  Rcpp::NumericVector domain_probs;
  std::map<int,  int> domainitem_counts;
}; // List and tuples not working here so using struct

// Bayes parameters  (random variables)
class BayesParameter {
public:
  Rcpp::NumericVector class_pi; // Prior for classes
  Rcpp::IntegerVector classes; // Class for each observation
  std::vector<std::map<int,  DomainCount> > domains; // Domains and thetas
  
  public: // Inferred
    Rcpp::IntegerMatrix item2domainid; // Function of domains
    Rcpp::IntegerMatrix domains_accept; // Function of domains
    Rcpp::NumericMatrix class_loglik; // Log likelihood of seeing this pattern conditional on class and all parameters
    
public:
  void set_initial(Rcpp::NumericVector class_pi_in, Rcpp::IntegerVector classes_in, std::vector<std::map<int,  DomainCount> > domains_in, Hyperparameter& hparams);
  void set_initial(Rcpp::List list_bparam, Hyperparameter& hparams);
  float class_lprob(Rcpp::IntegerMatrix::ConstRow xobs, int xclass);
  Rcpp::NumericVector class_lprob(Rcpp::IntegerMatrix::ConstRow xobs);
  void set_class_loglik(const Rcpp::IntegerMatrix& x, bool reset=false);
  
public:
  void domain_resetCounts();
  void domain_addCount(Rcpp::IntegerMatrix::ConstRow xobs, int xclass);
  void domain_addCounts(const Rcpp::IntegerMatrix& x, bool reset_counts=true);
  void domain_resetCounts(std::vector<std::map<int,  DomainCount> >& domains);
  void domain_addCount(Rcpp::IntegerMatrix::ConstRow xobs, int xclass, std::vector<std::map<int,  DomainCount> >& domains);
  void domain_addCounts(const Rcpp::IntegerMatrix& x, const Rcpp::IntegerVector& classes, bool reset_counts, std::vector<std::map<int,  DomainCount> >& domains);
  
public:
  int nclass_calc() {return class_pi.size();};
  int nitems_calc() {return domains[0].begin()->second.nitems_calc();};
  Rcpp::IntegerMatrix item2domainid_calc(Hyperparameter& hparams);
  int domain_id_new(int class2domain_id, Hyperparameter& hparams);
  
public:
  Rcpp::NumericVector class_pi_args(Hyperparameter& hparams);
  void class_pi_next(Hyperparameter& hparams);
  void classes_next(const Rcpp::IntegerMatrix& x);
  void thetas_next(const Rcpp::IntegerMatrix& x, Hyperparameter& hparams);
  
public:
  static float domain_getloglik_x(const Rcpp::IntegerVector& pattern_counts, float theta_alpha);
  static float domain_getlik_domain(Hyperparameter& hparams);
  domainProposalOut domain_proposal(int class2domain_id, Hyperparameter& hparams);
  domainAcceptOut domain_accept(const Rcpp::IntegerMatrix& x, domainProposalOut& proposal, Hyperparameter& hparams);
  int domain_next(int class2domain_id, const Rcpp::IntegerMatrix& x, Hyperparameter& hparams);
  void domains_next(const Rcpp::IntegerMatrix& x, Hyperparameter& hparams);
};

//' @name BayesParameter::set_initial
//' @title BayesParameter::set_initial
//' @description Set all BayesParameter properties
//' See getStart_bayes_params(.) in dependentLCM.r for more details
//' @keywords internal
void BayesParameter::set_initial(Rcpp::NumericVector class_pi_in, Rcpp::IntegerVector classes_in, std::vector<std::map<int,  DomainCount> > domains_in, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::set_initial" << "\n";}
  class_pi = class_pi_in;
  classes = classes_in;
  domains = domains_in;
  item2domainid = item2domainid_calc(hparams);
  domains_accept = Rcpp::IntegerMatrix(hparams.domain_nproposals, hparams.nclass2domain);
  domains_accept.fill(-1);
  class_loglik = Rcpp::NumericMatrix(hparams.nclass, hparams.nobs); // initialize with all zeroes
}

//' @name BayesParameter::set_initial
//' @title BayesParameter::set_initial
//' @description Set all BayesParameter properties
//' As other set_initial(.) but supports lists for R compatability
//' @keywords internal
void BayesParameter::set_initial(Rcpp::List list_bparam, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::set_initial" << "\n";}
  Rcpp::NumericVector class_pi_in = list_bparam("class_pi");
  Rcpp::IntegerVector classes_in = list_bparam("classes");
  Rcpp::List list_domains_in = list_bparam("domains");
  std::vector<std::map<int,  DomainCount> > domains_in = DomainCount::list2domains(list_domains_in, hparams);
  set_initial(class_pi_in, classes_in, domains_in, hparams);
}

//' @name BayesParameter::class_lprob
//' @title BayesParameter::class_lprob
//' @description What is the (log) probability that this response pattern was generated by this class?
//' @param xobs The response pattern we are investigation.
//' @param xclass The class this observation is (assumed to be) in
//' @keywords internal
float BayesParameter::class_lprob(Rcpp::IntegerMatrix::ConstRow xobs, int xclass) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::class_lprob" << "\n";}
  float lprob = 1;
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end = domains[xclass].end();
  
  
  for (domain_iter = domains[xclass].begin(); domain_iter!=domain_end; ++domain_iter) {
    lprob += domain_iter->second.get_ltheta(xobs);
  }
  
  return lprob;
}

//' @name BayesParameter::class_lprob
//' @title BayesParameter::class_lprob
//' @description What is the (log) probability that this response pattern was generated by EACH class?
//' Assumes all parameters are fixed/known.
//' @param xobs The response pattern we are investigation.
//' @keywords internal
Rcpp::NumericVector BayesParameter::class_lprob(Rcpp::IntegerMatrix::ConstRow xobs) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::class_lprob" << "\n";}
  Rcpp::NumericVector lprobs = Rcpp::NumericVector(nclass_calc());
  
  for (int i=0; i < nclass_calc(); i++) {
    lprobs[i] = class_lprob(xobs, i);
  }
  return lprobs;
}

//' @name BayesParameter::set_class_loglik
//' @title BayesParameter::set_class_loglik
//' @description For each observation what is the (log) probability of seeing its response pattern under EACH class?
//' Updates corresponding property in BayesParameter instance.
//' Assumes all parameters are fixed/known.
//' @param x Matrix of responses
//' @param reset Whether to initialize the matrix dimensions before running
//' @keywords internal
void BayesParameter::set_class_loglik(const Rcpp::IntegerMatrix& x, bool reset) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::set_class_loglik" << "\n";}
  // One column per observation and row per class. Get conditional likelihood
  
  int xnrow = x.nrow();
  
  if (reset == true) {
    // Initialize matrix
    class_loglik = Rcpp::NumericMatrix(nclass_calc(), xnrow);
  }
  
  for (int i = 0; i < xnrow; i++) {
    class_loglik.column(i) = class_lprob(x.row(i));
  }
}

//' @name BayesParameter::domain_resetCounts
//' @title BayesParameter::domain_resetCounts
//' @description Reset the counts of all domains in given map
//' @keywords internal
void BayesParameter::domain_resetCounts(std::vector<std::map<int,  DomainCount> >& domains) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_resetCounts" << "\n";}
  int nclass = nclass_calc();
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  
  for (int iclass=0; iclass < nclass; iclass++) {
    domain_end = domains[iclass].end();
    for (domain_iter = domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      domain_iter->second.countReset();
    }
  }
}

//' @name BayesParameter::domain_resetCounts
//' @title BayesParameter::domain_resetCounts
//' @description Reset the counts in all of BayesParameter's domains
//' @keywords internal
void BayesParameter::domain_resetCounts() {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_resetCounts" << "\n";}
  domain_resetCounts(domains);
}

//' @name BayesParameter::domain_addCount
//' @title BayesParameter::domain_addCount
//' @description Add given observation to all the counts in all provided domains
//' @param xobs the response pattern we are counting
//' @param xlcass The class of this observation. Needed because different domains correspond to different classes.
//' @param domains map of domains we wish to add this pattern to
//' @keywords internal
void BayesParameter::domain_addCount(Rcpp::IntegerMatrix::ConstRow xobs, int xclass, std::vector<std::map<int,  DomainCount> >& domains) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_addCount" << "\n";}
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  domain_end = domains[xclass].end();
  for (domain_iter = domains[xclass].begin(); domain_iter!=domain_end; ++domain_iter) {
    domain_iter->second.countAdd(xobs);
  }
}

//' @name BayesParameter::domain_addCount
//' @title BayesParameter::domain_addCount
//' @description Add given observation to all the counts in all of this BayesParameter's domains
//' @param xobs the response pattern we are counting
//' @param xlcass The class of this observation. Needed because different domains correspond to different classes.
//' @keywords internal
void BayesParameter::domain_addCount(Rcpp::IntegerMatrix::ConstRow xobs, int xclass) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_addCount" << "\n";}
  domain_addCount(xobs, xclass, domains);
}

//' @name BayesParameter::domain_addCounts
//' @title BayesParameter::domain_addCounts
//' @description Add given observations to counts in all of this given domains
//' @param x Matrix of the response patterns
//' @param reset_counts True if we should set all counts to zero before counting x.
//' @param domains Map of domains holding counts
//' @description Assumes the classes of x correspond to the classes in BayesParameter.
//' @keywords internal
void BayesParameter::domain_addCounts(const Rcpp::IntegerMatrix& x, const Rcpp::IntegerVector& classes, bool reset_counts, std::vector<std::map<int,  DomainCount> >& domains) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_addCounts" << "\n";}
  
  if (reset_counts == true) {
    domain_resetCounts(domains);
  }
  int nobs = x.nrow();
  for (int obs_id=0; obs_id < nobs; obs_id++) {
    domain_addCount(x.row(obs_id), classes[obs_id], domains);
  }
}

//' @name BayesParameter::domain_addCounts
//' @title BayesParameter::domain_addCounts
//' @description Add given observations to counts in all domains in BayesParameter
//' @param x Matrix of the response patterns
//' @param reset_counts True if we should set all counts to zero before counting x.
//' @keywords internal
void BayesParameter::domain_addCounts(const Rcpp::IntegerMatrix& x, bool reset_counts) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_addCounts" << "\n";}
  domain_addCounts(x, classes, reset_counts, domains);
}

//' @name BayesParameter::item2domainid_calc
//' @title BayesParameter::item2domainid_calc
//' @description For each item find the id for the domain it belongs to.
//' If there are multiple class2domain, then multiple columns are provided (one for each)
//' @keywords internal
Rcpp::IntegerMatrix BayesParameter::item2domainid_calc(Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::item2domainid_calc" << "\n";}
  
  Rcpp::IntegerMatrix out = Rcpp::IntegerMatrix(hparams.nitem, hparams.nclass2domain);
  
  int iclass;
  Rcpp::IntegerVector iitems;
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  
  for (int iclass2domain=0; iclass2domain < hparams.nclass2domain; iclass2domain++) {
    iclass = which(hparams.class2domain == iclass2domain)[0];
    domain_end = domains[iclass].end();
    for (domain_iter = domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      iitems = domain_iter->second.items;
      for (int i=0; i < iitems.size(); i++) {
        out(iitems[i], iclass2domain) = domain_iter->first;
      }
    }
    
  }
  
  return out;
}

//' @name BayesParameter::domain_getloglik_x
//' @title BayesParameter::domain_getloglik_x
//' @description For a given domain, we want to know the probability of observing a series of responses.
//' We calculate this probability conditioning on parameters, but marginalizing over theta (and convert to log scale)
//' @param pattern_counts The counts for each pattern in this domain (vector).
//' @param theta_alpha Hyperparameter describing the Dirichlet concentration parameters for the theta prior.
//' @keywords internal
float BayesParameter::domain_getloglik_x(const Rcpp::IntegerVector& pattern_counts, float theta_alpha) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_getloglik_x" << "\n";}
  if (pattern_counts.size() == 0) {
    return 0; // log(1)
  }
  
  return lbeta(Rcpp::as<Rcpp::NumericVector> (pattern_counts) + theta_alpha);
}

//' @name BayesParameter::domain_getlik_domain
//' @title BayesParameter::domain_getlik_domain
//' @description Before observing data/other-parameters, how likely are we to put items into these particular domains?
//' Some choices of domains may be more likely than other based on prior.
//' @keywords internal
float BayesParameter::domain_getlik_domain(Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_getlik_domain" << "\n";}
  return 1; // Assume flat prior
}

//' @name BayesParameter::domain_id_new
//' @title BayesParameter::domain_id_new
//' @description Find an empty domain. Return that domain's ID
//' @keywords internal
int BayesParameter::domain_id_new(int class2domain_id, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_id_new" << "\n";}
  
  int domain_id = -1;
  int domain_class = which(hparams.class2domain == class2domain_id)[0];
  int nmax = minimum(hparams.ndomains, hparams.nitem+1);
  
  for (int i=0; i < nmax; i++) {
    if (domains[domain_class].count(i) == 0) {
      domain_id = i;
      break;
    }
  }
  
  if (domain_id == -1) {
    Rcpp::warning("BayesParameter::domain_id_new:: No empty domain found");
  }
  
  return domain_id;
}


/*****************************************************
 ****** MCMC
 *****************************************************/

//' @name BayesParameter::class_pi_args
//' @title BayesParameter::class_pi_args
//' @description Calculate the dirichlet parameters for the posterior of pi (pi used for classes)
//' @keywords internal
Rcpp::NumericVector BayesParameter::class_pi_args(Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::class_pi_args" << "\n";}
  
  Rcpp::NumericVector args = Rcpp::clone(hparams.classPi_alpha);  // Prevent overwrites. Maybe instead remove & from input
  
  Rcpp::IntegerVector::iterator classes_itr;
  Rcpp::IntegerVector::const_iterator classes_end = classes.end();
  for (Rcpp::IntegerVector::const_iterator classes_itr = classes.begin();
       classes_itr != classes_end;
       ++classes_itr) {
    args(*classes_itr) += 1;
  }
  
  return(args);
}

//' @name BayesParameter::class_pi_next
//' @title BayesParameter::class_pi_next
//' @description Do gibbs sampling to generate pi (pi used for classes)
//' @keywords internal
void BayesParameter::class_pi_next(Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::class_pi_next" << "\n";}
  
  Rcpp::NumericVector args = class_pi_args(hparams);
  class_pi = rDirichlet(args);
}

//' @name BayesParameter::classes_next
//' @title BayesParameter::classes_next
//' @description Do gibbs sampling to calculate the class of each observation
//' @keywords internal
void BayesParameter::classes_next(const Rcpp::IntegerMatrix& x) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::classes_next" << "\n";}
  
  set_class_loglik(x);
  
  int nrow = x.nrow();
  Rcpp::NumericVector class_args(nclass_calc());
  for (int i=0; i < nrow; i++) {
    class_args = Rcpp::log(class_pi) + class_loglik.column(i);
    class_args = class_args - Rcpp::max(class_args); // Divide by max for precision: [e^a,e^b,e^c] = e^a[1, e^(b-a), e^(c-a)]
    class_args = Rcpp::exp(class_args);
    class_args = class_args / Rcpp::sum(class_args);
    classes(i) = rCategorical(class_args);
  }
}

//' @name BayesParameter::thetas_next
//' @title BayesParameter::thetas_next
//' @description Do gibbs sampling to calculate thetas for each domain
//' @keywords internal
void BayesParameter::thetas_next(const Rcpp::IntegerMatrix& x, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::thetas_next" << "\n";}
  
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  DomainCount *idomain;
  Rcpp::NumericVector iprob;
  for (int iclass=0; iclass < hparams.nclass; iclass++) {
    domain_end = domains[iclass].end();
    for (domain_iter = domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      iprob = idomain->counts + hparams.theta_alpha;
      // iprob = iprob / Rcpp::sum(iprob);
      idomain->lthetas = Rcpp::log(rDirichlet(iprob));
    }
  }
}

//' @name BayesParameter::domain_proposal
//' @title BayesParameter::domain_proposal
//' @description Propose new domains in metropolis algorithm
//' Essentially we choose two domains at random and move items between them (oversimplification).
//' We take a nonempty domain. 1) With some probability we move one of its items to an empty domain. Otherwise we choose a second nonempty domain. 2) With some probability we swap a pair of items between the two nonempty domains. 3) Otherwise we move one item from the first domain into the second.
//' There are some restrictions including 1) We do not swap items if both domains are singletons (this would cause no change up to relabeling), and enforce maximum number of items per domain given in hparams.
//' @keywords internal
domainProposalOut BayesParameter::domain_proposal(int class2domain_id, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_proposal" << "\n";}
  domainProposalOut out;
  out.swap_type = -100;
  out.domain_id1 = -1;
  out.domain_id2 = -1;
  out.item1_old = -1;
  out.item2_old = -1;
  
  /*
   * Choose Domains
   */
  
  out.domain_classes = which(hparams.class2domain == class2domain_id); // maybe make fixed
  Rcpp::IntegerVector domains_nonempty = Rcpp::unique(item2domainid.column(class2domain_id));
  std::map<int, int> domainitem_counts = count_integers(item2domainid.column(class2domain_id));
  Rcpp::IntegerVector domains_chosen = Rcpp::sample(domains_nonempty, 2, false);
  
  out.domain_id1 = domains_chosen[0];
  // Decide whether to split
  if (map_get(domainitem_counts, out.domain_id1, 0) >= 2) {
    // Split possible
    if (R::runif(0,1) < hparams.domain_proposal_empty) {
      // Ok split
      out.swap_type = 2;
      out.domain_id2 = domain_id_new(class2domain_id, hparams);
    }
  }
  if (out.domain_id2 == -1) {
    // Did not split
    out.domain_id2 = domains_chosen[1];
  }
  // Save (examples of) original domains
  out.domain_old1 = &domains[out.domain_classes[0]][out.domain_id1];
  if (domains[out.domain_classes[0]].count(out.domain_id2) > 0) {
    out.domain_old2 = &domains[out.domain_classes[0]][out.domain_id2];
  } else {
    // keep empty
    out.domain_old2 = &BLANK_DOMAIN; // do not mmodify!
  }
  
  /*
   * Swap, empty, or transfer?
   */
  
  if (out.swap_type == -100) {
    // Swap not yet chosen
    if (out.domain_old2->ndomainitems_calc() >= hparams.domain_maxitems) {
      // Cannot transfer if capped. Must swap
      out.swap_type = 0;
    } else if ((out.domain_old1->ndomainitems_calc() <= 1) & (out.domain_old2->ndomainitems_calc() <= 1)) {
      // Swapping irrelevant if both domains are of size 1. Must transfer
      out.swap_type = 1;
    } else if (R::runif(0,1) < hparams.domain_proposal_swap) {
      // Swap chosen at random
      out.swap_type = 0;
    } else {
      // Transfer chosen at random
      out.swap_type = 1;
    }
  }
  
  /*
   * Choose items
   */
  
  out.item1_old = Rcpp::sample(out.domain_old1->items, 1)[0];
  if (out.swap_type == 0) {
    out.item2_old = Rcpp::sample(out.domain_old2->items, 1)[0];
  } // else item2_old = -1
  
  // Move chosen items
  Rcpp::IntegerVector items_new1 = Rcpp::clone(out.domain_old1->items);
  Rcpp::IntegerVector items_new2 = Rcpp::clone(out.domain_old2->items);
  if (out.item1_old > -1) {
    items_new1 = items_new1[items_new1 != out.item1_old];
    items_new2.push_back(out.item1_old);  //tk enforce ordering
  }
  if (out.item2_old > -1) {
    items_new2 = items_new2[items_new2 != out.item2_old];
    items_new1.push_back(out.item2_old);  //tk enforce ordering
  }
  // Save as domains
  if (items_new1.size() > 0) {
    out.domain_new1.set_initial(items_new1, hparams);
  } else {
    out.domain_new1.npatterns = 0; // Should already be defaulted but make sure
  }
  out.domain_new1.countReset(); // Should already be zero, but to make sure
  out.domain_new2.set_initial(items_new2, hparams);
  out.domain_new2.countReset();
  
  /*
   * Proposal Probabilities
   */
  
  if (out.swap_type == 0) {
    // swap. Probabilities the same (not necessarily 1, but the same)
    out.forwardProb = 1;
    out.backwardProb = 1;
  }
  
  if (out.swap_type == 1) {
    // transfer
    out.forwardProb = (
      1 / float(domains_nonempty.size()) // Choosing domain1
    * (1 - hparams.domain_proposal_empty) // No split
    * 1 / float(domains_nonempty.size()-1) // Choosing domain2
    * (1 - hparams.domain_proposal_swap) // No swap
    * 1 / float(out.domain_old1->ndomainitems_calc()) // this specific item
    );
    
    if (out.domain_old1->ndomainitems_calc() > 1) {
      // New domain1 is nonempty
      out.backwardProb = (
        1 / float(domains_nonempty.size()) // Choosing domain2 for #1
      * (1 - hparams.domain_proposal_empty) // No split
      * 1 / float(domains_nonempty.size()-1) // Choosing domain1 for #2
      * (1 - hparams.domain_proposal_swap) // No swap
      * 1 / float(out.domain_old2->ndomainitems_calc() + 1) // this specific item
      );
    } else {
      // New domain1 is empty
      out.backwardProb = (
        1 / float(domains_nonempty.size()-1) // Choosing domain2 for #1
      * hparams.domain_proposal_empty // Choose to split
      * 1 // Do not care about what specific empty domain we split into
      * 1 / float(out.domain_old2->ndomainitems_calc() + 1) // Choose this item to split off
      );
    }
  }
  if (out.swap_type == 2) {
    // split
    out.forwardProb = (
      1 / float(domains_nonempty.size()) // Choosing domain1
    * hparams.domain_proposal_empty // Choose to split
    * 1 // Do not care about what specific empty domain we split into
    * 1 / float(out.domain_old1->ndomainitems_calc()) // Choose this item to split off
    );
    
    out.backwardProb = (
      1 / float(domains_nonempty.size()+1) // Choose domain2 for #1
      * 1 / float(domains_nonempty.size()) // Choose domain1 for #2
      * (1-hparams.domain_proposal_swap) // No swap
      * 1 // Choose to transfer my only item
    );
  }
  return out;
}

//' @name BayesParameter::domain_accept
//' @title BayesParameter::domain_accept
//' @description In metroplis algorithm when choosing to update domains.
//' We examine the proposal and decide whether to accept it.
//' @keywords internal
domainAcceptOut BayesParameter::domain_accept(const Rcpp::IntegerMatrix& x, domainProposalOut& proposal, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_accept" << "\n";}
  
  domainAcceptOut out;
  
  out.domains_new.resize(hparams.nclass);
  
  int i;
  int iclass;
  for (i = 0; i < proposal.domain_classes.size(); i++) {
    iclass = proposal.domain_classes[i];
    if (proposal.domain_new1.ndomainitems_calc() > 0) {
      out.domains_new[iclass][proposal.domain_id1] = proposal.domain_new1.copy();
    }
    out.domains_new[iclass][proposal.domain_id2] = proposal.domain_new2.copy();
  }
  domain_addCounts(x, classes, true, out.domains_new); // do I really need to be true to reset?
  
  out.loglik_old = std::log(domain_getlik_domain(hparams));
  out.loglik_new = std::log(domain_getlik_domain(hparams));
  
  for (i = 0; i < proposal.domain_classes.size(); i++) {
    
    iclass = proposal.domain_classes[i];
    
    if (domains[iclass].count(proposal.domain_id1) > 0) { // should always be true
      out.loglik_old += domain_getloglik_x(domains[iclass][proposal.domain_id1].counts, hparams.theta_alpha);
    }
    if (domains[iclass].count(proposal.domain_id2) > 0) {
      out.loglik_old += domain_getloglik_x(domains[iclass][proposal.domain_id2].counts, hparams.theta_alpha);
    }
    if (out.domains_new[iclass].count(proposal.domain_id1) > 0) {
      out.loglik_new += domain_getloglik_x(out.domains_new[iclass][proposal.domain_id1].counts, hparams.theta_alpha);
    }
    if (out.domains_new[iclass].count(proposal.domain_id2) > 0) { // should always be true
      out.loglik_new += domain_getloglik_x(out.domains_new[iclass][proposal.domain_id2].counts, hparams.theta_alpha);
    }
  }
  
  out.unif = R::runif(0, 1);
  out.log_cutoff = out.loglik_new - out.loglik_old + std::log(proposal.forwardProb) - std::log(proposal.backwardProb);
  out.accept = int( out.log_cutoff > std::log(out.unif) );
  
  return out;
}

//' @name BayesParameter::domain_next
//' @title BayesParameter::domain_next
//' @description We use a metropolis algorithm to try to update domains
//' In essence we choose one item at random and either move it to another domain or swap it with an item from another domain
//' @keywords internal
int BayesParameter::domain_next(int class2domain_id, const Rcpp::IntegerMatrix& x, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domain_next" << "\n";}
  
  /***
   *** Propose
   ***/
  
  domainProposalOut proposal =  domain_proposal(class2domain_id, hparams);
  
  /***
   *** Terminate early if...
   ***/
  
  int accept = 0;
  if (proposal.domain_id1 == proposal.domain_id2) {
    // No change
    accept = 2;
    return accept;
  }
  if (proposal.domain_old1->ndomainitems_calc() + proposal.domain_old2->ndomainitems_calc() <= 1) {
    // No change up to reordering of domains. Do not bother changing
    accept = 3;
  }
  if ((proposal.domain_new1.ndomainitems_calc() > hparams.domain_maxitems)
        | (proposal.domain_new2.ndomainitems_calc() > hparams.domain_maxitems)) {
    // Over the max items per domain, reject
    accept = -2;
  } else if (false) {
    // Identifiability restrictions violated. tk need to implement
    accept = -2;
  }
  if (accept != 0) {
    return accept;
  }
  
  /***
   *** Get Likelihood
   ***/
  
  domainAcceptOut accept_info = domain_accept(x,proposal, hparams);
  accept = accept_info.accept;
  
  /***
   *** Apply changes
   ***/
  
  if (accept == 0) {
    return accept; // Change nothing
  }
  
  if (proposal.item1_old > -1) {
    item2domainid(proposal.item1_old, class2domain_id) = proposal.domain_id2;
  }
  if (proposal.item2_old > -1) {
    item2domainid(proposal.item2_old, class2domain_id) = proposal.domain_id1;
  }
  int iclass;
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  for (int i = 0; i < proposal.domain_classes.size(); i++) {
    iclass = proposal.domain_classes[i];
    domain_end = accept_info.domains_new[iclass].end();
    for (domain_iter = accept_info.domains_new[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      domains[iclass][domain_iter->first] = domain_iter->second;
    }
    if (proposal.domain_new1.ndomainitems_calc() == 0) {
      domains[iclass].erase(proposal.domain_id1);
    }
    if (proposal.domain_new2.ndomainitems_calc() == 0) {
      domains[iclass].erase(proposal.domain_id2);
    }
  }
  return accept; // Is true
}

//' @name BayesParameter::domains_next
//' @title BayesParameter::domains_next
//' @description Use metropolis algorithm to update domains. Repeat # of times set by hparams.
//' @keywords internal
void BayesParameter::domains_next(const Rcpp::IntegerMatrix& x, Hyperparameter& hparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesParameter::domains_next" << "\n";}
  
  for (int iclass2domain=0; iclass2domain < hparams.nclass2domain; iclass2domain++) {
    for (int i=0; i < hparams.domain_nproposals; i++) {
      domains_accept(i, iclass2domain) = domain_next(iclass2domain, x, hparams);
    }
  }
  
}


/*****************************************************
 ****** Archive
 *****************************************************/


// Used to record the MCMC results
class Archive {
  
public:
  Rcpp::NumericMatrix class_pi;
  Rcpp::IntegerMatrix classes;
  std::vector<Rcpp::IntegerMatrix> domains_id; // Row for iter, classid,  clique_id, pattern_id
  std::vector<Rcpp::IntegerMatrix> domains_patterns; // Row per item
  std::vector<Rcpp::NumericVector> domains_lprobs;
  std::vector<Rcpp::IntegerMatrix> domains_accept;
  std::vector<Rcpp::NumericMatrix> class_loglik;
  int next_itr;
  int maxitr;
  
public:
  void set_initial(int nclasses, int nobs, int nitems, int maxiter_in);
  void add(BayesParameter& aparams);
  void domains2mat(BayesParameter& params, int itr, Rcpp::IntegerMatrix& out_domains_id, Rcpp::IntegerMatrix& out_domains_patterns, Rcpp::NumericVector& out_domains_lprobs);
};

//' @name Archive::set_initial
//' @title Archive::set_initial
//' @description Initializes the archive (namely reserving memory)
//' @keywords internal
void Archive::set_initial(int nclasses, int nobs, int nitems, int maxiter_in) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "Archive::set_initial" << "\n";}
  next_itr = 0;
  maxitr = maxiter_in;
  
  class_pi = Rcpp::NumericMatrix(nclasses, maxiter_in);
  classes = Rcpp::IntegerMatrix(nobs, maxiter_in);
  
  domains_id.resize(maxiter_in);
  domains_patterns.resize(maxiter_in);
  domains_lprobs.resize(maxiter_in);
  domains_accept.resize(maxiter_in);
  class_loglik.resize(maxiter_in);
}

//' @name Archive::domains2mat
//' @title Archive::domains2mat
//' @description Helper function for Archive::add
//' Used to convert the domains from map form to matrix form for storage
//' Although we could make a deep copy of the domain map each time, this would be unproductive because we need it in matrix form later for R. Therfore we convert to matrix.
//' @keywords internal
void Archive::domains2mat(BayesParameter& params, int itr, Rcpp::IntegerMatrix& out_domains_id, Rcpp::IntegerMatrix& out_domains_patterns, Rcpp::NumericVector& out_domains_lprobs) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "Archive::domains2mat" << "\n";}
  
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  int iclass;
  DomainCount *idomain;
  int i_npatterns;
  int idomain_id;
  int ithis_pattern_id;
  int i;
  // misc
  int npatterns;
  int nclass = params.nclass_calc();
  int nitems = params.nitems_calc();
  
  
  // Set Size
  npatterns = 0;
  for (iclass=0; iclass < nclass; iclass++) {
    domain_end = params.domains[iclass].end();
    for (domain_iter = params.domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      npatterns += idomain->npatterns;
    }
  }
  out_domains_id = Rcpp::IntegerMatrix(4, npatterns); // Row for iter, classid,  domain_id, pattern_id
  out_domains_patterns = Rcpp::IntegerMatrix(nitems, npatterns);
  out_domains_lprobs = Rcpp::NumericVector(npatterns);
  
  
  ithis_pattern_id = 0;
  for (iclass=0; iclass < nclass; iclass++) {
    domain_end = params.domains[iclass].end();
    for (domain_iter = params.domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      i_npatterns = idomain->npatterns;
      idomain_id = domain_iter->first;
      for (i = 0; i < i_npatterns; i++) {
        // Save info
        out_domains_id(0, ithis_pattern_id) = itr;
        out_domains_id(1, ithis_pattern_id) = iclass;
        out_domains_id(2, ithis_pattern_id) = idomain_id;
        out_domains_id(3, ithis_pattern_id) = i;
        out_domains_patterns.column(ithis_pattern_id) = idomain->id2pattern(i);
        out_domains_lprobs(ithis_pattern_id) = idomain->lthetas(i);
        
        ithis_pattern_id += 1;  //Increment
      }
    }
  }
  
  
  //return std::make_tuple(out_domains_id, out_domains_patterns, out_domains_lprobs); //std::tuple<Rcpp::IntegerMatrix, Rcpp::IntegerMatrix, Rcpp::NumericVector>
}

//' @name Archive::add
//' @title Archive::add
//' @description Add latest bayes parameters to the archive
//' @keywords internal
void Archive::add(BayesParameter& aparams) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "Archive::add" << "\n";}
  
  // Iterations are added by column for speed
  // Currently assumes fixed number of iterations. For flexible iterations std::list allows for easier extension.
  
  if (next_itr >= maxitr) {
    Rcpp::warning("Archive::add:: Max storage reached");
    return; // Exit Early. Maybe in future resize or error, but not necessary now
  }
  
  class_pi.column(next_itr) = aparams.class_pi;
  classes.column(next_itr) = aparams.classes;
  Archive::domains2mat(aparams, next_itr
                         , domains_id[next_itr], domains_patterns[next_itr], domains_lprobs[next_itr]);
  domains_accept[next_itr] = Rcpp::clone(aparams.domains_accept);
  class_loglik[next_itr] = Rcpp::clone(aparams.class_loglik);
  next_itr += 1;
}


/*****************************************************
 ****** BayesContainer
 *****************************************************/

// Houses everything: hyperparameters, bayes parameters, archive, data, etc.
class BayesContainer {
public:
  Hyperparameter hparams;
  BayesParameter params;
  Rcpp::IntegerMatrix x;
  Archive archive;
public:
  void set_initial(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list, int maxitr);
  void run(int niter);
  void run_init(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list, int nitr);
};

//' @name BayesContainer::set_initial
//' @title BayesContainer::set_initial
//' @description Sets all properties of BayesContainer
//' @keywords internal
void BayesContainer::set_initial(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list, int maxitr) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesContainer::set_initial" << "\n";}
  x = x_in;
  hparams.set_hparams(hparams_list);
  hparams.set_dataInfo(x_in);
  params.set_initial(params_list, hparams);
  archive.set_initial(hparams.nclass, hparams.nobs, hparams.nitem, maxitr);
}

//' @name BayesContainer::run
//' @title BayesContainer::run
//' @description Does #nitr MCMC steps on all bayes parameters
//' @keywords internal
void BayesContainer::run(int niter) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesContainer::run" << "\n";}
  for (int i=0; i < niter; i++) {
    
    if ((hparams.steps_active["thetas"]==true)
          | (hparams.steps_active["domains"]==true)) {
      params.domain_addCounts(x, true);
    }
    if (hparams.steps_active["domains"]==true) {
      params.domains_next(x, hparams);
    }
    if (hparams.steps_active["thetas"]==true) {
      params.thetas_next(x, hparams);
    }
    if (hparams.steps_active["class_pi"]==true) {
      params.class_pi_next(hparams);
    }
    if (hparams.steps_active["classes"]==true) {
      params.classes_next(x); // always last because it recalculates likelihood
    }
    archive.add(params);
  }
}

//' @name BayesContainer::run_init
//' @title BayesContainer::run_init
//' @description First initializes, and then does #nitr MCMC steps on all bayes parameters
//' @keywords internal
void BayesContainer::run_init(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list
                                , int nitr) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "BayesContainer::run_init" << "\n";}
  set_initial(x_in, hparams_list, params_list, nitr);
  run(nitr);
}


/*****************************************************
 ****** Public Functions
 *****************************************************/

//' @name dependentLCM_fit_cpp
//' @title dependentLCM_fit_cpp
//' @description Does MCMC simulations for dependent LCM model
//' For more details see dependentLCM_fit in dependentLCM.r
// [[Rcpp::export]]
Rcpp::List dependentLCM_fit_cpp(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list
                                  , int nitr) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "dependentLCM_fit_cpp" << "\n";}
  BayesContainer bcontainer;
  bcontainer.run_init(x_in, hparams_list, params_list, nitr);
  return Rcpp::List::create(
    Rcpp::Named("class_pi") = bcontainer.archive.class_pi
    , Rcpp::Named("classes") = bcontainer.archive.classes
    , Rcpp::Named("domains_id") = wrap(bcontainer.archive.domains_id)
    , Rcpp::Named("domains_patterns") = wrap(bcontainer.archive.domains_patterns)
    , Rcpp::Named("domains_lprobs") = wrap(bcontainer.archive.domains_lprobs)
    , Rcpp::Named("next_itr") = bcontainer.archive.next_itr
    , Rcpp::Named("maxitr") = bcontainer.archive.maxitr
    , Rcpp::Named("domains_accept") = bcontainer.archive.domains_accept
    , Rcpp::Named("class_loglik") = bcontainer.archive.class_loglik
  );
}

//' @name id2pattern
//' @title id2pattern
//' @description Convert pattern id back into its original pattern vector. Then repeat each id in xpattern.
//' PatternID = Sum PatternVec[i] * [product^(i-1) mapvec[j]]
//' @param xpattern Vector of ids each representing a vector
//' @param mapvec Vector for converting id to/from vector. 
//' mapvec is of size equal to the pattern vector length.
//' mapvec values of 0 produce a PatternVec values of -1 (consider -1 to be NA).
//' @examples Suppose we have a vector of length 4 which we wish to conver to an ID (this function reverses this process).
//' Suppose the first three indexes take values 0:1 and the last takes values 0:3.
//' The appropriate mapvec in this situation would be [2,2,2,3] equal to the number of unique values in each position.
//' If we have a vector of [1,0,1,2] the corresponding id is 1*(1) + 0*(2) + 1*(2*2) + 2*(2*2*2) = 21.
//' If we have a vector id of 12 we cand find a corresponding vector of [0,0,1,1].
// [[Rcpp::export]]
Rcpp::IntegerMatrix id2pattern(const Rcpp::IntegerVector& xpattern, const Rcpp::IntegerVector& mapvec) {
  if (TROUBLESHOOT==1){Rcpp::Rcout << "id2pattern" << "\n";}
  int npatterns = xpattern.size();
  int nmapvec = mapvec.size();
  Rcpp::IntegerMatrix unmapped_mat = Rcpp::IntegerMatrix(nmapvec, npatterns);
  
  for (int i = 0; i < npatterns; i++) {
    unmapped_mat.column(i) = id2pattern(xpattern(i), mapvec);
  }
  
  return unmapped_mat; 
}

