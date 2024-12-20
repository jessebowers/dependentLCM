// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "troubleshoot.h"

/*****************************************************
 ****** UTILITIES
 *****************************************************/


//' @name colMax
//' @title colMax
//' @description Get the max of each column column of matrix
//' @keywords internal
Rcpp::IntegerVector colMax(const Rcpp::IntegerMatrix& x) {
  TROUBLE_START(("colMax"));
  Rcpp::IntegerVector max = x(0, Rcpp::_);
  
  for (int irow=1; irow < x.nrow(); irow++) {
    for (int jcol=0; jcol < x.ncol(); jcol++) {
      if (x(irow, jcol) > max(jcol)) {
        max(jcol) = x(irow, jcol);
      }
    }
  }
  
  TROUBLE_END; return max;
}

//' @name sumLogs
//' @title sumLogs
//' @description Calculate log(sum(exp(x))) precisely
//' @keywords internal
template <typename T>
double sumLogs(const T& x) {
  double xmax = max(x);
  return xmax + log(sum(exp(x-xmax)));
}

//' @name rDirichlet
//' @title rDirichlet
//' @description Generate random values from dirichlet distribution
//' @param deltas vector of dirichlet concentration parameters
//' @keywords internal
Rcpp::NumericVector rDirichlet(const Rcpp::NumericVector& deltas) {
  TROUBLE_START(("rDirichlet"));
  int C = deltas.size();
  Rcpp::NumericVector Xgamma(C);
  
  // generating gamma(deltac,1)
  for (int c = 0; c < C; c++) {
    Xgamma(c) = R::rgamma(deltas(c), 1.0);
    //Xgamma(c) = Rcpp::rgamma(1, deltas(c), scale = 1.0);
  }
  TROUBLE_END; return (Xgamma / sum(Xgamma));
}

//' @name rCategorical
//' @title rCategorical
//' @description Generate random values from a polytomous categorical distribution
//' @param probs Vector of probabilities of each category from 0 to probs.size()-1. Should sum to 1.
//' @keywords internal
int rCategorical(const Rcpp::NumericVector& probs) {
  TROUBLE_START(("rCategorical"));
  int n = probs.size();
  float unif = R::runif(0, 1);
  
  float cutoff = 0;
  for (int i = 0; i < n; i++) {
    cutoff += probs(i);
    if (unif < cutoff) {
      TROUBLE_END; return i;
    }
  }
  
  TROUBLE_END; return (probs.size()-1); // unif~1, or sum(probs)<<1
}

//' @name count_unique
//' @title count_unique
//' @description Count number of unique values in vector
//' @keywords internal
int count_unique(const Rcpp::IntegerVector& x) {
  TROUBLE_START(("count_unique"));
  std::unordered_set<int> xset(x.begin(), x.end());
  TROUBLE_END; return xset.size();
}

//' @name lbeta
//' @title lbeta
//' @description Calculate the beta function on log scale
//' Log(Beta(alphas)) = Log([product Gamma(alpha_i)] / Gamma(sum(alphas)))
//' @keywords internal
float lbeta(const Rcpp::NumericVector& alpha) {
  TROUBLE_START(("lbeta"));
  float log_gamma_total = std::lgamma(Rcpp::sum(alpha));
  float log_gammas = Rcpp::sum(Rcpp::lgamma(alpha));
  
  TROUBLE_END; return (log_gammas - log_gamma_total);
}

//' @name which
//' @title which
//' @description Give the (integer) indices where vector is true
//' @keywords internal
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x) {
  TROUBLE_START(("which"));
  int n = x.size();
  std::list<int> out; // linked list for fast append
  
  for(int i = 0; i < n; i++) {
    if (x[i]) { // If x is true
      out.push_back(i);
    }
  }
  
  TROUBLE_END; return Rcpp::wrap(out);
}


//' @name count_integers
//' @title count_integers
//' @description For each unique value of x, count the number of times that value appears
//' @keywords internal
std::map<int,  int> count_integers(const Rcpp::IntegerVector& x) {
  TROUBLE_START(("count_integers"));
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
  
  TROUBLE_END; return counts_map;
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
  TROUBLE_START(("map_get"));
  
  typename std::map<K,V>::const_iterator iter = map.find( key );
  if ( iter == map.end() ) {
    TROUBLE_END; return defaultvalue;
  }
  else {
    TROUBLE_END; return iter->second;
  }
}

//' @name minimum
//' @title minimum
//' @description Calculate the minimum of two values
//' @keywords internal
template <typename T>
T minimum(const T x1, const T x2) {
  TROUBLE_START(("minimum"));
  if (x1 < x2) {
    TROUBLE_END; return x1;
  } else {
    TROUBLE_END; return x2;
  }
}

//' @name id2pattern
//' @title id2pattern
//' @description Convert pattern id to pattern vector
//' See other instance of id2pattern(.) for details
//' @keywords internal
Rcpp::IntegerVector id2pattern(int xpattern, const Rcpp::IntegerVector& mapvec) {
  TROUBLE_START(("id2pattern"));
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
  TROUBLE_END; return unmapped_vec;
}


//' @name insertSorted
//' @title insertSorted
//' @description Modifies x. Assumes x is sorted. Add new_value to x in correct position based on sort
//' @param x Sorted vector of numbers
//' @param new_value new number we wish to insert into x
//' @keywords internal
void insertSorted(Rcpp::IntegerVector& x, int new_value) {
  TROUBLE_START(("insertSorted"));
  int n = x.size();
  int i;
  for (i = 0; i < n; i++) {
    if (new_value <= x(i)) {
      break;
    }
  }
  x.insert(i, new_value);
  TROUBLE_END; return;
}

//' @name mmult
//' @title mmult
//' @description Multiply two matrixes
//' @keywords internal
Rcpp::IntegerMatrix mmult(Rcpp::IntegerMatrix& m1, Rcpp::IntegerMatrix& m2) {
  TROUBLE_START(("mmult"));
  Rcpp::Environment base("package:base");
  Rcpp::Function mat_Mult = base["%*%"]; // Steals from R::%*%.
  TROUBLE_END; return mat_Mult(m1, m2);
}

//' @name equal_to_adjmat
//' @title equal_to_adjmat
//' @description Convert equivalence classes into a (graph theory style) adjacency matrix.
//' All elements of the same equivalence class are considered neighbors.
//' @param eqclass_vec Vector describing the equivalence classes. 
//' Each index represents a separate item and indexes with the same value are in the same equivalence class.
//' @keywords internal
Rcpp::IntegerMatrix equal_to_adjmat(Rcpp::IntegerVector eqclass_vec) {
  TROUBLE_START(("equal_to_adjmat"));
  int n = eqclass_vec.size();
  Rcpp::IntegerMatrix adjmat = Rcpp::IntegerMatrix(n, n);
  
  for (int i = 0; i < n; i++) {
    adjmat.row(i) = (eqclass_vec == eqclass_vec(i));
    // Can we speed up by processing all of the same class at once?
  }
  
  TROUBLE_END; return adjmat;
}

//' @name helper_compare_adjmat
//' @title helper_compare_adjmat
//' @description Helper function. Returns true if (m1>0) == (m2>0) in all cells.
//' In other words returns true if both adjacency matrixes have the same connections (ignoring # of possible routes)
//' @keywords internal
bool helper_compare_adjmat(Rcpp::IntegerMatrix& m1, Rcpp::IntegerMatrix& m2) {
  TROUBLE_START(("helper_compare_adjmat"));
  int nrow = m1.nrow(); // m2 assumed to be of same size
  int ncol = m1.ncol(); // m2 assumed to be of same size
  
  int is_same = true;
  for (int irow=0; irow < nrow; irow++) {
    for (int icol=0; icol < ncol; icol++) {
      if ((m1(irow, icol) > 0) != (m2(irow, icol) > 0)) {
        is_same = false;
        break;
      }
    }
  }
  
  TROUBLE_END; return is_same;
}

//' @name adjmat_to_equal
//' @title adjmat_to_equal
//' @description Take a (graph theory) adjacency matrix and build equivalence classes.
//' Two nodes are in the same class if there is a path between them.
//' Returns vector. Each position represents a node and nodes with the same value are part of the same equivalence class.
//' Implicitally assumes that edges are bidirectional (a path in one direction leads to a path in the other)
//' @param adjmat Matrix identifying which nodes are neighbors (are linked).
//' @param maxitr Integer giving the maximum number of attempts to connect two nodes.
//' @keywords internal
Rcpp::IntegerVector adjmat_to_equal(Rcpp::IntegerMatrix adjmat, int maxitr = 100) {
  TROUBLE_START(("adjmat_to_equal"));
  
  int nitems = adjmat.nrow();
  
  // Find all linked nodes
  Rcpp::IntegerMatrix adjmat_new;
  for (int i=0; i < maxitr; i++) { // max just to prevent infinite loops
    adjmat_new = mmult(adjmat, adjmat);
    if (helper_compare_adjmat(adjmat, adjmat_new)==true) {
      // no changes
      break;
    }
    adjmat = adjmat_new;
  }
  
  // Convert to equivalence class vector
  Rcpp::IntegerVector equal_classes = Rcpp::IntegerVector(nitems);
  for (int iitem = 0; iitem < nitems; iitem++) {
    for (int jitem = 0; jitem < nitems; jitem++) {
      if (adjmat(iitem, jitem) > 0) {
        equal_classes(iitem) = jitem; // Should be same across equivalence class
        break;
      }
    }
  }
  
  TROUBLE_END; return equal_classes;
}

//' @name product
//' @title product
//' @description Multiply all elements of x together
//' @keywords internal
int product(Rcpp::IntegerVector x) {
  TROUBLE_START(("product"));
  int n = x.size();
  int agg = 1;
  for (int i=0; i<n; i++) {
    agg *= x[i];
  }
  TROUBLE_END; return agg;
}

//' @name powl
//' @title powl
//' @description calculate x^p power on integers
//' @keywords internal
int powl(int x, int p) {
  TROUBLE_START(("powl"));
  TROUBLE_END; return int(std::pow(double(x), double(p))+0.5);
}

//' @name is_identifiable
//' @title is_identifiable
//' @description Check if choice of domains is generically identifiable
//' Uses greedy algorithm. May give false negatives in some cases, but is quick and deterministic
//' See Allman paper (DOI:10.1214/09-AOS689 Theorem 4.) for criteria used: min(patterns1,nclass)+min(patterns2,nclass)+min(patterns3,nclass) > 2*nclass+2
//' @param item2superdomainid Vector describing which items must be grouped together
//' @param nclass Number of classes
//' @param item_nlevels Vector with levels of items
//' @keywords internal
bool is_identifiable(const Rcpp::IntegerVector& item2superdomainid, int nclass, const Rcpp::IntegerVector& item_nlevels) {
  TROUBLE_START(("is_identifiable"));
  
  int i;
  int domain_id;
  int nitems = item2superdomainid.size();
  int goal = 2*nclass + 2;
  
  // count patterns
  std::map<int, int> pattern_counts_map;
  for (i=0; i < nitems; i++) {
    domain_id = item2superdomainid(i);
    if (pattern_counts_map.count(domain_id)==0) {
      pattern_counts_map[domain_id] = product(item_nlevels[item2superdomainid==domain_id]);
      pattern_counts_map[domain_id] = minimum(pattern_counts_map[domain_id], nclass); // Beyond nclass not relevent. Don't necessarily need to do this
    }
  }
  
  // convert map to vector
  int ndomains = pattern_counts_map.size();
  Rcpp::IntegerVector pattern_counts = Rcpp::IntegerVector(ndomains);
  std::map<int,int>::iterator map_itr;
  std::map<int,int>::const_iterator map_end = pattern_counts_map.end();
  i = 0;
  for (map_itr = pattern_counts_map.begin(); map_itr != map_end; map_itr++) {
    pattern_counts[i] = map_itr->second; // Add number of patterns
    i += 1;
  }
  pattern_counts.sort(true); // descending
  
  // greedy search
  Rcpp::IntegerVector tripart = Rcpp::IntegerVector(3, 0);
  int best_index;
  int best_value;
  int best_diff;
  int new_value;
  int new_diff;
  int tripart_sum = 0;
  for (i=0; i<ndomains; i++) {
    
    // Find best group for i'th domain
    best_index = 0;
    best_value = 0;
    best_diff = 0;
    for (int j=0; j<3; j++) {
      
      if (tripart[j] > 0) {
        new_value = tripart[j] *  pattern_counts[i];
      } else {
        new_value = 0 + pattern_counts[i];
      }
      new_value = minimum(new_value, nclass);
      new_diff = new_value - tripart[j];
      
      if (new_diff > best_diff) {
        best_index = j;
        best_value = new_value;
        best_diff = new_diff;
      }
      
    }
    
    // Add to chosen group
    tripart[best_index] = best_value;
    tripart_sum += best_diff; // same as tripart.sum()
    
    if (tripart_sum >= goal) {
      break; // no need to continue
    }
    
  }
  
  TROUBLE_END; return (tripart_sum >= goal);
}

//' @name print_itr
//' @title print_itr
//' @description If itr is divisible by print_itrs then print progress to terminal
//' @param itr Number of iterations completed
//' @param nclass Number of classes.
//' @param print_itrs How frequently to print progress.
//' @keywords internal
void print_itr(const int itr, const int print_itrs) {
  TROUBLE_START(("print_itr"));
  
  if (itr % print_itrs == 0) {
    Rcpp::Rcout << "." << itr;
  }
  
  TROUBLE_END;
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
  float domain_proposal_empty;
  int domain_nproposals;
  int domain_maxitems;
  float theta_alpha;
  Rcpp::LogicalVector steps_active;
  std::string domain_theta_prior_type;
  Rcpp::NumericMatrix domainPriorKnowledgeLog_mat;
  // Data Info
  Rcpp::IntegerVector item_nlevels;
  int nobs;
  // Inferred. Saved for speed
  int nitem;
  int nclass2domain;
  int nitr;
  Rcpp::IntegerVector save_itrs;
  std::list<std::pair<int,int>> domainPriorKnowledgeLog_pairs;
  int print_itrs;
  
public:
  int nclass2domain_calc() const {return count_unique(class2domain);};
  int nitem_calc() const {return item_nlevels.size();};
  void set_hparams(int ndomains_in, int nclass_in, const Rcpp::IntegerVector& class2domain_in, const Rcpp::NumericVector& classPi_alpha_in, int domain_maxitems_in, float theta_alpha_in, float domain_proposal_empty_in, int domain_nproposals_in, Rcpp::LogicalVector steps_active_in, std::string domain_theta_prior_type_in, int nitr_in, Rcpp::IntegerVector& save_itrs_in, const Rcpp::NumericMatrix& domainPriorKnowledgeLog_mat_in, const int print_itrs_in);
  void set_hparams(Rcpp::List hparams_in);
  void set_dataInfo(const Rcpp::IntegerMatrix& x);
  static std::list<std::pair<int,int>> make_pairs(const Rcpp::NumericMatrix& domainPriorKnowledgeLog_mat_in);
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
  , int domain_maxitems_in
  , float theta_alpha_in
  , float domain_proposal_empty_in
  , int domain_nproposals_in
  , Rcpp::LogicalVector steps_active_in
  , std::string domain_theta_prior_type_in
  , int nitr_in
  , Rcpp::IntegerVector& save_itrs_in
  , const Rcpp::NumericMatrix& domainPriorKnowledgeLog_mat_in
  , const int print_itrs_in) {
  TROUBLE_START(("Hyperparameter::set_hparams #V1"));
  ndomains = ndomains_in;
  nclass = nclass_in;
  class2domain = class2domain_in;
  classPi_alpha = classPi_alpha_in;
  domain_maxitems = domain_maxitems_in;
  theta_alpha = theta_alpha_in;
  domain_proposal_empty = domain_proposal_empty_in;
  domain_nproposals = domain_nproposals_in;
  steps_active =  steps_active_in;
  domain_theta_prior_type = domain_theta_prior_type_in;
  nitr = nitr_in;
  save_itrs = save_itrs_in;
  domainPriorKnowledgeLog_mat = domainPriorKnowledgeLog_mat_in;
  print_itrs = print_itrs_in;
  
  // Inferred
  nclass2domain = nclass2domain_calc();
  domainPriorKnowledgeLog_pairs = make_pairs(domainPriorKnowledgeLog_mat)
  TROUBLE_END;
}

//' @name Hyperparameter::set_hparams
//' @title Hyperparameter::set_hparams
//' @description Set hyperparameter values
//' Same as other set_hparams but with list compatability for R
//' @param hparams_in list containing all arguments for Hyperparameter::set_hparams
//' @keywords internal
void Hyperparameter::set_hparams(Rcpp::List hparams_in) {
  TROUBLE_START(("Hyperparameter::set_hparams #V2"));
  
  int ndomains = hparams_in("ndomains");
  int nclass = hparams_in("nclass");
  Rcpp::IntegerVector class2domain = hparams_in("class2domain");
  Rcpp::NumericVector classPi_alpha = hparams_in("classPi_alpha");
  int domain_maxitems = hparams_in("domain_maxitems");
  float theta_alpha = hparams_in("theta_alpha");
  float domain_proposal_empty = hparams_in("domain_proposal_empty");
  float domain_nproposals = hparams_in("domain_nproposals");
  Rcpp::LogicalVector steps_active = hparams_in("steps_active");
  std::string domain_theta_prior_type = hparams_in("domain_theta_prior_type");
  int nitr = hparams_in("nitr");
  Rcpp::IntegerVector save_itrs = hparams_in("save_itrs");
  Rcpp::NumericMatrix domainPriorKnowledgeLog_mat = hparams_in("domainPriorKnowledgeLog");
  int print_itrs = hparams_in("print_itrs");
  
  set_hparams(ndomains, nclass, class2domain, classPi_alpha, domain_maxitems, theta_alpha, domain_proposal_empty, domain_nproposals, steps_active, domain_theta_prior_type, nitr, save_itrs, domainPriorKnowledgeLog_mat, print_itrs);
  
  TROUBLE_END;
}

//' @name Hyperparameter::set_dataInfo
//' @title Hyperparameter::set_dataInfo
//' @description Use the raw data we are modeling to set certain hparams settings (e.g. set number of items)
//' Assumptions: That there are no empty levels especially at end
//' @keywords internal
void Hyperparameter::set_dataInfo(const Rcpp::IntegerMatrix& x) {
  TROUBLE_START(("Hyperparameter::set_dataInfo"));
  item_nlevels = colMax(x) + 1; // length(0:n) = n+1
  nobs = x.nrow();
  nitem = nitem_calc();
  TROUBLE_END;
}

//' @name Hyperparameter::make_pairs
//' @title Hyperparameter::make_pairs
//' @description Find every nonzero index of domainPriorKnowledgeLog_mat_in. Only searches upper triangle
//' @return One column per nonzero value. Two rows indicating the row and column position respectively.
//' @keywords internal
 std::list<std::pair<int,int>> Hyperparameter::make_pairs(const Rcpp::NumericMatrix& domainPriorKnowledgeLog_mat_in) {
  TROUBLE_START(("Hyperparameter::set_dataInfo"));
  
  const double epsilon = 1e-8;
  
  // find nonzero pairs
  std::list<std::pair<int,int>> indices;
  for (int irow=0; irow < domainPriorKnowledgeLog_mat_in.nrow(); irow++) {
    for (int icol=irow+1; icol < domainPriorKnowledgeLog_mat_in.ncol(); icol++) {
      if (std::abs(domainPriorKnowledgeLog_mat_in(irow, icol)) > epsilon) {
        indices.push_back(std::pair<int,int>{irow, icol});
      }
    }
  }
  
  TROUBLE_END; return indices;
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
  int npatterns = 0;
  Rcpp::IntegerVector counts; // number of times each pattern appears in the data
  
public:
  void set_initial(Rcpp::IntegerVector& items_in, const Hyperparameter& hparams, const Rcpp::NumericVector& lthetas_in = Rcpp::NumericVector(0));
  void set_initial(Rcpp::List list_domain, const Hyperparameter& hparams);
  static std::vector<std::map<int,  DomainCount> > list2domains(Rcpp::List list_list_domains, const Hyperparameter& hparams);
  void set_pattern2id_map(const Rcpp::IntegerVector& item_nlevels);
  
public:
  template <typename vectype> int pattern2id(const vectype& xobs);
  Rcpp::IntegerVector id2pattern(int id);
  int id2altid(int id, const Hyperparameter& hparams);
  int itemsid_calc();
  int ndomainitems_calc() {return items.size();}; // number of items in this domain
  int nitems_calc() {return pattern2id_map.size();}; // number of items in the data
  
public:
  double get_ltheta(const Rcpp::IntegerMatrix::ConstRow& xobs); // maybe switch to template
  float getloglik_marginal(const Hyperparameter& hparams);
  
public:
  void countReset();
  void countAdd(const Rcpp::IntegerMatrix::ConstRow& xobs);
  void countSubtract(const Rcpp::IntegerMatrix::ConstRow& xobs);
  void reduce_items(Rcpp::IntegerVector items_new, const Hyperparameter& hparams);
  void drop_item(int item, const Hyperparameter& hparams);
  
public:
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
void DomainCount::set_initial(Rcpp::IntegerVector& items_in, const Hyperparameter& hparams, const Rcpp::NumericVector& lthetas_in) {
  TROUBLE_START(("DomainCount::set_initial #V1"));
  items = items_in;
  
  set_pattern2id_map(hparams.item_nlevels);
  
  if (lthetas_in.size() > 0) {
    lthetas = lthetas_in;
  } else {
    lthetas = Rcpp::NumericVector(npatterns);
  }
  
  countReset();
  
  TROUBLE_END;
}

//' @name DomainCount::set_pattern2id_map
//' @title DomainCount::set_pattern2id_map
//' @description Create 'conversion vector' for converting (vector) response pattern to representative ID
//' See id2pattern(.) for more details
//' Side effect. Sets npatterns
//' @keywords internal
void DomainCount::set_pattern2id_map(const Rcpp::IntegerVector& item_nlevels) {
  TROUBLE_START(("DomainCount::set_pattern2id_map"));
  
  pattern2id_map = Rcpp::IntegerVector(item_nlevels.size(), 0);
  int ndomainitems = ndomainitems_calc();
  
  if (ndomainitems==0) {
    npatterns = 0;
    TROUBLE_END; return;
  }
  
  int cumprod_current = 1;
  int iitem;
  for (int i = 0; i < ndomainitems; i++) {
    iitem = items[i];
    pattern2id_map[iitem] = cumprod_current;
    cumprod_current *= item_nlevels[iitem];
  }
  
  npatterns = cumprod_current; // Should match npatterns = lthetas.size() and product(item_nlevels[theseItems])
  TROUBLE_END;
}

//' @name DomainCount::set_initial
//' @title DomainCount::set_initial
//' @description As other set_initial() but with list compatibilities for R
//' @param list_domain List of all arguments for set_initial()
//' @param hparams hyperparmeters
//' @keywords internal
void DomainCount::set_initial(Rcpp::List list_domain, const Hyperparameter& hparams) {
  TROUBLE_START(("DomainCount::set_initial #V2"));
  Rcpp::NumericVector lthetas_in =  Rcpp::log(list_domain["thetas"]);
  Rcpp::IntegerVector items_in = list_domain["items"];
  set_initial(items_in, hparams, lthetas_in);
  TROUBLE_END;
}

//' @name DomainCount::reduce_items
//' @title DomainCount::reduce_items
//' @description Efficiently remove items from domain (without needing to re-calculate counts)
//' @param items_new Resulting items in your domain. Assumed(!!) to be subset of current items.
//' @keywords internal
void DomainCount::reduce_items(Rcpp::IntegerVector items_new, const Hyperparameter& hparams) {
  TROUBLE_START(("DomainCount::reduce_items"));
  
  DomainCount domain_orig = copy();
  set_initial(items_new, hparams);
  for (int i=0; i < domain_orig.npatterns; i++) { // implicit, ignores ndomains ndomainitems_calc()==0
    counts[pattern2id(domain_orig.id2pattern(i))] += domain_orig.counts[i];
  }
  TROUBLE_END;
}

//' @name DomainCount::drop_item
//' @title DomainCount::drop_item
//' @description Efficiently remove item from domain (without needing to re-calculate counts)
//' @param item item you wish to remove
//' @keywords internal
void DomainCount::drop_item(int item, const Hyperparameter& hparams) {
  TROUBLE_START(("DomainCount::drop_item"));
  
  Rcpp::IntegerVector items_new = items[items!=item];
  reduce_items(items_new, hparams);
  
  TROUBLE_END;
}

//' @name DomainCount::pattern2id
//' @title DomainCount::pattern2id
//' @description Convert (vector) response pattern to representative ID
//' See id2pattern(.) for more details
//' @keywords internal
template <typename vectype>
int DomainCount::pattern2id(const vectype& xobs) {
  TROUBLE_START(("DomainCount::pattern2id"));
  TROUBLE_END; return Rcpp::sum(xobs * pattern2id_map);
}

//' @name DomainCount::get_ltheta
//' @title DomainCount::get_ltheta
//' @description Look up the log-probability of seeing this response pattern in this domain
//' @param xobs vector of the FULL response pattern (not just the items in this domain)
//' @keywords internal
double DomainCount::get_ltheta(const Rcpp::IntegerMatrix::ConstRow& xobs) {
  TROUBLE_START(("DomainCount::get_ltheta"));
  double out;
  if (npatterns > 0) {
    out = lthetas(pattern2id(xobs));
  } else {
    out = 0; // log 1
  }
  TROUBLE_END; return out;
}

//' @name DomainCount::id2pattern
//' @title DomainCount::id2pattern
//' @description Convert pattern ID to original response pattern
//' See utilities ::id2pattern(.) for more details
//' Differs from ::id2pattern(.) in that pattern2id_map is a cumulative compared to mapvec
//' Assumes items are in same order as pattern2id_map;
//' @keywords internal
Rcpp::IntegerVector DomainCount::id2pattern(int id) {
  TROUBLE_START(("DomainCount::id2pattern #V1"));
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
  
  TROUBLE_END; return pattern;
}

//' @name DomainCount::countReset
//' @title DomainCount::countReset
//' @description Set counts to zero. 
//' Counts measure the number of times each pattern appears in the data.
//' @keywords internal
void DomainCount::countReset() {
  TROUBLE_START(("DomainCount::id2pattern #V2"));
  counts = Rcpp::IntegerVector(npatterns, 0);
  TROUBLE_END;
}

//' @name DomainCount::countAdd
//' @title DomainCount::countAdd
//' @description Add one observation to this domain's counts.
//' Counts measure the number of times each pattern appears in the data.
//' @param xobs One FULL response pattern (not just the items in this domain)
//' @keywords internal
void DomainCount::countAdd(const Rcpp::IntegerMatrix::ConstRow& xobs) {
  TROUBLE_START(("DomainCount::countAdd"));
  counts[pattern2id(xobs)] += 1;
  TROUBLE_END;
}

//' @name DomainCount::countSubtract
//' @title DomainCount::countSubtract
//' @description Removes one observation to this domain's counts.
//' @param xobs One FULL response pattern (not just the items in this domain)
//' @keywords internal
void DomainCount::countSubtract(const Rcpp::IntegerMatrix::ConstRow& xobs) {
  TROUBLE_START(("DomainCount::countSubtract"));
  counts[pattern2id(xobs)] -= 1;
  TROUBLE_END;
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
std::vector<std::map<int,  DomainCount> > DomainCount::list2domains(Rcpp::List list_list_domains, const Hyperparameter& hparams) {
  TROUBLE_START(("DomainCount::list2domains"));
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
  
  TROUBLE_END; return domains;
}

//' @name DomainCount::copy
//' @title DomainCount::copy
//' @description Creates a deep copy of this domain
//' @keywords internal
DomainCount DomainCount::copy() {
  TROUBLE_START(("DomainCount::copy"));
  DomainCount newDomain;
  newDomain.lthetas = Rcpp::clone(lthetas);
  newDomain.items = Rcpp::clone(items);
  newDomain.pattern2id_map = Rcpp::clone(pattern2id_map);
  newDomain.npatterns = npatterns;
  newDomain.counts = Rcpp::clone(counts);
  TROUBLE_END; return newDomain;
}

//' @name DomainCount::itemsid_calc
//' @title DomainCount::itemsid_calc
//' @description Convert domain.items into an integer id describing the items in this domain
//' @keywords internal
int DomainCount::itemsid_calc() {
  TROUBLE_START(("DomainCount::itemsid_calc"));
  int n = ndomainitems_calc();
  int sum = 0;
  for (int i=0; i < n; i++) {
    sum += powl(2,items[i]); // not optimized for speed. improve by storing vector of 2^k in hparams. Maybe save this to avoid recalculating too
  }
  TROUBLE_END; return sum;
}

//' @name DomainCount::getloglik_marginal
//' @title DomainCount::getloglik_marginal
//' @description For a given domain, we want to know the probability of observing a series of responses.
//' We calculate this probability conditioning on parameters, but marginalizing over theta (and convert to log scale)
//' For domain_theta_prior_type="patternadjusted", we also include partial prior for domain
//' @param hparams hyperparameters
//' @keywords internal
float DomainCount::getloglik_marginal(const Hyperparameter& hparams) {
  TROUBLE_START(("DomainCount::getloglik_marginal"));
  if (npatterns == 0) {
    TROUBLE_END; return 0; // log(1)
  }
  
  Rcpp::NumericVector theta_alpha = hparams.theta_alpha + Rcpp::NumericVector(npatterns);
  float ldenominator=0;
  if (hparams.domain_theta_prior_type != "patternadjusted") {
    ldenominator=lbeta(theta_alpha);
  } // if hparams.domain_theta_prior_type=="patternadjusted", applies partial prior for domain by setting ldenominator=0
  float loglik = (
    lbeta(Rcpp::as<Rcpp::NumericVector> (counts) + theta_alpha)
    - ldenominator
  );
  TROUBLE_END; return loglik;
}



//' @name DomainCount::print
//' @title DomainCount::print
//' @description Prints this domain (used mainly for troubleshooting purposes)
//' @keywords internal
void DomainCount::print() {
  TROUBLE_START(("DomainCount::print"));
  Rcpp::Rcout << "domain.lthetas:" << lthetas << "\n";
  Rcpp::Rcout << "domain.items:" << items << "\n";
  Rcpp::Rcout << "domain.pattern2id_map:" << pattern2id_map << "\n";
  Rcpp::Rcout << "domain.npatterns:" << npatterns << "\n";
  Rcpp::Rcout << "domain.counts:" << counts << "\n";
  TROUBLE_END;
}

// GLOBAL CONSTANTS - DO NOT MODIFY PROGRAMATICALLY!
DomainCount BLANK_DOMAIN; // Example blank domain
float fINF = std::numeric_limits<float>::infinity();
double dINF = std::numeric_limits<double>::infinity();


/*****************************************************
 ****** Helpers for Bayes Parameters
 *****************************************************/


// Output for BayesParameter::domain_proposal(.)
struct domainProposalOut {
  int class2domain_id;
  int swap_type;
  int domain0_id;
  int domain1_id;
  Rcpp::IntegerVector items0_old;
  Rcpp::IntegerVector items1_old;
  Rcpp::IntegerVector items0_new;
  Rcpp::IntegerVector items1_new;
  Rcpp::IntegerMatrix item2domainid_new;
  float forwardProb;
  float backwardProb;
  Rcpp::IntegerVector domain_classes;
  std::vector<std::map<int,  DomainCount> > domains_new;
  int is_ok;
};

// Output for BayesParameter::domain_accept(.)
struct domainAcceptOut {
  float loglik_old;
  float loglik_new;
  float unif;
  float log_cutoff;
  bool accept;
};


/*****************************************************
 ****** Bayes Parameters
 *****************************************************/


// Bayes parameters  (random variables)
class BayesParameter {
public:
  Rcpp::NumericVector class_pi; // Prior for classes
  Rcpp::IntegerVector classes; // Class for each observation
  std::vector<std::map<int,  DomainCount> > domains; // Domains and thetas
  
  public: // Inferred namely from domains
    Rcpp::IntegerMatrix item2domainid; // Function of domains
    Rcpp::IntegerMatrix domains_accept; // Function of domains
    Rcpp::NumericMatrix class_loglik; // Log likelihood of seeing this pattern conditional on class and all parameters
    Rcpp::NumericMatrix class_loglik_collapsed; // Log likelihood of seeing this pattern on domains and previous class membership
    
public:
  void set_initial(Rcpp::NumericVector class_pi_in, Rcpp::IntegerVector classes_in, std::vector<std::map<int,  DomainCount> > domains_in, const Hyperparameter& hparams);
  void set_initial(Rcpp::List list_bparam, const Hyperparameter& hparams);
  float class_lprob(const Rcpp::IntegerMatrix::ConstRow& xobs, int xclass);
  Rcpp::NumericVector class_lprob(const Rcpp::IntegerMatrix::ConstRow& xobs);
  void set_class_loglik(const Rcpp::IntegerMatrix& x, bool reset=false);
  void set_class_loglik_collapsed(const Rcpp::IntegerMatrix& x, const Hyperparameter &hparams);
  
public:
  void domain_resetCounts();
  void domain_addCount(const Rcpp::IntegerMatrix::ConstRow& xobs, int xclass);
  void domain_addCounts(const Rcpp::IntegerMatrix& x, bool reset_counts=true);
  void domain_resetCounts(std::vector<std::map<int,  DomainCount> >& domains);
  void domain_addCount(const Rcpp::IntegerMatrix::ConstRow& xobs, int xclass, std::vector<std::map<int,  DomainCount> >& domains);
  void domain_addCounts(const Rcpp::IntegerMatrix& x, const Rcpp::IntegerVector& classes, bool reset_counts, std::vector<std::map<int,  DomainCount> >& domains);
  
public:
  int nclass_calc() {return class_pi.size();};
  int nitems_calc() {return domains[0].begin()->second.nitems_calc();};
  int nobs_calc() {return classes.size();};
  Rcpp::IntegerVector nobs_classes_calc();
  Rcpp::IntegerMatrix item2domainid_calc(const Hyperparameter& hparams);
  int domain_id_new(int class2domain_id, const Hyperparameter& hparams);
  static Rcpp::IntegerVector get_superdomains(const Rcpp::IntegerMatrix& item2domainid);
  static bool is_identifiable(Rcpp::IntegerVector& item2superdomainid, const Hyperparameter& hparams);
  bool is_identifiable(domainProposalOut& proposal, const Hyperparameter& hparams);
  
public:
  Rcpp::NumericVector class_pi_args(const Hyperparameter& hparams);
  void class_pi_next(const Hyperparameter& hparams);
  void classes_next(const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams, const bool is_collapsed);
  void thetas_next(const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams);
  
public:
  static float domain_prior(Rcpp::IntegerVector& item2domainid_vec, const Hyperparameter& hparams);
  domainProposalOut domain_proposal(int class2domain_id, const Hyperparameter& hparams, int maxitr=5);
  domainAcceptOut domain_accept(const Rcpp::IntegerMatrix& x, domainProposalOut& proposal, const Hyperparameter& hparams);
  int domain_next(int class2domain_id, const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams);
  void domains_next(const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams);
};

//' @name BayesParameter::set_initial
//' @title BayesParameter::set_initial
//' @description Set all BayesParameter properties
//' See getStart_bayes_params(.) in dependentLCM.r for more details
//' @keywords internal
void BayesParameter::set_initial(Rcpp::NumericVector class_pi_in, Rcpp::IntegerVector classes_in, std::vector<std::map<int,  DomainCount> > domains_in, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::set_initial #V1"));
  class_pi = class_pi_in;
  classes = classes_in;
  domains = domains_in;
  item2domainid = item2domainid_calc(hparams);
  domains_accept = Rcpp::IntegerMatrix(hparams.domain_nproposals, hparams.nclass2domain);
  domains_accept.fill(-1);
  if (hparams.steps_active["likelihood"] == true) {
    class_loglik = Rcpp::NumericMatrix(hparams.nclass, hparams.nobs); // initialize with all zeroes
  }
  if (hparams.steps_active["class_collapse"] == true) {
    class_loglik_collapsed = Rcpp::NumericMatrix(hparams.nclass, hparams.nobs); // initialize with all zeroes
  }
  TROUBLE_END;
}

//' @name BayesParameter::set_initial
//' @title BayesParameter::set_initial
//' @description Set all BayesParameter properties
//' As other set_initial(.) but supports lists for R compatability
//' @keywords internal
void BayesParameter::set_initial(Rcpp::List list_bparam, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::set_initial #V2"));
  Rcpp::NumericVector class_pi_in = list_bparam("class_pi");
  Rcpp::IntegerVector classes_in = list_bparam("classes");
  Rcpp::List list_domains_in = list_bparam("domains");
  std::vector<std::map<int,  DomainCount> > domains_in = DomainCount::list2domains(list_domains_in, hparams);
  set_initial(class_pi_in, classes_in, domains_in, hparams);
  TROUBLE_END;
}

//' @name BayesParameter::nobs_classes_calc
//' @title BayesParameter::nobs_classes_calc
//' @description How many observations are in each class?
Rcpp::IntegerVector BayesParameter::nobs_classes_calc() {
  int nclass = nclass_calc();
  Rcpp::IntegerVector class_nobs (nclass);
  for (int iclass=0; iclass<nclass; iclass++) {
    class_nobs[iclass] = Rcpp::sum(domains[iclass].begin()->second.counts);
  }
  return class_nobs;
};

//' @name BayesParameter::class_lprob
//' @title BayesParameter::class_lprob
//' @description What is the (log) probability that this response pattern was generated by this class?
//' @param xobs The response pattern we are investigation.
//' @param xclass The class this observation is (assumed to be) in
//' @keywords internal
float BayesParameter::class_lprob(const Rcpp::IntegerMatrix::ConstRow& xobs, int xclass) {
  TROUBLE_START(("BayesParameter::class_lprob #V1"));
  float lprob = 0;
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end = domains[xclass].end();
  
  
  for (domain_iter = domains[xclass].begin(); domain_iter!=domain_end; ++domain_iter) {
    lprob += domain_iter->second.get_ltheta(xobs);
  }
  
  TROUBLE_END; return lprob;
}

//' @name BayesParameter::class_lprob
//' @title BayesParameter::class_lprob
//' @description What is the (log) probability that this response pattern was generated by EACH class?
//' Assumes all parameters are fixed/known.
//' @param xobs The response pattern we are investigation.
//' @keywords internal
Rcpp::NumericVector BayesParameter::class_lprob(const Rcpp::IntegerMatrix::ConstRow& xobs) {
  TROUBLE_START(("BayesParameter::class_lprob #V2"));
  Rcpp::NumericVector lprobs = Rcpp::NumericVector(nclass_calc());
  
  for (int i=0; i < nclass_calc(); i++) {
    lprobs[i] = class_lprob(xobs, i);
  }
  TROUBLE_END; return lprobs;
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
  TROUBLE_START(("BayesParameter::set_class_loglik"));
  // One column per observation and row per class. Get conditional likelihood
  
  int xnrow = x.nrow();
  
  if (reset == true) {
    // Initialize matrix
    class_loglik = Rcpp::NumericMatrix(nclass_calc(), xnrow);
  }
  
  for (int i = 0; i < xnrow; i++) {
    class_loglik.column(i) = class_lprob(x.row(i));
  }
  TROUBLE_END;
}

//' @name BayesParameter::set_class_loglik_collapsed
//' @title BayesParameter::set_class_loglik_collapsed
//' @description For each observation what is the collapsed (log) probability of this observation belonging to this class? Collapsed on theta and pi.
//' Updates corresponding property in BayesParameter instance. Only correct up to proportionality.
//' Assumes all parameters are fixed/known.
//' @param x Matrix of responses
//' @param hparams Hyperparameter values
//' @keywords internal
void BayesParameter::set_class_loglik_collapsed(const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::set_class_loglik_collapsed"));
  // One column per observation and row per class
  // influence via pi
  Rcpp::NumericVector nobs_classes = Rcpp::as<Rcpp::NumericVector>(nobs_classes_calc());
  Rcpp::NumericMatrix log_pi_collapsed = Rcpp::NumericMatrix(hparams.nclass, hparams.nclass);
  for (int iclass=0; iclass<hparams.nclass; iclass++) {
    
    // get dirichlet parameters
    log_pi_collapsed.column(iclass) = nobs_classes + hparams.classPi_alpha;
    log_pi_collapsed(iclass, iclass) -= 1; // do not include iobs'th observation itself
    
    // get expected value
    log_pi_collapsed.column(iclass) = (
      Rcpp::log(log_pi_collapsed.column(iclass))
      - std::log(Rcpp::sum(log_pi_collapsed.column(iclass)))
    );
  }
  for (int iobs = 0; iobs < hparams.nobs; iobs++) {
    class_loglik_collapsed.column(iobs) = log_pi_collapsed.column(classes[iobs]);
  }
    
  
  // influence via theta
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  Rcpp::NumericVector icounts;
  Rcpp::NumericVector ilog_theta_collapsed;
  Rcpp::NumericVector ilog_theta_collapsed1;
  for (int iclass = 0; iclass < hparams.nclass; iclass++) {
    
    domain_end = domains[iclass].end();
    for (domain_iter = domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      
      // get dirichlet parameters
      icounts = (
        Rcpp::as<Rcpp::NumericVector>(domain_iter->second.counts)
        + hparams.theta_alpha
       );
      
      // get expected values
      ilog_theta_collapsed = (
        Rcpp::log(icounts) 
        - std::log(Rcpp::sum(icounts))
      );
      ilog_theta_collapsed1 = (
        Rcpp::log(icounts-1)
        - std::log(Rcpp::sum(icounts) - 1)
      );  // subtract one to not count iobs'th observation itself
      for (int iobs = 0; iobs < hparams.nobs; iobs++) {
        if (iclass == classes[iobs]) {
          class_loglik_collapsed(iclass, iobs) += ilog_theta_collapsed1[domain_iter->second.pattern2id(x.row(iobs))];
        } else {
          class_loglik_collapsed(iclass, iobs) += ilog_theta_collapsed[domain_iter->second.pattern2id(x.row(iobs))];
        }
      }
    }
  }
  
  TROUBLE_END;
}

//' @name BayesParameter::domain_resetCounts
//' @title BayesParameter::domain_resetCounts
//' @description Reset the counts of all domains in given map
//' @keywords internal
void BayesParameter::domain_resetCounts(std::vector<std::map<int,  DomainCount> >& domains) {
  TROUBLE_START(("BayesParameter::domain_resetCounts #V1"));
  int nclass = nclass_calc();
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  
  for (int iclass=0; iclass < nclass; iclass++) {
    domain_end = domains[iclass].end();
    for (domain_iter = domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      domain_iter->second.countReset();
    }
  }
  TROUBLE_END;
}

//' @name BayesParameter::domain_resetCounts
//' @title BayesParameter::domain_resetCounts
//' @description Reset the counts in all of BayesParameter's domains
//' @keywords internal
void BayesParameter::domain_resetCounts() {
  TROUBLE_START(("BayesParameter::domain_resetCounts #V2"));
  domain_resetCounts(domains);
  TROUBLE_END;
}

//' @name BayesParameter::domain_addCount
//' @title BayesParameter::domain_addCount
//' @description Add given observation to all the counts in all provided domains
//' @param xobs the response pattern we are counting
//' @param xlcass The class of this observation. Needed because different domains correspond to different classes.
//' @param domains map of domains we wish to add this pattern to
//' @keywords internal
void BayesParameter::domain_addCount(const Rcpp::IntegerMatrix::ConstRow& xobs, int xclass, std::vector<std::map<int,  DomainCount> >& domains) {
  TROUBLE_START(("BayesParameter::domain_addCount #V1"));
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  domain_end = domains[xclass].end();
  for (domain_iter = domains[xclass].begin(); domain_iter!=domain_end; ++domain_iter) {
    domain_iter->second.countAdd(xobs);
  }
  TROUBLE_END;
}

//' @name BayesParameter::domain_addCount
//' @title BayesParameter::domain_addCount
//' @description Add given observation to all the counts in all of this BayesParameter's domains
//' @param xobs the response pattern we are counting
//' @param xlcass The class of this observation. Needed because different domains correspond to different classes.
//' @keywords internal
void BayesParameter::domain_addCount(const Rcpp::IntegerMatrix::ConstRow& xobs, int xclass) {
  TROUBLE_START(("BayesParameter::domain_addCount #V2"));
  domain_addCount(xobs, xclass, domains);
  TROUBLE_END;
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
  TROUBLE_START(("BayesParameter::domain_addCounts #V1"));
  
  if (reset_counts == true) {
    domain_resetCounts(domains);
  }
  int nobs = x.nrow();
  for (int obs_id=0; obs_id < nobs; obs_id++) {
    domain_addCount(x.row(obs_id), classes[obs_id], domains);
  }
  TROUBLE_END;
}

//' @name BayesParameter::domain_addCounts
//' @title BayesParameter::domain_addCounts
//' @description Add given observations to counts in all domains in BayesParameter
//' @param x Matrix of the response patterns
//' @param reset_counts True if we should set all counts to zero before counting x.
//' @keywords internal
void BayesParameter::domain_addCounts(const Rcpp::IntegerMatrix& x, bool reset_counts) {
  TROUBLE_START(("BayesParameter::domain_addCounts #V2"));
  domain_addCounts(x, classes, reset_counts, domains);
  TROUBLE_END;
}

//' @name BayesParameter::item2domainid_calc
//' @title BayesParameter::item2domainid_calc
//' @description For each item find the id for the domain it belongs to.
//' If there are multiple class2domain, then multiple columns are provided (one for each)
//' @keywords internal
Rcpp::IntegerMatrix BayesParameter::item2domainid_calc(const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::item2domainid_calc"));
  
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
  
  TROUBLE_END; return out;
}

//' @name BayesParameter::domain_prior
//' @title BayesParameter::domain_prior
//' @description Before observing data/other-parameters, how likely are we to put items into these particular domains?
//' Note that for domain_theta_prior_type="patternadjusted" some of the prior is also given in DomainCount::getloglik_marginal()
//' Some choices of domains may be more likely than other based on prior.
//' @keywords internal
float BayesParameter::domain_prior(Rcpp::IntegerVector& item2domainid_vec, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::domain_prior"));
  
  if (hparams.domain_theta_prior_type=="niave") {
    TROUBLE_END; return 0; // log(1)
  }
  
  int ndomains_nonempty = count_unique(item2domainid_vec);
  
  float loglik = (
    std::lgamma(hparams.ndomains + 1) // permuting groups
    - std::lgamma(hparams.ndomains - ndomains_nonempty + 1) // permuting groups
    - hparams.nitem * std::log(hparams.ndomains) // denominator
  );
  
  std::list<std::pair<int,int>>::const_iterator iter;
  std::list<std::pair<int,int>>::const_iterator iter_end = hparams.domainPriorKnowledgeLog_pairs.end();
  for (iter=hparams.domainPriorKnowledgeLog_pairs.begin(); iter!=iter_end; iter++) {
    if (item2domainid_vec[iter->first] == item2domainid_vec[iter->second]) {
      loglik += hparams.domainPriorKnowledgeLog_mat(iter->first, iter->second);
    }
  }
  
  TROUBLE_END; return loglik; // Assume flat prior
}

//' @name BayesParameter::get_superdomains
//' @title BayesParameter::get_superdomains
//' @description Merge overlapping domains from different class2domainid
//' @param item2domainid Each colum describes what items must be grouped together for this item2domainid
//' @keywords internal
Rcpp::IntegerVector BayesParameter::get_superdomains(const Rcpp::IntegerMatrix& item2domainid) {
  TROUBLE_START(("BayesParameter::get_superdomains"));
  int nitems = item2domainid.nrow();
  int nclass2domain = item2domainid.ncol();
  int i;
  
  if (nclass2domain == 1) {
    // Nothing to merge. Return item to domain associations
    TROUBLE_END; return Rcpp::IntegerVector(item2domainid.column(0));
  }
  
  // Find connections across class2domain
  Rcpp::IntegerMatrix adjmat = Rcpp::IntegerMatrix(nitems, nitems);
  for (i=0; i < nclass2domain; i++) {
    adjmat += equal_to_adjmat(item2domainid.column(i));
    // Can we make faster by ignoring singleton domains?
  }
  
  // Merge linked nodes
  Rcpp::IntegerVector item2superdomainid = adjmat_to_equal(adjmat, nitems);
  
  TROUBLE_END; return item2superdomainid;
}

//' @name BayesParameter::is_identifiable
//' @title BayesParameter::is_identifiable
//' @description Check if choice of domains is generically identifiable
//' Uses greedy algorithm. May fail in some cases, but is deterministic (if bad then consistently conservative for that choice of domains)
//' See Allman paper (DOI:10.1214/09-AOS689 Theorem 4.) for criteria used: min(patterns1,nclass)+min(patterns2,nclass)+min(patterns3,nclass) > 2*nclass+2
//' @param item2superdomainid Vector describing which items must be grouped together
//' @param hparams hyperparameters
//' @keywords internal
bool BayesParameter::is_identifiable(Rcpp::IntegerVector& item2superdomainid, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::is_identifiable #V1"));
  
  bool identifiable = ::is_identifiable(item2superdomainid, hparams.nclass, hparams.item_nlevels)
    
    TROUBLE_END; return identifiable;
}

//' @name BayesParameter::is_identifiable
//' @title BayesParameter::is_identifiable
//' @description Check if choice of domains is generically identifiable under proposal
//' @param proposal The candidate changes to domains we are evaluating
//' @param hparams hyperparameters
//' @keywords internal
bool BayesParameter::is_identifiable(domainProposalOut& proposal, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::is_identifiable #V2"));
  
  // Check identifiability of proposal
  Rcpp::IntegerVector item2superdomainid_proposed = get_superdomains(proposal.item2domainid_new);
  bool identifiable = is_identifiable(item2superdomainid_proposed, hparams);
  
  TROUBLE_END; return identifiable ;
}

//' @name BayesParameter::domain_id_new
//' @title BayesParameter::domain_id_new
//' @description Find an empty domain. Return that domain's ID
//' @keywords internal
int BayesParameter::domain_id_new(int class2domain_id, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::domain_id_new"));
  
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
  
  TROUBLE_END; return domain_id;
}


/*****************************************************
 ****** MCMC
 *****************************************************/

//' @name BayesParameter::class_pi_args
//' @title BayesParameter::class_pi_args
//' @description Calculate the dirichlet parameters for the posterior of pi (pi used for classes)
//' @keywords internal
Rcpp::NumericVector BayesParameter::class_pi_args(const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::class_pi_args"));
  
  Rcpp::NumericVector args = Rcpp::clone(hparams.classPi_alpha);  // Prevent overwrites. Maybe instead remove & from input
  
  Rcpp::IntegerVector::iterator classes_itr;
  Rcpp::IntegerVector::const_iterator classes_end = classes.end();
  for (classes_itr = classes.begin();
       classes_itr != classes_end;
       ++classes_itr) {
    args(*classes_itr) += 1;
  }
  
  TROUBLE_END; return args;
}

//' @name BayesParameter::class_pi_next
//' @title BayesParameter::class_pi_next
//' @description Do gibbs sampling to generate pi (pi used for classes)
//' @keywords internal
void BayesParameter::class_pi_next(const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::class_pi_next"));
  
  Rcpp::NumericVector args = class_pi_args(hparams);
  class_pi = rDirichlet(args);
  TROUBLE_END;
}

//' @name BayesParameter::classes_next
//' @title BayesParameter::classes_next
//' @description Do gibbs sampling to calculate the class of each observation
//' @keywords internal
void BayesParameter::classes_next(const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams, const bool is_collapsed) {
  TROUBLE_START(("BayesParameter::classes_next"));
  
  Rcpp::NumericMatrix *iclass_loglik;
  Rcpp::NumericVector class_args_init(nclass_calc());
  Rcpp::NumericMatrix class_args_init_base;
  if (is_collapsed==false) {
    iclass_loglik = &class_loglik;
    class_args_init = Rcpp::log(class_pi);
  } else {
    iclass_loglik = &class_loglik_collapsed;
    // leave class_args_init at zero
  }
  
  int nrow = x.nrow();
  Rcpp::NumericVector class_args(nclass_calc());
  for (int i=0; i < nrow; i++) {
    class_args = (
      class_args_init
      + (*iclass_loglik).column(i));
    class_args = class_args - Rcpp::max(class_args); // Divide by max for precision: [e^a,e^b,e^c] = e^a[1, e^(b-a), e^(c-a)]
    class_args = Rcpp::exp(class_args);
    class_args = class_args / Rcpp::sum(class_args);
    classes(i) = rCategorical(class_args);
  }
  TROUBLE_END;
}

//' @name BayesParameter::thetas_next
//' @title BayesParameter::thetas_next
//' @description Do gibbs sampling to calculate thetas for each domain
//' @keywords internal
void BayesParameter::thetas_next(const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::thetas_next"));
  
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  DomainCount *idomain;
  Rcpp::NumericVector iconcentration;
  for (int iclass=0; iclass < hparams.nclass; iclass++) {
    domain_end = domains[iclass].end();
    for (domain_iter = domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      iconcentration = idomain->counts + hparams.theta_alpha;
      idomain->lthetas = Rcpp::log(rDirichlet(iconcentration));
    }
  }
  TROUBLE_END;
}

//' @name BayesParameter::domain_proposal
//' @title BayesParameter::domain_proposal
//' @description Propose new domains in metropolis algorithm
//' Essentially we choose two domains at random and move items between them (oversimplification).
//' We take a nonempty domain. 1) With some probability we move one of its items to an empty domain. Otherwise we choose a second nonempty domain. 2) With some probability we swap a pair of items between the two nonempty domains. 3) Otherwise we move one item from the first domain into the second.
//' There are some restrictions including 1) We do not swap items if both domains are singletons (this would cause no change up to relabeling), and enforce maximum number of items per domain given in hparams.
//' @keywords internal
domainProposalOut BayesParameter::domain_proposal(int class2domain_id, const Hyperparameter& hparams, int maxitr) {
  TROUBLE_START(("BayesParameter::domain_proposal"));
  
  std::map<int, int> domainitem_counts = count_integers(item2domainid.column(class2domain_id));
  
  /*
   * Choose Domains
   */
  
  Rcpp::IntegerVector domain_classes = which(hparams.class2domain == class2domain_id); // maybe make fixed
  Rcpp::IntegerVector domains_nonempty = Rcpp::unique(item2domainid.column(class2domain_id));
  Rcpp::IntegerVector domains_chosen;
  if (domains_nonempty.size() >= 2) {
    // always do this except in edge case
    domains_chosen = Rcpp::sample(domains_nonempty, 2, false);
  } else {
    // only single domain (edge case when identifiability is off)
    domains_chosen = {Rcpp::sample(domains_nonempty, 1, false)[0], -1};
  }
  
  Rcpp::IntegerVector items0_old = which(item2domainid.column(class2domain_id) == domains_chosen[0]);
  int items0_old_size = items0_old.size();
  
  /*
   * Decide type
   */
  
  int swap_type;
  if (domains_chosen[1] == -1) {
    swap_type = 2; // split
  } else if (items0_old_size <= 1) {
    // If size 1 then nothing to split
    swap_type = 1; // mix
  } else if (R::runif(0,1) < hparams.domain_proposal_empty) {
    swap_type = 2; // split
  } else {
    swap_type = 1; // mix
  }
  
  if (swap_type==2) { // split
    domains_chosen[1] = -1;
  }
  
  if (domains_chosen[1] == -1) {
    // pick new domain
    domains_chosen[1] = domain_id_new(class2domain_id, hparams);
  }
  
  /*
   * Get old item info
   */
  
  Rcpp::IntegerVector items1_old = which(item2domainid.column(class2domain_id) == domains_chosen[1]);
  
  int items1_old_size = items1_old.size();
  int items_old_size = items0_old_size + items1_old_size;
  
  Rcpp::IntegerVector items_old (items_old_size);
  Rcpp::IntegerVector mask_old (items_old_size);
  for (int i=0; i < items_old_size; i++) {
    if (i < items0_old_size) {
      items_old[i] = items0_old[i];
      mask_old[i] = 0;
    } else {
      // items0_old_size <= i < items0_old_size+items1_old_size
      items_old[i] = items1_old[i-items0_old_size];
      mask_old[i] = 1;
    }
  }
  
  
  /*
   * Move Items
   */
  
  Rcpp::IntegerVector items0_new;
  Rcpp::IntegerVector items1_new;
  
  // get new mask
  int is_ok;
  Rcpp::IntegerVector mask_new (items_old_size);
  if (items_old_size <= 2) {
    // special case. Separated strictly for speed
    mask_new = Rcpp::clone(mask_old);
    int i = (int)(R::runif(0, 1) < 0.5);
    mask_new[i] = 1-mask_new[i]; // flip
    is_ok = 0; // Unknown
  } else {
    // items_old_size > 2
    int diff;
    is_ok = -1; // bad
    for (int i=0; i<maxitr; i++) {
      for (int item_pos=0; item_pos<items_old_size; item_pos++) {
        mask_new[item_pos] = (R::runif(0, 1) < 0.5);
      }
      
      
      diff = Rcpp::sum(mask_new-mask_old);
      if (
          (diff == 0)
        | (diff == items_old_size)
      ) {
        ; // no change. Loop and try again
      } else if (
          (Rcpp::sum(mask_new==0) > hparams.domain_maxitems)
        | (Rcpp::sum(mask_new==1) > hparams.domain_maxitems)
      ) {
        ; // Too many items. Loop and try again
      } else {
        is_ok = 1; // good
        break; // Success, quit loop
      }
    }
  }
  
  items0_new = items_old[mask_new==0];
  items1_new = items_old[mask_new==1];
  std::sort(items0_new.begin(),items0_new.end());
  std::sort(items1_new.begin(),items1_new.end());
  
  // item2domainid
  Rcpp::IntegerMatrix item2domainid_new = Rcpp::clone(item2domainid);
  int domain_id;
  Rcpp::IntegerVector items;
  for (int i=0; i<2; i++) {
    domain_id = domains_chosen[i];
    if (i == 0) {
      items = items0_new;
    } else if (i == 1) {
      items = items1_new;
    }
    for (int j=0; j < items.size(); j++) {
      item2domainid_new(items[j], class2domain_id) = domain_id;
    }
  }
  
  /*
   * Proposal Probabilities
   */
  
  float forwardProb;
  float backwardProb;
  int domains_nonempty_size = domains_nonempty.size();
  bool empty_created = (
    (items0_new.size() == 0)
    | (items1_new.size() == 0)
  );
  
  if (swap_type == 1) {
    // mix
    if (empty_created == false) {
      // equally likely in either direction
      forwardProb = (
        ((items0_old.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
      + ((items1_old.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
      );
      backwardProb = (
        ((items0_new.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
        + ((items1_new.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
      );
    } else {
      forwardProb = (
        1 / float(domains_nonempty_size) // Choosing domain1
      * 1 / float(domains_nonempty_size-1) // Choosing domain2
      * (
          ((items0_old.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
      + ((items1_old.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
      ) // No split (need to consider both directions)
      * 1 // every mix equally likely up to proportionality
      );
      backwardProb = (
        1 / float(domains_nonempty_size-1) // Choosing domain
        * hparams.domain_proposal_empty // split
        * 1 // every mix equally likely up to proportionality
      );
    }
  } else if (swap_type == 2) {
    // split
    forwardProb = (
      1 / float(domains_nonempty_size) // Choosing domain1
    * hparams.domain_proposal_empty // split
    * 1 // every split equally likely up to proportionality
    );
    backwardProb = (
      1 / float(domains_nonempty_size+1) // Choosing domain1
      * 1 / float(domains_nonempty_size) // Choosing domain2
      * (
          ((items0_new.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
      + ((items1_new.size() > 1) ? (1 - hparams.domain_proposal_empty) : 1)
      ) // No split (need to consider both directions)
      * 1 // every mix equally likely up to proportionality
    );
  }
  
  
  /*
   * As domain
   */
  
  DomainCount domain0_new;
  DomainCount domain1_new;
  domain0_new.set_initial(items0_new, hparams);
  domain1_new.set_initial(items1_new, hparams);
  
  std::vector<std::map<int,  DomainCount> > domains_new;
  domains_new.resize(hparams.nclass);
  int iclass;
  for (int i = 0; i < domain_classes.size(); i++) {
    iclass = domain_classes[i];
    if (domain0_new.ndomainitems_calc() > 0) {
      domains_new[iclass][domains_chosen[0]] = domain0_new.copy();
    }
    if (domain1_new.ndomainitems_calc() > 0) {
      domains_new[iclass][domains_chosen[1]] = domain1_new.copy();
    }
  }
  
  /*
   * Write output
   */
  
  domainProposalOut  proposal;
  proposal.class2domain_id = class2domain_id;
  proposal.swap_type = swap_type;
  proposal.domain0_id = domains_chosen[0];
  proposal.domain1_id = domains_chosen[1];
  proposal.items0_old = items0_old;
  proposal.items1_old = items1_old;
  proposal.items0_new = items0_new;
  proposal.items1_new = items1_new;
  proposal.item2domainid_new = item2domainid_new;
  proposal.forwardProb = forwardProb;
  proposal.backwardProb = backwardProb;
  proposal.domain_classes = domain_classes;
  proposal.domains_new = domains_new;
  proposal.is_ok = is_ok;
  
  TROUBLE_END; return proposal;
}

//' @name BayesParameter::domain_accept
//' @title BayesParameter::domain_accept
//' @description In metroplis algorithm when choosing to update domains.
//' We examine the proposal and decide whether to accept it.
//' ASSUMPTIONS domain_addCounts has been run on proposal.domains_new
//' @keywords internal
domainAcceptOut BayesParameter::domain_accept(const Rcpp::IntegerMatrix& x, domainProposalOut& proposal, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::domain_accept"));
  
  domainAcceptOut out;
  
  Rcpp::IntegerVector item2domainid_vec;
  item2domainid_vec = item2domainid.column(proposal.class2domain_id);
  item2domainid_vec = Rcpp::clone(item2domainid_vec); // do we need to explcitly clone?
  out.loglik_old = domain_prior(item2domainid_vec, hparams);
  item2domainid_vec = proposal.item2domainid_new.column(proposal.class2domain_id);
  item2domainid_vec = Rcpp::clone(item2domainid_vec); // do we need to explcitly clone?
  out.loglik_new = domain_prior(item2domainid_vec, hparams);
  
  
  int iclass;
  for (int i = 0; i < proposal.domain_classes.size(); i++) {
    
    iclass = proposal.domain_classes[i];
    
    if (domains[iclass].count(proposal.domain0_id) > 0) {
      out.loglik_old += domains[iclass][proposal.domain0_id].getloglik_marginal(hparams);
    }
    if (domains[iclass].count(proposal.domain1_id) > 0) {
      out.loglik_old += domains[iclass][proposal.domain1_id].getloglik_marginal(hparams);
    }
    if (proposal.domains_new[iclass].count(proposal.domain0_id) > 0) {
      out.loglik_new += proposal.domains_new[iclass][proposal.domain0_id].getloglik_marginal(hparams);
    }
    if (proposal.domains_new[iclass].count(proposal.domain1_id) > 0) {
      out.loglik_new += proposal.domains_new[iclass][proposal.domain1_id].getloglik_marginal(hparams);
    }
  }
  
  out.unif = R::runif(0, 1);
  out.log_cutoff = out.loglik_new - out.loglik_old + std::log(proposal.backwardProb) - std::log(proposal.forwardProb);
  out.accept = int( out.log_cutoff > std::log(out.unif) );
  
  
  TROUBLE_END; return out;
}

//' @name BayesParameter::domain_next
//' @title BayesParameter::domain_next
//' @description We use a metropolis algorithm to try to update domains
//' In essence we choose one item at random and either move it to another domain or swap it with an item from another domain
//' @keywords internal
int BayesParameter::domain_next(int class2domain_id, const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::domain_next"));
  
  /***
   *** Propose
   ***/
  
  domainProposalOut proposal =  domain_proposal(class2domain_id, hparams);
  
  /***
   *** Terminate early if...
   ***/
  
  int accept = 0;
  if (proposal.domain0_id == proposal.domain1_id) {
    // No change
    accept = 2;
    TROUBLE_END; return accept;
  }
  if (proposal.is_ok == -1) {
    // No change up to reordering of domains. Do not bother changing
    // (also checks too many items in a given domain, partially. future might give these separate ids)
    accept = -1;
  }
  if ((proposal.items1_new.size() > hparams.domain_maxitems)
        | (proposal.items1_new.size() > hparams.domain_maxitems)) {
    // Over the max items per domain, reject
    accept = -2;
  } else if (hparams.steps_active["identifiable"]==true) {
    if (is_identifiable(proposal, hparams)==false) {
      // Identifiability restrictions violated
      accept = -2;
    }
  }
  if (accept != 0) {
    TROUBLE_END; return accept;
  }
  
  /***
   *** Get Likelihood
   ***/
  
  domain_addCounts(x, classes, true, proposal.domains_new); // true not strictly necessary
  domainAcceptOut accept_info = domain_accept(x,proposal, hparams);
  accept = accept_info.accept;
  
  /***
   *** Apply changes
   ***/
  
  if (accept == 0) {
    TROUBLE_END; return accept; // Change nothing
  }
  
  item2domainid = proposal.item2domainid_new; 
  int iclass;
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  for (int i = 0; i < proposal.domain_classes.size(); i++) {
    iclass = proposal.domain_classes[i];
    domain_end = proposal.domains_new[iclass].end();
    for (domain_iter = proposal.domains_new[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      domains[iclass][domain_iter->first] = domain_iter->second;
    }
    if (proposal.items0_new.size() == 0) {
      domains[iclass].erase(proposal.domain0_id);
    }
    if (proposal.items1_new.size() == 0) {
      domains[iclass].erase(proposal.domain1_id);
    }
  }
  TROUBLE_END; return accept; // Is true
}

//' @name BayesParameter::domains_next
//' @title BayesParameter::domains_next
//' @description Use metropolis algorithm to update domains. Repeat # of times set by hparams.
//' @keywords internal
void BayesParameter::domains_next(const Rcpp::IntegerMatrix& x, const Hyperparameter& hparams) {
  TROUBLE_START(("BayesParameter::domains_next"));
  
  for (int i=0; i < hparams.domain_nproposals; i++) {
    for (int iclass2domain=0; iclass2domain < hparams.nclass2domain; iclass2domain++) {
      domains_accept(i, iclass2domain) = domain_next(iclass2domain, x, hparams);
    }
  }
  TROUBLE_END;
}


/*****************************************************
 ****** Archive
 *****************************************************/


// Used to record the MCMC results
class Archive {
  
public:
  Rcpp::NumericMatrix class_pi;
  Rcpp::IntegerMatrix classes;
  Rcpp::IntegerMatrix class_counts;
  std::vector<Rcpp::IntegerMatrix> domains_id; // Row for iter, classid,  domain_id, pattern_id, items_id
  std::vector<Rcpp::NumericVector> domains_lprobs;
  std::vector<Rcpp::IntegerMatrix> domains_accept;
  std::vector<Rcpp::NumericMatrix> class_loglik;
  std::vector<Rcpp::NumericMatrix> class_loglik_collapsed;
  int next_itr;
public:
  Rcpp::NumericVector obsLik;
  Rcpp::NumericVector obsLogLik;
  Rcpp::NumericVector obsLogLik2;
  Rcpp::NumericVector itrLogLik;
  int nitrLik;
  
public:
  void set_initial(const Hyperparameter& hparams);
  void add(BayesParameter& aparams, const Hyperparameter& hparams);
  void domains2mat(BayesParameter& params, int itr, Rcpp::IntegerMatrix& out_domains_id, Rcpp::NumericVector& out_domains_lprobs);
  Rcpp::List toList();
};

//' @name Archive::set_initial
//' @title Archive::set_initial
//' @description Initializes the archive (namely reserving memory)
//' @keywords internal
void Archive::set_initial(const Hyperparameter& hparams) {
  TROUBLE_START(("Archive::set_initial"));
  next_itr = 0;
  
  int nitr_all = hparams.save_itrs["all"];
  
  class_pi = Rcpp::NumericMatrix(hparams.nclass, nitr_all);
  classes = Rcpp::IntegerMatrix(hparams.nobs, hparams.save_itrs["classes"]);
  
  class_counts = Rcpp::IntegerMatrix(hparams.nclass, hparams.nobs);
  
  domains_id.resize(nitr_all);
  domains_lprobs.resize(nitr_all);
  
  domains_accept.resize(hparams.save_itrs["domains_accept"]);
  std::fill(domains_accept.begin(), domains_accept.end()
              , Rcpp::IntegerMatrix(hparams.domain_nproposals, hparams.nclass2domain)
  ); // does initializing actually speed things up later? I considered removing the std::vector and then implicitly working with a 3d matrix by way of a super-long 1d rcpp:vector (and attr("dim")). But this was unstable
  
  class_loglik.resize(hparams.save_itrs["class_loglik"]);
  std::fill(class_loglik.begin(), class_loglik.end()
              , Rcpp::NumericMatrix(hparams.nclass, hparams.nobs));
  
  class_loglik_collapsed.resize(hparams.save_itrs["class_loglik_collapsed"]);
  std::fill(class_loglik_collapsed.begin(), class_loglik_collapsed.end()
              , Rcpp::NumericMatrix(hparams.nclass, hparams.nobs));
  
  
  int save_itrs_agg_loglik = hparams.save_itrs["agg_loglik"]; // direct call causing issues downstream
  if (save_itrs_agg_loglik>0) {
    obsLik = Rcpp::NumericVector(hparams.nobs);
    obsLogLik = Rcpp::NumericVector(hparams.nobs);
    itrLogLik = Rcpp::NumericVector(save_itrs_agg_loglik);
    obsLogLik2 = Rcpp::NumericVector(hparams.nobs);
  }
  nitrLik = 0;
  
  TROUBLE_END;
}

//' @name Archive::domains2mat
//' @title Archive::domains2mat
//' @description Helper function for Archive::add
//' Used to convert the domains from map form to matrix form for storage
//' Although we could make a deep copy of the domain map each time, this would be unproductive because we need it in matrix form later for R. Therfore we convert to matrix.
//' @keywords internal
void Archive::domains2mat(BayesParameter& params, int itr, Rcpp::IntegerMatrix& out_domains_id, Rcpp::NumericVector& out_domains_lprobs) {
  TROUBLE_START(("Archive::domains2mat"));
  
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  int iclass;
  DomainCount *idomain;
  int i_npatterns;
  int idomain_id;
  int iitems_id;
  int ithis_pattern_id;
  int i;
  // misc
  int npatterns;
  int nclass = params.nclass_calc();
  // int nitems = params.nitems_calc();
  
  
  // Set Size
  npatterns = 0;
  for (iclass=0; iclass < nclass; iclass++) {
    domain_end = params.domains[iclass].end();
    for (domain_iter = params.domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      npatterns += idomain->npatterns;
    }
  }
  out_domains_id = Rcpp::IntegerMatrix(5, npatterns); // Row for iter, classid,  domain_id, pattern_id, items_id
  out_domains_lprobs = Rcpp::NumericVector(npatterns);
  
  
  ithis_pattern_id = 0;
  for (iclass=0; iclass < nclass; iclass++) {
    domain_end = params.domains[iclass].end();
    for (domain_iter = params.domains[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      i_npatterns = idomain->npatterns;
      idomain_id = domain_iter->first;
      iitems_id = domain_iter->second.itemsid_calc();
      for (i = 0; i < i_npatterns; i++) {
        // Save info
        out_domains_id(0, ithis_pattern_id) = itr;
        out_domains_id(1, ithis_pattern_id) = iclass;
        out_domains_id(2, ithis_pattern_id) = idomain_id;
        out_domains_id(3, ithis_pattern_id) = i;
        out_domains_id(4, ithis_pattern_id) = iitems_id;
        out_domains_lprobs(ithis_pattern_id) = idomain->lthetas(i);
        
        ithis_pattern_id += 1;  //Increment
      }
    }
  }
  
  TROUBLE_END;
}

//' @name Archive::add
//' @title Archive::add
//' @description Add latest bayes parameters to the archive
//' @keywords internal
void Archive::add(BayesParameter& aparams, const Hyperparameter& hparams) {
  TROUBLE_START(("Archive::add"));
  
  // Iterations are added by column for speed
  // Currently assumes fixed number of iterations. For flexible iterations std::list would allow for easier extension, but slower.
  
  int inext_itr_all = next_itr - ( hparams.nitr - hparams.save_itrs["all"] );
  
  if (inext_itr_all < 0) {
    // too early, add nothing
    next_itr += 1;
    TROUBLE_END; return;
  }
  
  class_pi.column(inext_itr_all) = aparams.class_pi;
  Archive::domains2mat(aparams, next_itr
                         , domains_id[inext_itr_all], domains_lprobs[inext_itr_all]);
  
  int inext_itr;
  
  inext_itr = next_itr - ( hparams.nitr - hparams.save_itrs["classes"] );
  if (inext_itr >= 0) {
    classes.column(inext_itr) = aparams.classes;
  }
  
  inext_itr = next_itr - ( hparams.nitr - hparams.save_itrs["class_counts"] );
  if (inext_itr >= 0) {
    for (int iobs=0; iobs<hparams.nobs; iobs++) {
      class_counts(aparams.classes[iobs], iobs) += 1L;
    }
  }
  
  inext_itr = next_itr - ( hparams.nitr - hparams.save_itrs["domains_accept"] );
  if (inext_itr >= 0) {
    domains_accept[inext_itr] = Rcpp::clone(aparams.domains_accept);
  }
  
  inext_itr = next_itr - ( hparams.nitr - hparams.save_itrs["class_loglik"] );
  if (inext_itr >= 0) {
    class_loglik[inext_itr] = Rcpp::clone(aparams.class_loglik);
  }
  
  inext_itr = next_itr - ( hparams.nitr - hparams.save_itrs["class_loglik_collapsed"] );
  if (inext_itr >= 0) {
    class_loglik_collapsed[inext_itr] = Rcpp::clone(aparams.class_loglik_collapsed);
  }
  
  inext_itr = next_itr - ( hparams.nitr - hparams.save_itrs["agg_loglik"] );
  if (inext_itr >= 0) {
    Rcpp::NumericVector iobsLogLik (hparams.nobs);
    for (int iobs=0; iobs<hparams.nobs; iobs++) {
      iobsLogLik[iobs] = sumLogs(aparams.class_loglik.column(iobs) + log(aparams.class_pi));
    }
    
    obsLogLik += iobsLogLik;
    obsLogLik2 += iobsLogLik*iobsLogLik;
    obsLik += exp(iobsLogLik);
    itrLogLik[inext_itr] = sum(iobsLogLik);
    nitrLik += 1;
  }
  
  next_itr += 1;
  TROUBLE_END;
}

//' @name Archive::toList
//' @title Archive::add
//' @description Add latest bayes parameters to the archive
//' @keywords internal
Rcpp::List Archive::toList() {
  TROUBLE_START(("Archive::toList"));
  TROUBLE_END; return Rcpp::List::create(
      Rcpp::Named("class_pi") = class_pi
    , Rcpp::Named("classes") = classes
    , Rcpp::Named("class_counts") = class_counts
    , Rcpp::Named("domains_id") = wrap(domains_id)
    , Rcpp::Named("domains_lprobs") = wrap(domains_lprobs)
    , Rcpp::Named("domains_accept") = domains_accept
    , Rcpp::Named("class_loglik") = class_loglik
    , Rcpp::Named("class_loglik_collapsed") = class_loglik_collapsed
    , Rcpp::Named("obsLik") = obsLik
    , Rcpp::Named("obsLogLik") = obsLogLik
    , Rcpp::Named("obsLogLik2") = obsLogLik2
    , Rcpp::Named("itrLogLik") = itrLogLik
    , Rcpp::Named("nitrLik") = nitrLik
    , Rcpp::Named("troubleshooting") = TROUBLE_LIST
  );
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
public:
  void set_initial(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list);
  void next(const bool skip_class_update);
  Archive run();
  Archive run_init(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list);
};

//' @name BayesContainer::set_initial
//' @title BayesContainer::set_initial
//' @description Sets all properties of BayesContainer
//' @keywords internal
void BayesContainer::set_initial(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list) {
  TROUBLE_START(("BayesContainer::set_initial"));
  x = x_in;
  hparams.set_hparams(hparams_list);
  hparams.set_dataInfo(x_in);
  params.set_initial(params_list, hparams);
  TROUBLE_END;
}

//' @name BayesContainer::next
//' @title BayesContainer::next
//' @description Does a single MCMC step on all bayes parameters
//' @keywords internal
void BayesContainer::next(const bool skip_class_update) {
  TROUBLE_START(("BayesContainer::run"));
  
  if ((hparams.steps_active["class_collapse"]==true) & !skip_class_update) {
    params.set_class_loglik_collapsed(x, hparams);
  }
  
  if ((hparams.steps_active["classes"]==true) & !skip_class_update) {
    params.classes_next(x, hparams, hparams.steps_active["class_collapse"]);
  }
  
  if (
    (hparams.steps_active["thetas"]==true)
    | (hparams.steps_active["domains"]==true)
  ) {
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
  
  if (hparams.steps_active["likelihood"]) {
    params.set_class_loglik(x); // always last because it recalculates likelihood for above choices
  }
  
  TROUBLE_END;
}

//' @name BayesContainer::run
//' @title BayesContainer::run
//' @description Does hparams.nitr MCMC steps on all bayes parameters
//' @keywords internal
Archive BayesContainer::run() {
  TROUBLE_START(("BayesContainer::run"));
  
  Archive archive;
  archive.set_initial(hparams);
  
  archive.add(params, hparams); // save initial value (iteration 0)
  
  if (0 < hparams.nitr) {
    // In first iteration we update parameters but not classes (iteration 1)
    next(true);
    archive.add(params, hparams);
  }
  
  for (int i=2; i < hparams.nitr; i++) {
    next(false);
    archive.add(params, hparams);
    print_itr(i, hparams.print_itrs);
  }
  TROUBLE_END; return archive;
}

//' @name BayesContainer::run_init
//' @title BayesContainer::run_init
//' @description First initializes, and then does hparams.nitr MCMC steps on all bayes parameters
//' @keywords internal
Archive BayesContainer::run_init(const Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list) {
  TROUBLE_START(("BayesContainer::run_init"));
  set_initial(x_in, hparams_list, params_list);
  Archive archive = run();
  TROUBLE_END; return archive;
}


/*****************************************************
 ****** Public Functions
 *****************************************************/


//' @name dependentLCM_fit_cpp
//' @title dependentLCM_fit_cpp
//' @description Does MCMC simulations for domain LCM model.
//' This is a C++ script. Run dependentLCM_fit in dependentLCM.r to execute.
//' @param x_in Matrix of responses we are analyzing
//' @param hparams_list List of hyperparameter info. See getStart_hparams() in R.
//' @param params_list List of parameter info. See getStart_bayes_params() in R.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List dependentLCM_fit_cpp(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list) {
  TROUBLE_INIT;
  TROUBLE_START(("dependentLCM_fit_cpp"));
  BayesContainer bcontainer;
  Archive archive = bcontainer.run_init(x_in, hparams_list, params_list);
  TROUBLE_END; return archive.toList();
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
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix id2pattern(const Rcpp::IntegerVector& xpattern, const Rcpp::IntegerVector& mapvec) {
  int npatterns = xpattern.size();
  int nmapvec = mapvec.size();
  Rcpp::IntegerMatrix unmapped_mat = Rcpp::IntegerMatrix(nmapvec, npatterns);
  
  for (int i = 0; i < npatterns; i++) {
    unmapped_mat.column(i) = id2pattern(xpattern(i), mapvec);
  }
  
  return unmapped_mat; 
}

//' @name itemid2patterns
//' @title itemid2patterns
//' @description Take pattern_id/item_id combinations and form pattern.
//' @param pattern_ids ID reprsenting the value of 'filled in' items. (unfilled items get -1)
//' @param items_ids ID reprsenting which items (positions) should be 'filled'. 
//' @param item_nlevels How many possible values each item can take
//' @examples
//' \dontrun{
//' # Find the pattern of each row in domains from IDs
//' domain_patterns <- itemid2patterns(
//'   pattern_ids = dlcm_out$mcmc$domains$pattern_id
//'   , items_ids = dlcm_out$mcmc$domains$items_id
//'   , item_nlevels = dlcm_out$hparams$item_nlevels
//'   )
//' }
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector itemid2patterns(const Rcpp::IntegerVector& pattern_ids, const Rcpp::IntegerVector& items_ids, const Rcpp::IntegerVector& item_nlevels) {
  int n = pattern_ids.size();
  int nitems = item_nlevels.size();
  Rcpp::IntegerVector items_2s (nitems, 2);
  Rcpp::IntegerMatrix patterns (nitems, n);
  std::map<std::tuple<int,int>, Rcpp::IntegerVector> lookups;
  
  std::tuple<int,int> key;
  Rcpp::LogicalVector items;
  Rcpp::IntegerVector pattern (nitems);
  for (int i=0; i<n; i++) {
    key = std::make_pair(items_ids[i], pattern_ids[i]);
    if (lookups.count(key) > 0) {
      // use stored pattern
      patterns.column(i) = Rcpp::clone(lookups[key]);
    } else {
      // find pattern
      items = id2pattern(std::get<0>(key), items_2s)==1;
      pattern.fill(-1);
      pattern[items] = id2pattern(std::get<1>(key), item_nlevels[items]);
      lookups[key] = Rcpp::clone(pattern);
      patterns.column(i) = pattern;
    }
  }
  
  return patterns;
}

//' @name get_pattern2id_map
//' @title get_pattern2id_map
//' @description Create 'conversion vector' for converting (vector) response pattern to representative ID
//' The dot product of the pattern2id_map and the response pattern gives the pattern id. This does the opposite of id2pattern(.).
//' @param items_ids ID reprsenting which items (positions) should be 'filled'. 
//' @param item_nlevels How many possible values each item can take.
//' @param items_addend Should be -1 when run in R and 0 when run in C++. This is added to items_ids.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_pattern2id_map(const Rcpp::IntegerVector& item_nlevels, const Rcpp::IntegerVector& items, const int items_addend=-1) {
  TROUBLE_START(("get_pattern2id_map"));
  // Wraps DomainCount::set_pattern2id_map(.) for use in R.
  
  DomainCount idomain;
  idomain.items = items + items_addend;
  idomain.set_pattern2id_map(item_nlevels);
  
  TROUBLE_END; return idomain.pattern2id_map;
}



//' @name expSumLog
//' @title expSumLog
//' Same as log(sum(exp(x))), but adjusted to improve precision.
//' Assumes x are logged values. Calculates sum(e^x) and then converts back to log scale
//' Handle precision: log(e^a+e^b+e^c) = log(a * (1+e^(b-a)+e^(c-1))) = log(a) + log(1+e^(b-a)+e^(c-1)))
//' @param x numeric vector in log scale
//' @keywords internal
//' @export
// [[Rcpp::export]]
double expSumLog(const Rcpp::NumericVector& x) {
  return sumLogs(x);
}

//' @name is_identifiable_r
//' @title is_identifiable_r
//' @description Check if choice of domains is generically identifiable
//' Uses greedy algorithm. May give false negatives in some cases, but is quick and deterministic
//' @param item2domainid Integer matrix with one row per item. For homogeneous DLCM give a single column. In this column items with the same value are interpreted as being in the same domain. For heterogeneous DLCM, provide one column per class. In a single column, items with the same value are assumed to be in the same domain for that class. Putting every item into its own domain is analogous to checking generic identifiability of a traditional LCM.
//' @param nclass Number of classes
//' @param item_nlevels Vector with levels of each item
//' @export
// [[Rcpp::export]]
bool is_identifiable_r(const Rcpp::IntegerMatrix& item2domainid, int nclass, const Rcpp::IntegerVector& item_nlevels) {
  
  // group into conditionally indpendent blocks (i.e. for heterogeneous)
  Rcpp::IntegerVector item2superdomainid = BayesParameter::get_superdomains(item2domainid);
  
  bool identifiable = is_identifiable(item2superdomainid, nclass, item_nlevels);
  
  return identifiable;
}