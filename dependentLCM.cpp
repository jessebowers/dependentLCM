#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]


/* 
 * Open items:
 * Enforce identifiability
 * Allow skipping/randomization of steps (including archive)
 * Allow warmup without archiving
 * incorporate checks/consistency checks eg domain size, nclasses, etc
 * Make sure there is no unnecessary copying
 * Look at vector.push and note that this always creates a whole new vector which is slow. Though resizing upfront means only one copy
 * post-hoc handle label switching
 * DomainCount initialize counts in set_initial instead of having to set/fill later
 * Skip counts if symmetric change to DomainCounts (only for common alpha)
 * Apply retrictions:
 * 1) Max items per domain
 * 2) Classes with same domains
 * 3) Fix items to be in the same domain (post-hoc)
 * Make ndomains=0 as infinite domains
 * Do gibbs instead??
 * Try shuffling domains??
 * Save pattern id to archive and enforce ordering of items
 * reduce domains to just nonempty
 */


/*****************************************************
 ****** GLOBAL CONSTANTS
 *****************************************************/


bool IS_LOUD = false; // tkprint for troubleshooting


/*****************************************************
 ****** UTILITIES
 *****************************************************/


Rcpp::IntegerVector colMax(Rcpp::IntegerMatrix& x) {
  if (IS_LOUD) {Rcpp::Rcout << "colMax" << "\n";} // tkprint
  // Purpose: Get maximum value of each column of x
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

Rcpp::NumericVector rDirchlet(const Rcpp::NumericVector& deltas) {
  if (IS_LOUD) {Rcpp::Rcout << "rDirchlet" << "\n";} // tkprint
  // Purpose: Randomly generate one observatation of Dirchlet distribution
  int C = deltas.size();
  Rcpp::NumericVector Xgamma(C);
  
  // generating gamma(deltac,1)
  for (int c = 0; c < C; c++) {
    Xgamma(c) = R::rgamma(deltas(c), 1.0);  // tk maybe swap out for c function
    //Xgamma(c) = Rcpp::rgamma(1, deltas(c), scale = 1.0);
  }
  return Xgamma / sum(Xgamma);
}

int rCategorical(const Rcpp::NumericVector& probs) {
  if (IS_LOUD) {Rcpp::Rcout << "rCategorical" << "\n";} // tkprint
  // Purpose: Randomly generate one observatation of a categorical distibution
  
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


int count_unique(Rcpp::IntegerVector& x) {
  if (IS_LOUD) {Rcpp::Rcout << "count_unique" << "\n";} // tkprint
  std::unordered_set<int> xset(x.begin(), x.end());
  return xset.size();
}

Rcpp::IntegerVector rowVec(Rcpp::IntegerMatrix& mat, int i) {
  if (IS_LOUD) {Rcpp::Rcout << "rowVec" << "\n";} // tkprint
  Rcpp::IntegerVector irow = mat.row(i);  // Convert
  return irow;
}

float lbeta(const Rcpp::NumericVector& alpha) {
  if (IS_LOUD) {Rcpp::Rcout << "lbeta" << "\n";} // tkprint
  float log_gamma_total = std::lgamma(Rcpp::sum(alpha));
  float log_gammas = Rcpp::sum(Rcpp::lgamma(alpha));
  
  return (log_gammas - log_gamma_total);
}

float beta(const Rcpp::NumericVector& alpha) {
  if (IS_LOUD) {Rcpp::Rcout << "beta" << "\n";} // tkprint
  return std::exp(lbeta(alpha));
}

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x) {
  if (IS_LOUD) {Rcpp::Rcout << "which" << "\n";} // tkprint
  int n = x.size();
  std::list<int> out; // linked list for fast append
  
  for(int i = 0; i < n; i++) {
    if (x[i]) { // If x is true
      out.push_back(i);
    }
  }
  
  return Rcpp::wrap(out);
}


std::map<int,  int> count_integers(const Rcpp::IntegerVector& x) {
  
  std::map<int,  int> counts_map;
  int i;
  int nx = x.size();
  int ix;
  for (i=0; i < nx; i++) {
    ix = x[i];
    counts_map.insert({ix, 0});
    counts_map[ix] += 1;
  }
  
  // Rcpp::IntegerMatrix counts_mat = Rcpp::IntegerMatrix(Rcpp::no_init_matrix(counts_map.size(), 2));
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

template <typename K, typename V>
V map_get(const  std::map <K,V> & map, const K & key, const V & defaultvalue ) {
  typename std::map<K,V>::const_iterator iter = map.find( key );
  if ( iter == map.end() ) {
    return defaultvalue;
  }
  else {
    return iter->second;
  }
}

template <typename T>
T minimum(const T x1, const T x2) {
  if (x1 < x2) {
    return x1;
  } else {
    return x2;
  }
}


/*****************************************************
 ****** Hyperparameters
 *****************************************************/


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
  float domain_alpha_one_calc() {if (ndomains>0){return domain_alpha/float(ndomains);} else {return 0;};}; // treat 0 as infinity
  void set_hparams(int ndomains_in, int nclass_in, Rcpp::IntegerVector& class2domain_in, Rcpp::NumericVector& classPi_alpha_in, float domain_alpha_in, int domain_maxitems_in, float theta_alpha_in, float domain_proposal_empty_in, float domain_proposal_swap_in, int domain_nproposals_in);
  void set_hparams(Rcpp::List hparams_in);
  void set_dataInfo(Rcpp::IntegerMatrix& x);
  void print();
};

void Hyperparameter::set_hparams(
    int ndomains_in
  , int nclass_in
  , Rcpp::IntegerVector& class2domain_in
  , Rcpp::NumericVector& classPi_alpha_in
  , float domain_alpha_in
  , int domain_maxitems_in
  , float theta_alpha_in
  , float domain_proposal_empty_in
  , float domain_proposal_swap_in
  , int domain_nproposals_in) {
  if (IS_LOUD) {Rcpp::Rcout << "Hyperparameter::set_hparams" << "\n";} // tkprint
  
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
  
  // Inferred
  nclass2domain = nclass2domain_calc();
}

void Hyperparameter::set_hparams(Rcpp::List hparams_in) {
  if (IS_LOUD) {Rcpp::Rcout << "Hyperparameter::set_hparams" << "\n";} // tkprint
  // Same as other set_hparams but with list compatability for R
  
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
  
  set_hparams(ndomains, nclass, class2domain, classPi_alpha, domain_alpha, domain_maxitems, theta_alpha, domain_proposal_empty, domain_proposal_swap, domain_nproposals);
}


void Hyperparameter::set_dataInfo(Rcpp::IntegerMatrix& x) {
  if (IS_LOUD) {Rcpp::Rcout << "Hyperparameter::set_dataInfo" << "\n";} // tkprint
  // Purpose: Get information about x data
  // Assumptions: That there are no empty levels especially at end
  
  item_nlevels = colMax(x) + 1; // length(0:n) = n+1
  nobs = x.nrow();
  nitem = nitem_calc();
}

void Hyperparameter::print() {
  if (IS_LOUD) {Rcpp::Rcout << "Hyperparameter::print" << "\n";} // tkprint
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
}


/*****************************************************
 ****** Domains
 *****************************************************/


class Domain {
public:
  Rcpp::NumericVector thetas;
  Rcpp::IntegerVector items;
  Rcpp::IntegerVector pattern2id_map;
  int npatterns = 0; // tk make better
  
public:
  void set_initial(Rcpp::IntegerVector& items_in, Hyperparameter& hparams, const Rcpp::NumericVector& thetas_in = Rcpp::NumericVector(0));
  void set_initial(Rcpp::List list_domain, Hyperparameter& hparams);
  void set_pattern2id_map(Hyperparameter& hparams);
  int pattern2id(Rcpp::IntegerMatrix::Row xobs); //tk maybe switch to template
  Rcpp::IntegerVector id2pattern(int id);
  float get_theta(Rcpp::IntegerMatrix::Row xobs); //tk maybe switch to template
  
public:
  int ndomainitems_calc() {return items.size();};
  int nitems_calc() {return pattern2id_map.size();};
  void print();
  
  // Rcpp::IntegerVector item_nlevels(Hyperparameter hparams){return hparams.item_nlevels(items);}; tk errors out for some reason
};

void Domain::set_initial(Rcpp::IntegerVector& items_in, Hyperparameter& hparams, const Rcpp::NumericVector& thetas_in) {
  if (IS_LOUD) {Rcpp::Rcout << "Domain::set_initial" << "\n";} // tkprint
  items = items_in;
  set_pattern2id_map(hparams);
  
  if (thetas_in.size() > 0) {
    thetas = thetas_in;
  } else {
    thetas = static_cast<Rcpp::NumericVector>(Rcpp::no_init_vector(npatterns));
  }
  
  thetas = thetas_in;
}

void Domain::set_pattern2id_map(Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "Domain::set_pattern2id_map" << "\n";} // tkprint
  // Side effect. Sets npatterns
  pattern2id_map = static_cast<Rcpp::NumericVector>(Rcpp::no_init_vector(hparams.nitem));
  pattern2id_map.fill(0);
  
  int cumprod_current = 1;
  int iitem;
  int ndomainitems = ndomainitems_calc();
  for (int i = 0; i < ndomainitems; i++) {
    iitem = items[i];
    pattern2id_map[iitem] = cumprod_current;
    cumprod_current *= hparams.item_nlevels[iitem];
  }
  
  npatterns = cumprod_current; // Should match npatterns = thetas.size() and product(item_nlevels[theseItems])
}

void Domain::set_initial(Rcpp::List list_domain, Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "Domain::set_initial" << "\n";} // tkprint
  // Purpose: For r compatability convert list to domain
  Rcpp::NumericVector thetas_in =  list_domain["thetas"];
  Rcpp::IntegerVector items_in = list_domain["items"];
  set_initial(items_in, hparams, thetas_in);
}

int Domain::pattern2id(Rcpp::IntegerMatrix::Row xobs) {
  if (IS_LOUD) {Rcpp::Rcout << "Domain::pattern2id" << "\n";} // tkprint
  // Purpose: Which thetas index corresponds with these item values?
  return Rcpp::sum(xobs * pattern2id_map); 
}

float Domain::get_theta(Rcpp::IntegerMatrix::Row xobs) {
  if (IS_LOUD) {Rcpp::Rcout << "Domain::get_theta" << "| nitems" << items.size() << "| npatterns" << npatterns << "\n";} // tkprint
  // Purpose: Extract one theta value
  
  return thetas(pattern2id(xobs));
}

Rcpp::IntegerVector Domain::id2pattern(int id) {
  if (IS_LOUD) {Rcpp::Rcout << "Domain::id2pattern" << "\n";} // tkprint
  // Assumes items are in same order as pattern2id_map;
  Rcpp::IntegerVector pattern(nitems_calc(), -1);
  
  int i_item;
  int i_divisor;
  int i_value;
  int ndomainItems = ndomainitems_calc();
  for (int i=ndomainItems-1; i > -1; i--) {
    i_item = items[i];
    i_divisor = pattern2id_map[i_item];
    if (i_divisor == 0) {
      Rcpp::warning("Domain::id2pattern:: Zero divisor.");
      break; // should never happen. Error out gracefully.
    }
    i_value = id / i_divisor;  // Silently truncates by design
    pattern[i_item] = i_value;
    id = id - i_value * i_divisor;
  }
  
  return pattern;
}

void Domain::print() {
  if (IS_LOUD) {Rcpp::Rcout << "Domain::print" << "\n";} // tkprint
  Rcpp::Rcout << "domain.thetas:" << thetas << "\n";
  Rcpp::Rcout << "domain.items:" << items << "\n";
  Rcpp::Rcout << "domain.pattern2id_map:" << pattern2id_map << "\n";
  Rcpp::Rcout << "domain.npatterns:" << npatterns << "\n";
}


class DomainCount : public Domain {
  // Extends Domain to allow us to count observations matching each domain pattern
  
public: 
  Rcpp::IntegerVector counts;
  
public:
  void countReset();
  void countAdd(Rcpp::IntegerMatrix::Row xobs);
  static std::vector<std::map<int,  DomainCount> > list2thetas(Rcpp::List list_list_domains, Hyperparameter hparams);
  
public:
  DomainCount copy();
  void print();
};

void DomainCount::countReset() {
  if (IS_LOUD) {Rcpp::Rcout << "DomainCount::countReset" << "\n";} // tkprint
  counts = static_cast<Rcpp::IntegerVector>(Rcpp::no_init_vector(npatterns));
  counts.fill(0);
}

void DomainCount::countAdd(Rcpp::IntegerMatrix::Row xobs) {
  if (IS_LOUD) {Rcpp::Rcout << "DomainCount::countAdd" << "\n";} // tkprint
  counts[pattern2id(xobs)] += 1;
}

std::vector<std::map<int,  DomainCount> > DomainCount::list2thetas(Rcpp::List list_list_domains, Hyperparameter hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "DomainCount::list2thetas" << "\n";} // tkprint
  // Purpose: For r compatability convert list of lists of domains into thetas
  // Ignores list names
  // For now just at DomainCount level but theoretically the same at the domain level too. Just don't need it there and does not inheret well
  
  std::vector<std::map<int,  DomainCount> > thetas;
  thetas.resize(list_list_domains.length());
  
  for (int iclass=0; iclass < list_list_domains.length(); iclass++) {
    Rcpp::List iclass_domains = list_list_domains[iclass]; // slow
    for (int jdomain=0; jdomain < iclass_domains.length(); jdomain++) {
      DomainCount jjdomain;
      jjdomain.set_initial(iclass_domains[jdomain], hparams);
      thetas[iclass].insert({jdomain, jjdomain});
    }
  };
  
  return thetas;
};

DomainCount DomainCount::copy() {
  if (IS_LOUD) {Rcpp::Rcout << "DomainCount::copy" << "\n";} // tkprint
  // tk ideally define at Domain level too then inheret but we would need to be able to convert types which is a pain
  DomainCount newDomain;
  newDomain.thetas = Rcpp::clone(thetas);
  newDomain.items = Rcpp::clone(items);
  newDomain.pattern2id_map = Rcpp::clone(pattern2id_map);
  newDomain.npatterns = npatterns;
  newDomain.counts = Rcpp::clone(counts);
  return newDomain;
}

void DomainCount::print() {
  if (IS_LOUD) {Rcpp::Rcout << "DomainCount::print" << "\n";} // tkprint
  Domain::print();
  Rcpp::Rcout << "domain.counts:" << counts << "\n";
}

DomainCount BLANK_DOMAIN;



/*****************************************************
 ****** Bayes Parameters
 *****************************************************/

struct domainProposalOut {
  int swap_type; // 0 swap, 1 transfer to nonempty domain, 2 split into empty domain
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

struct domainAcceptOut {
  std::vector<std::map<int,  DomainCount> > thetas_new;
  float loglik_old;
  float loglik_new;
  float unif;
  float log_cutoff;
  bool accept;
};

struct domainPropProbOut {
  Rcpp::NumericVector domain_probs;
  std::map<int,  int> domainitem_counts;
}; // List and tuples not working here so using struct

void domainAcceptOut_print(domainAcceptOut x) {
  if (true) {Rcpp::Rcout << "domainAcceptOut_print" << "\n";} // tkprint
  Rcpp::Rcout << "domainAcceptOut.loglik_old" << x.loglik_old << "\n"; // tkprint
  Rcpp::Rcout << "domainAcceptOut.loglik_new" << x.loglik_new << "\n"; // tkprint
  Rcpp::Rcout << "domainAcceptOut.unif" << x.unif << "\n"; // tkprint
  Rcpp::Rcout << "domainAcceptOut.log_cutoff" << x.log_cutoff << "\n"; // tkprint
  Rcpp::Rcout << "domainAcceptOut.accept" << x.accept << "\n"; // tkprint
  if (true) {Rcpp::Rcout << "domainAcceptOut_print" << "| END" << "\n";} // tkprint
} // tkjmb delete

void domainProposalOut_print(domainProposalOut& x) {
  if (true) {Rcpp::Rcout << "domainProposalOut_print" << "\n";} // tkprint
  Rcpp::Rcout << "domainProposalOut.swap_type:" << x.swap_type << "\n";
  Rcpp::Rcout << "domainProposalOut.domain_id1:" << x.domain_id1 << "\n";
  Rcpp::Rcout << "domainProposalOut.domain_id2:" << x.domain_id2 << "\n";
  Rcpp::Rcout << "domainProposalOut.item1_old:" << x.item1_old << "\n";
  Rcpp::Rcout << "domainProposalOut.item2_old:" << x.item2_old << "\n";
  Rcpp::Rcout << "domainProposalOut.forwardProb:" << x.forwardProb << "\n";
  Rcpp::Rcout << "domainProposalOut.backwardProb:" << x.backwardProb << "\n";
  Rcpp::Rcout << "domainProposalOut.domain_classes:" << x.domain_classes << "\n";
  Rcpp::Rcout << "domainProposalOut.domain_old1:" << "\n";
  x.domain_old1->print();
  Rcpp::Rcout << "domainProposalOut.domain_old2:" << "\n";
  x.domain_old2->print();
  Rcpp::Rcout << "domainProposalOut.domain_new1:" << "\n";
  x.domain_new1.print();
  Rcpp::Rcout << "domainProposalOut.domain_new2:" << "\n";
  x.domain_new2.print();
  if (true) {Rcpp::Rcout << "domainProposalOut_print" << "| END" << "\n";} // tkprint
} // tkjmb delete

class BayesParameter {
public:
  Rcpp::NumericVector class_pi;
  Rcpp::IntegerVector classes;
  std::vector<std::map<int,  DomainCount> > thetas;
  Rcpp::IntegerMatrix item2domainid; // Function of thetas
  Rcpp::IntegerMatrix thetas_accept; // Function of thetas
  
public:
  void set_initial(Rcpp::NumericVector class_pi_in, Rcpp::IntegerVector classes_in, std::vector<std::map<int,  DomainCount> > thetas_in, Hyperparameter& hparams);
  void set_initial(Rcpp::List list_bparam, Hyperparameter& hparams);
  float class_prob(Rcpp::IntegerMatrix::Row xobs, int xclass);
  Rcpp::NumericVector class_prob(Rcpp::IntegerMatrix::Row xobs);
  
public:
  void theta_resetCounts();
  void theta_addCount(Rcpp::IntegerMatrix::Row xobs, int xclass);
  void theta_addCounts(Rcpp::IntegerMatrix& x, bool reset_counts=true);
  void theta_resetCounts(std::vector<std::map<int,  DomainCount> >& thetas);
  void theta_addCount(Rcpp::IntegerMatrix::Row xobs, int xclass, std::vector<std::map<int,  DomainCount> >& thetas);
  void theta_addCounts(Rcpp::IntegerMatrix& x, bool reset_counts, std::vector<std::map<int,  DomainCount> >& thetas);
  
public:
  int nclass_calc() {return class_pi.size();};
  int nitems_calc() {return thetas[0].begin()->second.nitems_calc();};
  Rcpp::IntegerMatrix item2domainid_calc(Hyperparameter hparams);
  int domain_id_new(int class2domain_id, Hyperparameter& hparams);
  
public:
  Rcpp::NumericVector class_pi_args(Hyperparameter& hparams);
  void class_pi_next(Hyperparameter& hparams);
  Rcpp::NumericVector class_args(Rcpp::IntegerMatrix::Row xobs);
  void classes_next(Rcpp::IntegerMatrix x);
  void thetas_next(Rcpp::IntegerMatrix& x, Hyperparameter& hparams);
  
public:
  static float domain_getloglik_x(Rcpp::IntegerVector& pattern_counts, float theta_alpha);
  static float domain_getlik_domain(Hyperparameter hparams);
  domainProposalOut domain_proposal(int class2domain_id, Hyperparameter& hparams);
  domainAcceptOut domain_accept(Rcpp::IntegerMatrix& x, domainProposalOut& proposal, Hyperparameter& hparams);
  int domain_next(int class2domain_id, Rcpp::IntegerMatrix& x, Hyperparameter hparams);
  void domains_next(Rcpp::IntegerMatrix& x, Hyperparameter& hparams);
};


void BayesParameter::set_initial(Rcpp::NumericVector class_pi_in, Rcpp::IntegerVector classes_in, std::vector<std::map<int,  DomainCount> > thetas_in, Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::set_initial" << "\n";} // tkprint
  class_pi = class_pi_in;
  classes = classes_in;
  thetas = thetas_in;
  item2domainid = item2domainid_calc(hparams);
  thetas_accept = Rcpp::IntegerMatrix(Rcpp::no_init_matrix(hparams.domain_nproposals, hparams.nclass2domain));
  thetas_accept.fill(-1);
};

void BayesParameter::set_initial(Rcpp::List list_bparam, Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::set_initial" << "\n";} // tkprint
  // Purpose: Same as other set_initial, but with list compatability for R
  Rcpp::NumericVector class_pi_in = list_bparam("class_pi");
  Rcpp::IntegerVector classes_in = list_bparam("classes");
  Rcpp::List list_thetas_in = list_bparam("thetas");
  std::vector<std::map<int,  DomainCount> > thetas_in = DomainCount::list2thetas(list_thetas_in, hparams);
  set_initial(class_pi_in, classes_in, thetas_in, hparams);
};

float BayesParameter::class_prob(Rcpp::IntegerMatrix::Row xobs, int xclass) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::class_prob" << "\n";} // tkprint
  // Purpose: Get probability of xclass generating xobs
  float prob = 1;
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end = thetas[xclass].end();
  
  
  for (domain_iter = thetas[xclass].begin(); domain_iter!=domain_end; ++domain_iter) {
    prob *= domain_iter->second.get_theta(xobs);
  }
  
  return prob;
};

Rcpp::NumericVector BayesParameter::class_prob(Rcpp::IntegerMatrix::Row xobs) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::class_prob" << "\n";} // tkprint
  // Purpose: Get probability of generating xobs for each class
  Rcpp::NumericVector probs = static_cast<Rcpp::NumericVector>(Rcpp::no_init_vector(nclass_calc()));
  
  for (int i=0; i < nclass_calc(); i++) {
    probs[i] = class_prob(xobs, i);
  }
  return probs;
};

void BayesParameter::theta_resetCounts(std::vector<std::map<int,  DomainCount> >& thetas) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::theta_resetCounts" << "\n";} // tkprint
  int nclass = nclass_calc();
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  
  for (int iclass=0; iclass < nclass; iclass++) {
    domain_end = thetas[iclass].end();
    for (domain_iter = thetas[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      domain_iter->second.countReset();
    }
  }
}

void BayesParameter::theta_resetCounts() {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::theta_resetCounts" << "\n";} // tkprint
  theta_resetCounts(thetas);
}

void BayesParameter::theta_addCount(Rcpp::IntegerMatrix::Row xobs, int xclass, std::vector<std::map<int,  DomainCount> >& thetas) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::theta_addCount" << "\n";} // tkprint
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  domain_end = thetas[xclass].end();
  for (domain_iter = thetas[xclass].begin(); domain_iter!=domain_end; ++domain_iter) {
    domain_iter->second.countAdd(xobs);
  }
}

void BayesParameter::theta_addCount(Rcpp::IntegerMatrix::Row xobs, int xclass) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::theta_addCount" << "\n";} // tkprint
  theta_addCount(xobs, xclass, thetas);
}

void BayesParameter::theta_addCounts(Rcpp::IntegerMatrix& x, bool reset_counts, std::vector<std::map<int,  DomainCount> >& thetas) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::theta_addCounts" << "\n";} // tkprint
  
  if (reset_counts == true) {
    theta_resetCounts(thetas);
  }
  
  int nobs = x.nrow();
  for (int obs_id=0; obs_id < nobs; obs_id++) {
    theta_addCount(x.row(obs_id), classes[obs_id], thetas);
  }
}

void BayesParameter::theta_addCounts(Rcpp::IntegerMatrix& x, bool reset_counts) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::theta_addCounts" << "\n";} // tkprint
  theta_addCounts(x, reset_counts, thetas);
}

Rcpp::IntegerMatrix BayesParameter::item2domainid_calc(Hyperparameter hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::item2domainid_calc" << "\n";} // tkprint
  
  Rcpp::IntegerMatrix out = Rcpp::IntegerMatrix(hparams.nitem, hparams.nclass2domain);
  
  int iclass;
  Rcpp::IntegerVector iitems;
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  
  for (int iclass2domain=0; iclass2domain < hparams.nclass2domain; iclass2domain++) {
    iclass = which(hparams.class2domain == iclass2domain)[0];
    domain_end = thetas[iclass].end();
    for (domain_iter = thetas[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      iitems = domain_iter->second.items;
      for (int i=0; i < iitems.size(); i++) {
        out(iitems[i], iclass2domain) = domain_iter->first;
      }
    }
    
  }
  
  return out;
}

float BayesParameter::domain_getloglik_x(Rcpp::IntegerVector& pattern_counts, float theta_alpha) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::domain_getloglik_x" << "\n";} // tkprint
  // log integral p(theta|domain) p(x|theta) dtheta
  
  if (pattern_counts.size() == 0) {
    return 0; // log(1)
  }
  
  return lbeta(Rcpp::as<Rcpp::NumericVector> (pattern_counts) + theta_alpha);
}

float BayesParameter::domain_getlik_domain(Hyperparameter hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::domain_getlik_domain" << "\n";} // tkprint
  // p(domain = k|other domains)
  return 1; // Assume flat prior
}

int BayesParameter::domain_id_new(int class2domain_id, Hyperparameter& hparams) {
  
  int domain_id = -1;
  int domain_class = which(hparams.class2domain == class2domain_id)[0];
  int nmax = minimum(hparams.ndomains, hparams.nitem+1);
  
  for (int i=0; i < nmax; i++) {
    if (thetas[domain_class].count(i) == 0) {
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


Rcpp::NumericVector BayesParameter::class_pi_args(Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::class_pi_args" << "\n";} // tkprint
  /**
   * Suppose alpha=R^n, pi ~ Dirchlet(alpha), classes ~ Categorical(pi)
   * Find pi|classes ~ Dirchlet(Parameters) and return Parameters
   **/
  
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

void BayesParameter::class_pi_next(Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::class_pi_next" << "\n";} // tkprint
  /**
   * Suppose alpha=R^n, pi ~ Dirchlet(alpha), classes ~ Categorical(pi)
   * Generate random sample from  pi|classes ~ Dirchlet(Parameters)
   * Gibbs sampling
   **/
  
  Rcpp::NumericVector args = class_pi_args(hparams);
  class_pi = rDirchlet(args);
}


Rcpp::NumericVector BayesParameter::class_args(Rcpp::IntegerMatrix::Row xobs) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::class_args" << "\n";} // tkprint
  /**
   * Purpose: For a single observation find the posterior classes probabilities
   **/
  
  Rcpp::NumericVector args = Rcpp::clone(class_pi);  // Prevent overwrites. Maybe instead remove & from input
  args = args * class_prob(xobs);
  // tk maybe set exact zeroes to std::numeric_limits<float>::denorm_min()
  args = args / Rcpp::sum(args);
  
  return(args);
}

void BayesParameter::classes_next(Rcpp::IntegerMatrix x) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::classes_next" << "\n";} // tkprint
  /**
   * Purpose: Gibbs sampling class
   **/
  
  int nrow = x.nrow();
  
  for (int i=0; i < nrow; i++) {
    classes(i) = rCategorical(class_args(x.row(i)));
  }
}

void BayesParameter::thetas_next(Rcpp::IntegerMatrix& x, Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::thetas_next" << "\n";} // tkprint
  /**
   * Purpose: Gibbs sampling of theta
   * Assumptions: theta_addCounts has been run providng accurate counts
   */
  
  std::map<int,  DomainCount>::iterator domain_iter;
  std::map<int,  DomainCount>::const_iterator domain_end;
  DomainCount *idomain;
  Rcpp::NumericVector iprob;
  for (int iclass=0; iclass < hparams.nclass; iclass++) {
    domain_end = thetas[iclass].end();
    for (domain_iter = thetas[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      iprob = idomain->counts + hparams.theta_alpha;
      iprob = iprob / Rcpp::sum(iprob);
      idomain->thetas = rDirchlet(iprob);
    }
  }
}

domainProposalOut BayesParameter::domain_proposal(int class2domain_id, Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::domain_proposal" << "\n";} // tkprint
  
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
  out.domain_old1 = &thetas[out.domain_classes[0]][out.domain_id1];
  if (thetas[out.domain_classes[0]].count(out.domain_id2) > 0) {
    out.domain_old2 = &thetas[out.domain_classes[0]][out.domain_id2];
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
    out.domain_new1.npatterns = 0; // tk better handling of this scenario
  }
  out.domain_new1.countReset(); // tk better handling
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
    // Rcpp::Rcout << "BayesParameter::domain_proposal:" <<"| domains_nonempty.size():" << domains_nonempty.size() << "| hparams.domain_proposal_empty:" << hparams.domain_proposal_empty << "| hparams.domain_proposal_swap:" << hparams.domain_proposal_swap << "| out.domain_old1->ndomainitems_calc():" << out.domain_old1->ndomainitems_calc() << "| out.domain_old2->ndomainitems_calc():" << out.domain_old2->ndomainitems_calc() << "| out.forwardProb:" << out.forwardProb << "\n"; // tktrouble
    
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

domainAcceptOut BayesParameter::domain_accept(Rcpp::IntegerMatrix& x, domainProposalOut& proposal, Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::domain_accept" << "\n";} // tkprint
  
  domainAcceptOut out;
  
  out.thetas_new.resize(hparams.nclass);
  
  int i;
  int iclass;
  for (i = 0; i < proposal.domain_classes.size(); i++) {
    iclass = proposal.domain_classes[i];
    if (proposal.domain_new1.ndomainitems_calc() > 0) {
      out.thetas_new[iclass][proposal.domain_id1] = proposal.domain_new1.copy();
    }
    out.thetas_new[iclass][proposal.domain_id2] = proposal.domain_new2.copy();
  } // tk move this over to proposal and only do counts here
  theta_addCounts(x, true, out.thetas_new); // tk do I really need to be true?
  
  out.loglik_old = std::log(domain_getlik_domain(hparams));
  out.loglik_new = std::log(domain_getlik_domain(hparams));
  
  for (i = 0; i < proposal.domain_classes.size(); i++) {
    
    iclass = proposal.domain_classes[i];
    
    if (thetas[iclass].count(proposal.domain_id1) > 0) { // should always be true
      out.loglik_old += domain_getloglik_x(thetas[iclass][proposal.domain_id1].counts, hparams.theta_alpha);
    }
    if (thetas[iclass].count(proposal.domain_id2) > 0) {
      out.loglik_old += domain_getloglik_x(thetas[iclass][proposal.domain_id2].counts, hparams.theta_alpha);
    }
    if (out.thetas_new[iclass].count(proposal.domain_id1) > 0) {
      out.loglik_new += domain_getloglik_x(out.thetas_new[iclass][proposal.domain_id1].counts, hparams.theta_alpha);
    }
    if (out.thetas_new[iclass].count(proposal.domain_id2) > 0) { // should always be true
      out.loglik_new += domain_getloglik_x(out.thetas_new[iclass][proposal.domain_id2].counts, hparams.theta_alpha);
    }
  }
  
  out.unif = R::runif(0, 1);
  out.log_cutoff = out.loglik_new - out.loglik_old + std::log(proposal.forwardProb) - std::log(proposal.backwardProb);
  out.accept = int( out.log_cutoff > std::log(out.unif) );
  
  return out;
}

int BayesParameter::domain_next(int class2domain_id, Rcpp::IntegerMatrix& x, Hyperparameter hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::domain_next" << "\n";} // tkprint
  // side effect. sets thetas.
  
  // tkjmb
  
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
    // domainProposalOut_print(proposal); // tktrouble 
    // domainAcceptOut_print(accept_info); // tktrouble 
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
    domain_end = accept_info.thetas_new[iclass].end();
    for (domain_iter = accept_info.thetas_new[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      thetas[iclass][domain_iter->first] = domain_iter->second; // tk confirm by reference
    }
    
    if (proposal.domain_new1.ndomainitems_calc() == 0) {
      thetas[iclass].erase(proposal.domain_id1);
    }
    if (proposal.domain_new2.ndomainitems_calc() == 0) {
      thetas[iclass].erase(proposal.domain_id2);
    }
  }
  
  return accept; // Is true
}

void BayesParameter::domains_next(Rcpp::IntegerMatrix& x, Hyperparameter& hparams) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesParameter::domains_next" << "\n";} // tkprint
  
  for (int iclass2domain=0; iclass2domain < hparams.nclass2domain; iclass2domain++) {
    for (int i=0; i < hparams.domain_nproposals; i++) {
      thetas_accept(i, iclass2domain) = domain_next(iclass2domain, x, hparams);
    }
  }
  
}


/*****************************************************
 ****** Archive
 *****************************************************/


class Archive {
  // Iterations are added by column for speed
  // tk maybe change vector to list for easier extension if we do not know maxitr
public:
  Rcpp::NumericMatrix class_pi;
  Rcpp::IntegerMatrix classes;
  std::vector<Rcpp::IntegerMatrix> thetas_id; // Row for iter, classid,  clique_id, pattern_id
  std::vector<Rcpp::IntegerMatrix> thetas_patterns; // Row per item
  std::vector<Rcpp::NumericVector> thetas_probs;
  std::vector<Rcpp::IntegerMatrix> thetas_accept;
  int next_itr;
  int maxitr;
  
public:
  void set_initial(int nclasses, int nobs, int nitems, int maxiter_in);
  void add(BayesParameter& aparams);
  void thetas2mat(BayesParameter& params, int itr, Rcpp::IntegerMatrix& out_thetas_id, Rcpp::IntegerMatrix& out_thetas_patterns, Rcpp::NumericVector& out_thetas_probs);
};

void Archive::set_initial(int nclasses, int nobs, int nitems, int maxiter_in) {
  if (IS_LOUD) {Rcpp::Rcout << "Archive::set_initial" << "\n";} // tkprint
  next_itr = 0;
  maxitr = maxiter_in;
  
  class_pi = static_cast<Rcpp::NumericMatrix>(Rcpp::no_init_matrix(nclasses, maxiter_in));
  classes = static_cast<Rcpp::IntegerMatrix>(Rcpp::no_init_matrix(nobs, maxiter_in));
  
  thetas_id.resize(maxiter_in);
  thetas_patterns.resize(maxiter_in);
  thetas_probs.resize(maxiter_in);
  thetas_accept.resize(maxiter_in);
}


void Archive::thetas2mat(BayesParameter& params, int itr, Rcpp::IntegerMatrix& out_thetas_id, Rcpp::IntegerMatrix& out_thetas_patterns, Rcpp::NumericVector& out_thetas_probs) {
  if (IS_LOUD) {Rcpp::Rcout << "Archive::thetas2mat" << "\n";} // tkprint
  
  
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
    domain_end = params.thetas[iclass].end();
    for (domain_iter = params.thetas[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      npatterns += idomain->npatterns;
    }
  }
  out_thetas_id = Rcpp::IntegerMatrix(Rcpp::no_init_matrix(4, npatterns)); // Row for iter, classid,  domain_id, pattern_id
  out_thetas_patterns = Rcpp::IntegerMatrix(Rcpp::no_init_matrix(nitems, npatterns));
  out_thetas_probs = Rcpp::NumericVector(Rcpp::no_init_vector(npatterns));
  
  
  ithis_pattern_id = 0;
  for (iclass=0; iclass < nclass; iclass++) {
    domain_end = params.thetas[iclass].end();
    for (domain_iter = params.thetas[iclass].begin(); domain_iter!=domain_end; ++domain_iter) {
      idomain = &(domain_iter->second);
      i_npatterns = idomain->npatterns;
      idomain_id = domain_iter->first;
      for (i = 0; i < i_npatterns; i++) {
        // Save info
        out_thetas_id(0, ithis_pattern_id) = itr;
        out_thetas_id(1, ithis_pattern_id) = iclass;
        out_thetas_id(2, ithis_pattern_id) = idomain_id;
        out_thetas_id(3, ithis_pattern_id) = i;
        out_thetas_patterns.column(ithis_pattern_id) = idomain->id2pattern(i);
        out_thetas_probs(ithis_pattern_id) = idomain->thetas(i);
        
        ithis_pattern_id += 1;  //Increment
      }
    }
  }
  
  
  //return std::make_tuple(out_thetas_id, out_thetas_patterns, out_thetas_probs); //std::tuple<Rcpp::IntegerMatrix, Rcpp::IntegerMatrix, Rcpp::NumericVector>
}

void Archive::add(BayesParameter& aparams) {
  if (IS_LOUD) {Rcpp::Rcout << "Archive::add" << "\n";} // tkprint
  
  if (next_itr >= maxitr) {
    Rcpp::warning("Archive::add:: Max storage reached");
    return; // Exit Early. Maybe in future resize or error, but not necessary now
  }
  
  class_pi.column(next_itr) = aparams.class_pi;
  classes.column(next_itr) = aparams.classes;
  Archive::thetas2mat(aparams, next_itr
                        , thetas_id[next_itr], thetas_patterns[next_itr], thetas_probs[next_itr]);
  thetas_accept[next_itr] = Rcpp::clone(aparams.thetas_accept);
  next_itr += 1;
}


/*****************************************************
 ****** BayesContainer
 *****************************************************/


class BayesContainer {
public:
  Hyperparameter hparams;
  BayesParameter params;
  Rcpp::IntegerMatrix x;
  Archive archive;
public:
  void set_initial(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list, int maxitr);
  void run(int niter);
  void run_once();
  void run_init(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list, int nitr);
};


void BayesContainer::set_initial(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list, int maxitr) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesContainer::set_initial" << "\n";} // tkprint
  x = x_in;
  hparams.set_hparams(hparams_list);
  hparams.set_dataInfo(x_in);
  params.set_initial(params_list, hparams);
  archive.set_initial(hparams.nclass, hparams.nobs, hparams.nitem, maxitr);
};


void BayesContainer::run(int niter) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesContainer::run" << "\n";} // tkprint
  
  for (int i=0; i < niter; i++) {
    if (IS_LOUD) {Rcpp::Rcout << "BayesContainer::run" << "| i" << i << "\n";} // tkprint
    // tk maybe optionally randomize order or add/omit steps
    params.class_pi_next(hparams);
    params.classes_next(x);
    params.theta_addCounts(x, true);
    params.domains_next(x, hparams); // tkjmb
    params.thetas_next(x, hparams);
    archive.add(params);
  }
};

void BayesContainer::run_init(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list
                                , int nitr) {
  if (IS_LOUD) {Rcpp::Rcout << "BayesContainer::run_init" << "\n";} // tkprint
  
  //tk optionallly set thetas automatically with gibbs??
  set_initial(x_in, hparams_list, params_list, nitr);
  run(nitr);
}


/*****************************************************
 ****** Public Functions
 *****************************************************/


// [[Rcpp::export]]
Rcpp::List dependentLCM_fit_cpp(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list
                                  , int nitr) {
  if (IS_LOUD) {Rcpp::Rcout << "dependentLCM_fit" << "\n";} // tkprint
  BayesContainer bcontainer;
  bcontainer.run_init(x_in, hparams_list, params_list, nitr);
  return Rcpp::List::create(
    Rcpp::Named("class_pi") = bcontainer.archive.class_pi
    , Rcpp::Named("classes") = bcontainer.archive.classes
    , Rcpp::Named("thetas_id") = wrap(bcontainer.archive.thetas_id)
    , Rcpp::Named("thetas_patterns") = wrap(bcontainer.archive.thetas_patterns)
    , Rcpp::Named("thetas_probs") = wrap(bcontainer.archive.thetas_probs)
    , Rcpp::Named("next_itr") = bcontainer.archive.next_itr
    , Rcpp::Named("maxitr") = bcontainer.archive.maxitr
    , Rcpp::Named("thetas_accept") = bcontainer.archive.thetas_accept
  );
}
