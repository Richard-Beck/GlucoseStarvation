#include <Rcpp.h>
using namespace Rcpp;

// A little struct to hold per‐frame data
struct FrameData {
  int    n;               // # objects
  double sum_max_prob;    // sum of pmax(...)
  double sum_junk;        // total prob_junk
  std::vector<double> pd, pa, pj;        // dead, alive, junk probs
  std::vector<int>     dead_idx, alive_idx; // sorted orders
};

// We store all frames in one vector, keyed by unique time.
// We wrap it in an XPtr so Rcpp can hold it across calls.
typedef std::vector<FrameData> Preproc;

// [[Rcpp::export]]
XPtr<Preproc> preprocess_dt(
    NumericVector time,
    NumericVector prob_alive,
    NumericVector prob_dead,
    NumericVector prob_junk
) {
  int N = time.size();
  // 1) figure out unique times and an index for each row
  NumericVector times = clone(time);
  std::vector<double> uniq = as< std::vector<double> >(times);
  std::sort(uniq.begin(), uniq.end());
  uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
  int F = uniq.size();
  
  // map time -> frame#
  std::unordered_map<double,int> idx_map;
  idx_map.reserve(F);
  for(int i=0;i<F;i++) idx_map[uniq[i]] = i;
  
  // allocate
  Preproc *out = new Preproc(F);
  // first pass: count objects per frame
  std::vector< std::vector<int> > rows(F);
  for(int i=0;i<N;i++){
    int f = idx_map[ time[i] ];
    rows[f].push_back(i);
  }
  
  // fill each FrameData
  for(int f=0; f<F; f++){
    FrameData &fd = (*out)[f];
    fd.n = rows[f].size();
    fd.pd.resize(fd.n);
    fd.pa.resize(fd.n);
    fd.pj.resize(fd.n);
    fd.dead_idx.resize(fd.n);
    fd.alive_idx.resize(fd.n);
    
    // gather
    for(int j=0; j<fd.n; j++){
      int r = rows[f][j];
      fd.pd[j] = prob_dead[r];
      fd.pa[j] = prob_alive[r];
      fd.pj[j] = prob_junk[r];
      fd.dead_idx[j] = j;
      fd.alive_idx[j] = j;
    }
    // compute sum_max_prob and sum_junk
    fd.sum_max_prob = 0;
    fd.sum_junk     = 0;
    for(int j=0;j<fd.n;j++){
      fd.sum_max_prob += std::max({ fd.pd[j], fd.pa[j], fd.pj[j] });
      fd.sum_junk     += fd.pj[j];
    }
    // sort dead_idx by pd desc
    std::sort(fd.dead_idx.begin(), fd.dead_idx.end(),
              [&](int a,int b){ return fd.pd[a] > fd.pd[b]; });
    // sort alive_idx by pa desc
    std::sort(fd.alive_idx.begin(), fd.alive_idx.end(),
              [&](int a,int b){ return fd.pa[a] > fd.pa[b]; });
  }
  
  // store unique times as attribute, so we can map later
  XPtr<Preproc> p(out, true);
  p.attr("times") = uniq;
  return p;
}

// simple n^2 cumulative‐penalty
static double cum_penalty(const IntegerVector &v){
  int n=v.size();
  double tot=0;
  for(int i=0;i<n-1;i++){
    for(int j=i+1;j<n;j++){
      double d = v[i] - v[j];
      if(d>0) tot += d;
    }
  }
  return tot;
}

// [[Rcpp::export]]
double calculate_total_cost(
    XPtr<Preproc> dt_ref,
    NumericVector target_time,
    IntegerVector N_target,
    IntegerVector D_target,
    double lambda_N_mono=1.0,
    double lambda_D_mono=1.0,
    double lambda_create=1.0
) {
  // --- 1. Monotonicity penalty ---
  int M = target_time.size();
  // ensure sorted by time
  IntegerVector idx = seq_len(M) - 1;
  std::sort(idx.begin(), idx.end(),
            [&](int a,int b){ return target_time[a] < target_time[b]; });
  IntegerVector Nvec(M), Dvec(M);
  for(int i=0;i<M;i++){
    Nvec[i] = N_target[ idx[i] ];
    Dvec[i] = D_target[ idx[i] ];
  }
  double cost_M = lambda_N_mono * cum_penalty(Nvec)
    + lambda_D_mono * cum_penalty(Dvec);
  
  // build map time->frame#
  NumericVector uniq = dt_ref.attr("times");
  int F = uniq.size();
  std::unordered_map<double,int> tmap;
  tmap.reserve(F);
  for(int i=0;i<F;i++) tmap[ uniq[i] ] = i;
  
  // --- 2. per-frame costs ---
  double sum_cost_D = 0, sum_cost_C = 0;
  for(int i=0;i<M;i++){
    double t = target_time[i];
    auto it = tmap.find(t);
    if(it==tmap.end()) {
      // if no data for that time: creation penalty = lambda_create*N_target, D and C accordingly
      sum_cost_C += lambda_create * N_target[i];
      continue;
    }
    FrameData &fd = (*dt_ref)[ it->second ];
    int n_tot = fd.n;
    int Nt = N_target[i];
    int Dt = D_target[i];
    int At = Nt - Dt;
    // creation
    double cost_C = lambda_create * std::max(0, Nt - n_tot);
    
    // select top-Dt by pd
    double sumD=0, sumA=0, sumJsel=0;
    std::vector<bool> used(n_tot,false);
    for(int k=0;k<std::min(Dt,n_tot);k++){
      int r = fd.dead_idx[k];
      sumD += fd.pd[r];
      sumJsel += fd.pj[r];
      used[r]=true;
    }
    // select top-At by pa among unused
    int gotA=0;
    for(int id: fd.alive_idx){
      if(gotA>=At) break;
      if(!used[id]){
        sumA += fd.pa[id];
        sumJsel += fd.pj[id];
        used[id]=true;
        gotA++;
      }
    }
    // junk remaining is fd.sum_junk - sumJsel
    double sumJrem = fd.sum_junk - sumJsel;
    double achieved = sumD + sumA + sumJrem;
    double cost_D = fd.sum_max_prob - achieved;
    
    sum_cost_C += cost_C;
    sum_cost_D += cost_D;
  }
  
  return cost_M + sum_cost_C + sum_cost_D;
}

