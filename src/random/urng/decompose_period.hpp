#ifndef flopoco_random_decompose_period_hpp
#define flopoco_random_decompose_period_hpp

#include <gmpxx.h>

#include <vector>
#include <stdint.h>
#include <limits.h>
#include <algorithm>
#include <iostream>

namespace flopoco
{
namespace random
{

/* NOTE : I've come to the conclusion this is a fairly flawed way of doing things, as
  it is not getting either the maximum period for the ensemble, or trying to maximise
  the minimum period in the ensemble. However, it's still much better to have very low
  pair-wise period GCDs, so for now we'll stick with it.

/*! We often want to generate some huge number of bits (e.g. r>=1000) , and may not have the
  correct components. We probably also don't want a connected structure of that size anyway.
  Using parallel generators of the same period is not a good idea, so what we'll do is select
  a set of distinct generators G_1..G_k, such that sum(r_1)...r_k) >= r, but also such that
  forall(i=1..k, n_i>=n_min and forall(j=1..k, gcd(n_i,n_j) >= p ) )
  for some small p and some reasonable n_min.

  If we are able to choose p=1, then the period of the ensemble generator would be prod(n_i^2-1,i=1..k),
  but there are often gaps in the table, making it difficult to find such solutions. In general for p>=1, then
  in the worst case the period is prod(n_i^2-1,i=1..k) / (p^k). This means the overall generator loses a bit of
  period, but is _much_ better than the overall generator only have period 2^n_min-1.
  
  It also allows us to cheat with seeding, as one max(n_1..n_k) bit seed can safely be used for
  all the generators, so we can retain the two-bit seeding interface.

  For practical purposes I've tended to enumerate even n when making tables, so we often
  need p>=3 or p>=7 to get reasonable solutions.

  \todo A nicety would be to find the solution with minimal lcm(p_1..p_k), but that doesn't
  seem warranted.
*/
class decompose_period{
  // Vector of (r,n) which are available.
  std::vector<std::pair<int,int> > m_generators;
  
  int m_maxN;   // Largest value of n available
  std::vector<int> m_gcds;  // Cache of gcd calculations
  
  int m_targetR;   // The overall number of output bits we would like
  int m_p;              // Maximum pair-wise gcd
  
  //! There are probably faster tricks for this, but this just does gcd(2^n1-1,2^n2-1)
  int GCD(int n1, int n2)
  {
    int key=m_maxN*n1+n2;
    int cached=m_gcds[key];
    
    if(cached==0){
      mpz_class p1=(mpz_class(1)<<n1)-1, p2=(mpz_class(1)<<n2)-1;
      mpz_class g;
      mpz_gcd(g.get_mpz_t(), p1.get_mpz_t(), p2.get_mpz_t());
      if(g>=INT_MAX)
        cached=INT_MAX;
      else
        cached=g.get_si();
    }
    
    return cached;
  }
  
  //! Find next generator g such that n>=n_min and forall(i=1..curr.size(), gcd(2^curr[0].first-1,2^n-1) >=p). If no such generator, return -1
  int FindNextG(const std::vector<int> &curr, int begin_g, int end_g){
    int g=begin_g;
    while(g<end_g){
      
      int n=m_generators[g].second;
      bool ok=true;
      for(int i=0;i<(int)curr.size();i++){
        if(GCD(curr[i], n)>m_p){
          ok=false;
          break;
        }
      }
      if(ok)
        return g;
      g++;
    }
    return -1;
  }
  
  int m_bestTotalR;
  std::vector<int> m_bestSolution;

  bool BuildSolution(std::vector<int> &curr, int total_r)
  {
    bool found=false;
    if(total_r >= m_targetR){
      if(total_r < m_bestTotalR){
        m_bestSolution=curr;
        m_bestTotalR=total_r;
        std::cerr<<"bestTotalR="<<m_bestTotalR<<"\n";
        found=m_bestTotalR==m_targetR;
      }
    }else{
      int g=curr.size()==0 ? 0 : curr.back();
      while(true){
        g=FindNextG(curr, g+1, m_generators.size());
        if(g==-1)
          return found;
        curr.push_back(g);
        found=found || BuildSolution(curr, total_r+m_generators[g].first);
        curr.pop_back();
        if(m_bestTotalR==m_targetR)
          return found;
      }
    }
    return found;
  }

  bool Search()
  {
    std::vector<int> working;
    return BuildSolution(working, 0);
  }
 
  decompose_period(int target_r,  const std::vector<std::pair<int,int> > &generators)
    : m_generators(generators)
    , m_targetR(target_r)
    , m_bestTotalR(INT_MAX)
  {    
    mpz_class pp(1);
    
    std::sort(m_generators.begin(), m_generators.end());
    std::reverse(m_generators.begin(), m_generators.end());
    
    while(true){
      m_p=pp.get_si();
      
      m_maxN=1;
      for(int i=0;i<(int)m_generators.size();i++){
        m_maxN=std::max(m_maxN, m_generators[i].second+1);
      }
      m_gcds.resize(m_maxN*m_maxN);
      
      if(Search())
        return;
      
      mpz_nextprime(pp.get_mpz_t(), pp.get_mpz_t());
      if(pp>=13){
        if(m_bestTotalR!=INT_MAX)
          return;
      }
        
    }
  }
public:
  static std::vector<std::pair<int,int> > Execute(int target_r,   const std::vector<std::pair<int,int> > &generators)
  {
    decompose_period d(target_r, generators);
    std::vector<std::pair<int,int> > res(d.m_bestSolution.size());
    for(int i=0;i<res.size();i++){
      res[i]=d.m_generators[d.m_bestSolution[i]];
    }
    return res;
  }
  
  static std::pair<mpz_class,mpz_class> period(const std::vector<std::pair<int,int> > & parts)
  {
    mpz_class p_max=1, p_got=1;
      for(int i=0;i<parts.size();i++){
        p_max=(p_max<<parts[i].second)-p_max;
        mpz_class curr=(mpz_class(1)<<parts[i].second)-1;
        mpz_lcm(p_got.get_mpz_t(), p_got.get_mpz_t(), curr.get_mpz_t());
      }
    return std::make_pair(p_max,p_got);
  }
};

}; // random
}; // flopoco

#endif
