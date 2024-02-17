// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <gtl/include/gtl/phmap.hpp>

using namespace Rcpp;

// mappings
typedef gtl::flat_hash_map<std::string, int> KmerIndexMap;

// [[Rcpp::export]]
std::vector<int> get_kmer_indices(
  const std::vector<std::string> &all_kmers, 
  const std::vector<std::string> &kmers_extracted){
  
  // Create the hashmap
  KmerIndexMap kmer_index_map;
  for(int i = 0; i < all_kmers.size(); i++){
      kmer_index_map[all_kmers[i]] = i + 1;  // map each k-mer to its index (1-based)
  }
  
  // Look up the indices of kmers_extracted
  std::vector<int> indices(kmers_extracted.size());
  for(int i = 0; i < kmers_extracted.size(); i++){
      auto it = kmer_index_map.find(kmers_extracted[i]);
      if(it != kmer_index_map.end()){
          indices[i] = it->second;
      } else {
          indices[i] = NA_INTEGER;  // return NA if the k-mer is not found
      }
  }
  
  return indices;
}