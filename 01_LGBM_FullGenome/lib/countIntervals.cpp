// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>

using namespace Rcpp;

/**
 * Calculate number of breaks in a genomic interval.
 * @param bins_start start of bin interval
 * @param bins_end end of bin interval
 * @param start_bp start of feature interval
 * @param end_bp end of feature interval
 * @return Sum of breaks in each interval.
*/
// [[Rcpp::export]]
std::vector<int> countIntervalForFeatures(
    std::vector<int> &bins_start,
    std::vector<int> &bins_end,
    std::vector<int> &start_bp,
    std::vector<int> &end_bp
    ){
    int n_res = bins_start.size();
    int n_bp = start_bp.size();

    // sort based on start_bp and then end_bp for the two-pointer approach
    std::vector<int> indices(n_bp);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(
        indices.begin(), indices.end(),
        [&](int a, int b){
            return start_bp[a] < start_bp[b] || 
            (start_bp[a] == start_bp[b] && end_bp[a] < end_bp[b]);
        }
    );
    
    // pre-allocate space for count vector
    std::vector<int> count(n_res);

    int start_ptr = 0;
    // using a two-pointer approach to efficiently count
    for(int i = 0; i < n_res; i++){
        // move the start pointer to the first breakpoint whose end >= bins_start[i]
        // purpose: move the start_ptr forward, skipping any breakpoints that are 
        //          entirely before the current start of the bin. The start_ptr stops  
        //          moving when it hits a breakpoint that might overlap with the  
        //          start of the bin or lies completely within it.
        while(start_ptr < n_bp && end_bp[indices[start_ptr]] < bins_start[i]){
            start_ptr++;
        }

        int end_ptr = start_ptr;

        // move the end pointer while its start is <= bins_end[i]
        while(end_ptr < n_bp && start_bp[indices[end_ptr]] <= bins_end[i]){
            end_ptr++;
        }

        count[i] = end_ptr - start_ptr;
    }

    return count;
}

/**
 * Calculate number of breaks in a genomic interval.
 * @param bins_start start of bin interval
 * @param bins_end end of bin interval
 * @param start_bp start of break interval
 * @return Sum of breaks in each interval.
*/
// [[Rcpp::export]]
std::vector<int> countIntervalForBreaks(
    std::vector<int> &bins_start,
    std::vector<int> &bins_end,
    std::vector<int> &start_bp
    ){
    int n_res = bins_start.size();
    int n_bp = start_bp.size();

    // sort start_bp for the two-pointer approach
    std::sort(start_bp.begin(), start_bp.end());

    // pre-allocate space for count vector
    std::vector<int> count(n_res);

    int start_ptr = 0;
    // using a two-pointer approach to efficiently count
    for(int i = 0; i < n_res; i++){
        int end_ptr = start_ptr;

        // move the start pointer to the first breakpoint >= bins_start[i]
        while(start_ptr < n_bp && start_bp[start_ptr] < bins_start[i]){
            start_ptr++;
        }

        // ove the end pointer to the first breakpoint > bins_end[i]
        while(end_ptr < n_bp && start_bp[end_ptr] <= bins_end[i]){
            end_ptr++;
        }

        count[i] = end_ptr - start_ptr;
    }

    return count;
}