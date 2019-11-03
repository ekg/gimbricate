#include <vector>
#include <string>
#include <utility>
#include <cctype>
#include <sstream>
#include <algorithm>
#include "gssw.h"

namespace gimbricate {

std::vector<std::pair<uint64_t, char>> split_cigar(const std::string& cigar_str);
uint64_t cigar_length(const std::vector<std::pair<uint64_t, char>>& cigar);
std::string compress_cigar(const std::string& cigar);
std::string flip_cigar(const std::string& cigar);
std::string overlap_cigar(gssw_graph_mapping* gm);
std::string simple_cigar(gssw_graph_mapping* gm);
void mapping_boundaries(gssw_graph_mapping* gm,
                        uint64_t& query_start, uint64_t& query_end,
                        uint64_t& target_start, uint64_t& target_end,
                        uint64_t& num_matches);


}
