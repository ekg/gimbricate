#include <vector>
#include <string>
#include <utility>
#include <cctype>
#include <sstream>
#include "gssw.h"

namespace gimbricate {

std::vector<std::pair<uint64_t, char>> split_cigar(const std::string& cigar_str);
uint64_t cigar_length(const std::vector<std::pair<uint64_t, char>>& cigar);
std::string graph_cigar(gssw_graph_mapping* gm);

}
