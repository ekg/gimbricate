#include <vector>
#include <string>
#include <cstring>
#include <utility>
#include <cctype>
#include <sstream>
#include <algorithm>

namespace gimbricate {

std::vector<std::pair<uint64_t, char>> split_cigar(const std::string& cigar_str);
uint64_t cigar_length(const std::vector<std::pair<uint64_t, char>>& cigar);
std::string compress_cigar(const std::string& cigar);
std::string flip_cigar(const std::string& cigar);

char* edlib_alignment_to_cigar(const unsigned char* const alignment,
                               const int alignmentLength,
                               uint64_t& refAlignedLength,
                               uint64_t& qAlignedLength,
                               uint64_t& matches,
                               uint64_t& mismatches,
                               uint64_t& insertions,
                               uint64_t& deletions,
                               uint64_t& softclips);

}
