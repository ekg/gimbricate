#pragma once

#include <string>
#include <iostream>
#include "edlib.h"
#include "cigar.hpp"

namespace gimbricate {

/// align the last l bp of x to the start of y, returning the cigar for this region
std::string align_ends(const std::string& seq_x_full,
                       const std::string& seq_y_full,
                       const uint64_t& length,
                       bool seq_x_forward, bool seq_y_forward,
                       uint64_t& query_start, uint64_t& query_end,
                       uint64_t& target_start, uint64_t& target_end,
                       uint64_t& num_matches,
                       std::string& basic_cigar,
                       bool perfectOverlaps);


}
