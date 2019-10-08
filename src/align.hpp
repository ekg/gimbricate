#pragma once

#include <string>
#include <iostream>
#include "gssw.h"
#include "cigar.hpp"

namespace gimbricate {

/// align the last l bp of x to the start of y, returning the cigar for this region
std::string align_ends(const std::string& seq_x_full, const std::string& seq_y_full, const uint64_t& length);

}
