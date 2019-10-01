#pragma once

#include <string>
#include "gssw.h"

namespace gimbricate {

/// align the last l bp of x to the start of y, returning the cigar for this region
std::string align_ends(const std::string& seq_x, const std::string& seq_y, const uint64_t& length);

}
