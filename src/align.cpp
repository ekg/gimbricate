#include "align.hpp"
#include <iostream>
#include <fstream>

using std::endl;
namespace gimbricate {

std::string align_ends(const std::string& seq_x_full, const std::string& seq_y_full,
                       const uint64_t& length,
                       bool seq_x_forward, bool seq_y_forward,
                       uint64_t& query_start, uint64_t& query_end,
                       uint64_t& target_start, uint64_t& target_end,
                       uint64_t& num_matches, std::string& basic_cigar) {

    // default parameters for genome sequence alignment
    std::string seq_x;
    uint64_t seq_x_offset = 0;
    if (length > seq_x_full.size()) {
        seq_x = seq_x_full;
    } else {
        seq_x_offset = seq_x_full.size()-length;
        seq_x = seq_x_full.substr(seq_x_offset, length);
    }
    std::string seq_y;
    uint64_t seq_y_offset = 0;
    if (length > seq_y_full.size()) {
        seq_y = seq_y_full;
    } else {
        seq_y = seq_y_full.substr(0, length);
    }
    
    //gssw_print_graph_mapping(gm, stdout);
    basic_cigar = std::to_string(length)+'M';

    query_start = 0;
    query_end = length;
    target_start = seq_x_full.length()-length;
    target_end = seq_x_full.length();
    num_matches = length;

    // correct back to the coordinates of the full sequence
//    query_start += seq_y_offset;
//    query_end += seq_y_offset;
//    target_start += seq_x_offset;
//    target_end += seq_x_offset;

    // then invert the coordinates in the case that the pieces are inverted
    if (!seq_y_forward) {
        uint64_t i = seq_y_full.length() - query_start;
        uint64_t j = seq_y_full.length() - query_end;
        query_start = j;
        query_end = i;
    }
    if (!seq_x_forward) {
        uint64_t i = seq_x_full.length() - target_start;
        uint64_t j = seq_x_full.length() - target_end;
        target_start = j;
        target_end = i;
    }

    return std::to_string(length)+'M';
}

}
