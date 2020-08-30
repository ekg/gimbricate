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
                       uint64_t& num_matches, std::string& basic_cigar, bool perfectOverlaps) {

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
    
    if (perfectOverlaps) {

        basic_cigar = std::to_string(length)+'M';

        query_start = 0;
        query_end = length;
        target_start = seq_x_full.length()-length;
        target_end = seq_x_full.length();
        num_matches = length;

        // invert the coordinates in the case that the pieces are inverted
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

    else {

        EdlibAlignResult result = edlibAlign(seq_y.c_str(), seq_y.length(),
                                             seq_x.c_str(), seq_x.length(),
                                             edlibNewAlignConfig(-1, // no edit limit
                                                                 EDLIB_MODE_NW, // global alignment
                                                                 EDLIB_TASK_PATH,
                                                                 NULL, 0));
        if (result.status != EDLIB_STATUS_OK) {
            std::cerr << "[gimbricate::align] error: alignment failure" << std::endl
                      << "query\t" << seq_y << std::endl
                      << "target\t" << seq_x << std::endl;
            exit(1);
        }

        // get the cigar
        // get the starting and ending locations
        // boom
        uint64_t target_aligned_length=0, query_aligned_length=0;
        uint64_t mismatches=0, insertions=0, deletions=0, softclips=0;

        char* cigar = edlib_alignment_to_cigar(result.alignment,
                                               result.alignmentLength,
                                               target_aligned_length,
                                               query_aligned_length,
                                               num_matches,
                                               mismatches,
                                               insertions,
                                               deletions,
                                               softclips);

        query_start = seq_y_offset;
        query_end = query_aligned_length + seq_y_offset;
        target_start = seq_x_offset + result.startLocations[0];
        target_end = seq_x_offset + result.endLocations[0];

        //printf("edit_distance('hello', 'world!') = %d\n", result.editDistance);
        edlibFreeAlignResult(result);

        std::string final_cigar = std::string(cigar);
        free(cigar);
        basic_cigar = final_cigar;

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

        return final_cigar;
    }


}

}
