#include "align.hpp"

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
        seq_y_offset = 0;
        seq_y = seq_y_full.substr(0, length);
    }
    
    int8_t match = 1, mismatch = 4;
    uint8_t gap_open = 6, gap_extension = 1;

	/* This table is used to transform nucleotide letters into numbers. */
    int8_t* nt_table = gssw_create_nt_table();
    
	// initialize scoring matrix for genome sequences
	//  A  C  G  T	N (or other ambiguous code)
	//  2 -2 -2 -2 	0	A
	// -2  2 -2 -2 	0	C
	// -2 -2  2 -2 	0	G
	// -2 -2 -2  2 	0	T
	//	0  0  0  0  0	N (or other ambiguous code)
    int8_t* mat = gssw_create_score_matrix(match, mismatch);

    gssw_node* node = (gssw_node*)gssw_node_create(NULL, 1, seq_x.c_str(), nt_table, mat);
    gssw_graph* graph = gssw_graph_create(1);
    gssw_graph_add_node(graph, node);

    int8_t fl_bonus = 5; // try 63 to force full length alignment
    gssw_graph_fill(graph, seq_y.c_str(), nt_table, mat, gap_open, gap_extension, fl_bonus, fl_bonus, 15, 2, true);

    gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                    seq_y.c_str(),
                                                    seq_y.length(),
                                                    nt_table,
                                                    mat,
                                                    gap_open,
                                                    gap_extension,
                                                    fl_bonus, fl_bonus);

    
    //gssw_print_graph_mapping(gm, stdout);
    std::string final_cigar = overlap_cigar(gm);
    basic_cigar = simple_cigar(gm);

    // determine the start and end of the alignment in the target and query
    mapping_boundaries(gm,
                       query_start, query_end,
                       target_start, target_end,
                       num_matches);

    // correct back to the coordinates of the full sequence
    query_start += seq_y_offset;
    query_end += seq_y_offset;
    target_start += seq_x_offset;
    target_end += seq_x_offset;

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

    // clean up
    gssw_graph_mapping_destroy(gm);

    // note that nodes which are referred to in this graph are destroyed as well
    gssw_graph_destroy(graph);

    free(nt_table);
	free(mat);

    return final_cigar;
}

}
