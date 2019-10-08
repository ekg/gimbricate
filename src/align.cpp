#include "align.hpp"

namespace gimbricate {

std::string align_ends(const std::string& seq_x_full, const std::string& seq_y_full, const uint64_t& length) {
    // default parameters for genome sequence alignment
    std::string seq_x = (length > seq_x_full.size() ? seq_x_full : seq_x_full.substr(seq_x_full.size()-length, length));
    std::string seq_y = seq_y_full.substr(0, length);
    
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

    int8_t fl_bonus = 63;
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
    gssw_graph_mapping_destroy(gm);
    
    // note that nodes which are referred to in this graph are destroyed as well
    gssw_graph_destroy(graph);

    free(nt_table);
	free(mat);

    return final_cigar;
}

}
