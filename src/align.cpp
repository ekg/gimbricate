#include "align.hpp"

namespace gimbricate {

std::string align_ends(const std::string& seq_x, const std::string& seq_y, const uint64_t& length) {
        // default parameters for genome sequence alignment
    int8_t match = 1, mismatch = 4;
    uint8_t gap_open = 6, gap_extension = 1;
    // from Mengyao's example about the importance of using all three matrices in traceback.
    // int32_t l, m, k, match = 2, mismatch = 1, gap_open = 2, gap_extension = 1;

    gssw_sse2_disable();

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

    gssw_node* nodes[1];
    nodes[0] = (gssw_node*)gssw_node_create(NULL, 1, seq_y.c_str(), nt_table, mat);
    
    gssw_graph* graph = gssw_graph_create(1);
    //memcpy((void*)graph->nodes, (void*)nodes, 1*sizeof(gssw_node*));
    graph->size = 1;
    gssw_graph_add_node(graph, nodes[0]);
    
    gssw_graph_fill(graph, seq_x.c_str(), nt_table, mat, gap_open, gap_extension, 10, 10, 15, 2, true);
    
    gssw_node** pinning_nodes = (gssw_node**) malloc(sizeof(gssw_node*));
    pinning_nodes[0] = nodes[0];
    
    gssw_graph_mapping* gmp = gssw_graph_trace_back_pinned (graph,
                                                            seq_x.c_str(),
                                                            seq_x.size(),
                                                            pinning_nodes,
                                                            1,
                                                            nt_table,
                                                            mat,
                                                            gap_open,
                                                            gap_extension,
                                                            10, 10);
    
    printf("Optimal pinned mapping:\n");
    gssw_print_graph_mapping(gmp, stdout);
    gssw_graph_mapping_destroy(gmp);
    
    free(pinning_nodes);
    
    // note that nodes which are referred to in this graph are destroyed as well
    gssw_graph_destroy(graph);

    free(nt_table);
	free(mat);

    return "";
}

}
