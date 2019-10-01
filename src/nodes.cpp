#include "nodes.hpp"

namespace gimbricate {

node_index::node_index (const std::vector<std::string>& data) {
	uint nthreads = get_thread_count();
	// mphf takes as input a c++ range. A std::vector is already a c++ range
	double gammaFactor = 2.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )
	//build the mphf
	bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(data.size(),data,nthreads,gammaFactor);
}

node_index::~node_index(void) {
    delete bphf;
}

uint64_t node_index::get_id(const std::string& key) {
    return bphf->lookup(key);
}

}
