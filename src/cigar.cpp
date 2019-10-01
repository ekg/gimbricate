#include "cigar.hpp"

namespace gimbricate {

std::vector<std::pair<uint64_t, char> > split_cigar(const std::string& cigar_str) {
    std::vector<std::pair<uint64_t, char> > cigar;
    std::string number;
    char type = '\0';
    // strings go [Number][Type] ...
    for (auto& c : cigar_str) {
        if (std::isdigit(c)) {
            if (type == '\0') {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(std::make_pair(atoi(number.c_str()), type));
                number.clear();
                type = '\0';
                number += c;
            }
        } else {
            type += c;
        }
    }
    if (!number.empty() && type != '\0') {
        cigar.push_back(std::make_pair(std::atoi(number.c_str()), type));
    }
    return cigar;
}

uint64_t cigar_length(const std::vector<std::pair<uint64_t, char>>& cigar) {
    uint64_t len = 0;
    for (auto& elem : cigar) {
        len += elem.first;
    }
    return len;
}

std::string graph_cigar(gssw_graph_mapping* gm) {
    std::stringstream s;
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int to_pos = 0;
    int from_pos = gm->position;
    s << from_pos << '@';
    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (int j=0; j < l; ++j, ++e) {
            s << e->length << e->type;
        }
    }
    return s.str();
}

}
