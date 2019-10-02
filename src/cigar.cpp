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

std::string compress_cigar(const std::string& cigar) {
    std::vector<std::pair<uint64_t, char>> compressed;
    for (auto& elem : split_cigar(cigar)) {
        if (compressed.empty()) {
            compressed.push_back(elem);
        } else {
            if (elem.second == compressed.back().second) {
                compressed.back().first += elem.first;
            } else {
                compressed.push_back(elem);
            }
        }
    }
    std::stringstream s;
    for (auto& elem : compressed) {
        s << elem.first << elem.second;
    }
    return s.str();
}

std::string overlap_cigar(gssw_graph_mapping* gm) {
    std::stringstream s;
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    //int to_pos = 0;
    int from_pos = gm->position;
    if (from_pos) s << from_pos << 'D';
    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (int j=0; j < l; ++j, ++e) {
            char type = e->type;
            switch (type) {
            case 'S': type = 'I'; break;
            case 'X': type = 'M'; break;
            default: break;
            }
            // trim leading insertions as these just change the overlap length
            // checking that from_pos is 0 ensures that we didn't already write D into the cigar
            if (!from_pos && type == 'I') continue;
            s << e->length << type;
        }
    }
    return compress_cigar(s.str());
}

}
