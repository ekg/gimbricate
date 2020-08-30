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

std::string flip_cigar(const std::string& cigar) {
    auto cigar_split = split_cigar(cigar);
    std::reverse(cigar_split.begin(), cigar_split.end());
    std::stringstream s;
    for (auto& c : cigar_split) {
        s << c.first << c.second;
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
            //if (!from_pos && type == 'I') continue;
            s << e->length << type;
        }
    }
    //return compress_cigar(s.str());
    return s.str(); //compress_cigar(s.str());
}

std::string simple_cigar(gssw_graph_mapping* gm) {
    std::stringstream s;
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    for (int i = 0; i < gc->length; ++i, ++nc) {
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (int j=0; j < l; ++j, ++e) {
            if (e->type == 'S') continue;
            s << e->length << e->type;
        }
    }
    return compress_cigar(s.str());
}

void mapping_boundaries(gssw_graph_mapping* gm,
                        uint64_t& query_start, uint64_t& query_end,
                        uint64_t& target_start, uint64_t& target_end,
                        uint64_t& num_matches) {
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int64_t to_pos = 0;
    query_start = 0; // wait for our first soft clip to decide where we start
    int64_t from_pos = gm->position;
    target_start = from_pos;
    for (int64_t i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first,
        gssw_cigar* c = nc->cigar;
        uint64_t l = c->length;
        gssw_cigar_element* e = c->elements;
        for (uint64_t j=0; j < l; ++j, ++e) {
            char type = e->type;
            switch (type) {
            case 'I':
            case 'S':
                if (to_pos == 0) {
                    to_pos += e->length;
                    query_start = to_pos;
                }
                break;
            case 'D': from_pos += e->length; break;
            case 'M':
            case 'X':
                from_pos += e->length;
                to_pos += e->length;
                num_matches += e->length;
                break;
            default: break;
            }
        }
    }
    query_end = to_pos;
    target_end = from_pos;
}

char* edlib_alignment_to_cigar(const unsigned char* const alignment,
                               const int alignmentLength,
                               uint64_t& refAlignedLength,
                               uint64_t& qAlignedLength,
                               uint64_t& matches,
                               uint64_t& mismatches,
                               uint64_t& insertions,
                               uint64_t& deletions,
                               uint64_t& softclips) {

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};

    std::vector<char>* cigar = new std::vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
            case 'I':
                // assume that starting and ending insertions are softclips
                if (i == alignmentLength || cigar->empty()) {
                    softclips += numOfSameMoves;
                } else {
                    insertions += numOfSameMoves;
                }
                qAlignedLength += numOfSameMoves;
                break;
            case '=':
                matches += numOfSameMoves;
                qAlignedLength += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            case 'X':
                mismatches += numOfSameMoves;
                qAlignedLength += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            case 'D':
                deletions += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            default:
                break;
            }
                  
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    delete cigar;
                    return 0;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0);  // Null character termination.
    char* cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}

}
