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

std::string remove_soft_clips(const std::string& cigar) {
    auto cigar_split = split_cigar(cigar);
    std::stringstream s;
    for (auto& c : cigar_split) {
        if (c.second != 'S') {
            s << c.first << c.second;
        }
    }
    return s.str();
}

uint64_t count_matches(const std::string& cigar) {
    auto cigar_split = split_cigar(cigar);
    uint64_t matches = 0;
    for (auto& c : cigar_split) {
        if (c.second == '=') {
            matches += c.first;
        }
    }
    return matches;
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
