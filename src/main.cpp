#include <iostream>
#include <cstdint>
#include <cstdio>
#include <string>
#include "exists.hpp"
#include "threads.hpp"
#include "args.hxx"
#include "gssw.h"
#include "gfakluge.hpp"
#include "nodes.hpp"
#include "cigar.hpp"
#include "dna.hpp"
#include "align.hpp"
#include "paf.hpp"
//#include <limits>

using namespace gimbricate;

int main(int argc, char** argv) {
    args::ArgumentParser parser("gimbricate: recompute GFA link overlaps");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in_file(parser, "FILE", "use this GFA FILE as input", {'g', "gfa-in"});
    //args::ValueFlag<std::string> gfa_out_file(parser, "FILE", "write transformed GFA to FILE", {'o', "gfa-out"});
    args::ValueFlag<std::string> fasta_out_file(parser, "FILE", "write renamed sequences to FASTA FILE", {'f', "fasta-out"});
    args::ValueFlag<std::string> paf_out_file(parser, "FILE", "write GFA overlap alignments to this PAF FILE", {'p', "paf-out"});
    args::ValueFlag<uint64_t> read_cov_min(parser, "N", "require this many supporting reads in the RC tag to keep a node", {'c', "read-coverage"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<double> expand_align(parser, "N", "expand the alignment length by this ratio (default 2.0)", {'e', "expand-align"});
    args::Flag no_rename(parser, "no-rename", "don't rename sequences to have increasing non-0 integer ids", {'n', "no-rename"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    size_t n_threads = args::get(num_threads);
    if (n_threads) {
        omp_set_num_threads(args::get(num_threads));
    } else {
        omp_set_num_threads(1);
    }

    if (!args::get(gfa_in_file).empty() && !file_exists(args::get(gfa_in_file))) {
        std::cerr << "[gimbricate] ERROR: input sequence file " << args::get(gfa_in_file) << " does not exist" << std::endl;
        return 2;
    }
    /*
    if (args::get(gfa_out_file).empty() && args::get(paf_out_file).empty() && args::get(fasta_out_file).empty()) {
        std::cerr << "[gimbricate] ERROR: please specify at least one output, either GFA (-o), PAF (-p), and/or FASTA (-f)" << std::endl;
        return 4;
    }
    */

    double scale_align = args::get(expand_align) ? args::get(expand_align) : 1.0;

    char* filename = (char*) args::get(gfa_in_file).c_str();
    //std::cerr << "filename is " << filename << std::endl;
    gfak::GFAKluge gg;
    //double version = gg.detect_version_from_file(filename);
    //std::cerr << version << " be version" << std::endl;
    //assert(version == 1.0);
    /*
      uint64_t num_nodes = 0;
      gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
      ++num_nodes;
      });
      graph_t graph(num_nodes+1); // include delimiter
    */

    
    bool write_fasta = !args::get(fasta_out_file).empty();
    bool write_paf = !args::get(paf_out_file).empty();
    std::ofstream fasta_out, paf_out;
    if (write_paf) paf_out.open(args::get(paf_out_file));
    if (write_fasta) fasta_out.open(args::get(fasta_out_file));

    std::vector<std::string> seq_names;
    gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            seq_names.push_back(s.name);
        });
    // build the pmhf
    node_index nidx(seq_names);
    // store the seqs for random access by node name
    std::vector<std::string> seqs(seq_names.size());
    std::vector<uint64_t> ids(seq_names.size());
    seq_names.clear();
    std::cout << "H\tVN:Z:1.0" << std::endl;
    auto min_read_cov = args::get(read_cov_min);
    bool rename_seqs = !args::get(no_rename);
    uint64_t id = 1;
    gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            bool to_filter = false;
            if (min_read_cov != 0) {
                for (auto& opt : s.opt_fields) { //.push_back(o);
                    if (opt.key == "RC"
                        && std::stoul(opt.val) < min_read_cov) {
                        to_filter = true;
                    }
                }
            }
            // if we're not filtered
            if (!to_filter) {
                // save the sequence and write GFA lines here to output
                uint64_t nidx_id = nidx.get_id(s.name);
                seqs[nidx_id] = s.sequence;
                if (rename_seqs) {
                    ids[nidx_id] = id;
                    s.name = std::to_string(id);
                    ++id;
                }
                std::cout << s.to_string_1() << std::endl;
                if (write_fasta) {
                    fasta_out << ">" << s.name << std::endl
                              << s.sequence << std::endl;
                }
            }
        });
    gg.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
            auto cigar = split_cigar(e.alignment);
            uint64_t len = std::round(cigar_length(cigar) * scale_align);
            uint64_t source_nid = nidx.get_id(e.source_name);
            uint64_t sink_nid = nidx.get_id(e.sink_name);
            // flip around double inversions
            /*
            if (e.source_orientation_forward == e.sink_orientation_forward
                && !e.source_orientation_forward) {
                e.source_orientation_forward = true;
                e.sink_orientation_forward = true;
                std::swap(source_nid, sink_nid);
            }
            */
            std::string source_seq = seqs[source_nid];
            std::string sink_seq = seqs[sink_nid];
            if (rename_seqs) {
                e.source_name = std::to_string(ids[source_nid]);
                e.sink_name = std::to_string(ids[sink_nid]);
            }
            // if it's not filtered
            if (!(source_seq.empty() || sink_seq.empty())) {
                if (!e.source_orientation_forward) reverse_complement_in_place(source_seq);
                if (!e.sink_orientation_forward) reverse_complement_in_place(sink_seq);
                uint64_t query_start=0, query_end=0, target_start=0, target_end=0, num_matches=0;
                std::string basic_cigar;
                e.alignment = align_ends(source_seq, sink_seq, len,
                                         e.source_orientation_forward,
                                         e.sink_orientation_forward,
                                         query_start, query_end,
                                         target_start, target_end,
                                         num_matches, basic_cigar);
                std::cout << e.to_string_1() << std::endl;
                if (write_paf) {
                    paf_row_t paf;
                    paf.query_sequence_name = e.sink_name;
                    paf.query_sequence_length = sink_seq.length();
                    paf.query_start = query_start;
                    paf.query_end = query_end;
                    paf.query_target_same_strand = e.source_orientation_forward == e.sink_orientation_forward;
                    paf.target_sequence_name = e.source_name;
                    paf.target_sequence_length = source_seq.length();
                    paf.target_start = target_start;
                    paf.target_end = target_end;
                    paf.num_matches = num_matches;
                    paf.alignment_block_length = 0;
                    paf.mapping_quality = 100;
                    paf.cigar = (!paf.query_target_same_strand ? flip_cigar(basic_cigar) : basic_cigar);
                    paf_out << paf << std::endl;
                }
            }
        });

    return(0);
}
