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
//#include <limits>

using namespace gimbricate;

int main(int argc, char** argv) {
    args::ArgumentParser parser("gimbricate: recompute GFA link overlaps");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "use this GFA FILE as input", {'g', "gfa"});
    args::ValueFlag<uint64_t> read_cov_min(parser, "N", "require this many supporting reads in the RC tag to keep a node", {'c', "read-coverage"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
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

    if (!args::get(gfa_in).empty() && !file_exists(args::get(gfa_in))) {
        std::cerr << "[gimbricate] ERROR: input sequence file " << args::get(gfa_in) << " does not exist" << std::endl;
        return 2;
    }

    char* filename = (char*) args::get(gfa_in).c_str();
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
            }
        });
    gg.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
            //if (e.source_name.empty()) return;
            //handlegraph::handle_t a = graph->get_handle(stol(e.source_name), !e.source_orientation_forward);
            //handlegraph::handle_t b = graph->get_handle(stol(e.sink_name), !e.sink_orientation_forward);
            //graph->create_edge(a, b);
            //std::cerr << "alignment! " << e.alignment << std::endl;
            auto cigar = split_cigar(e.alignment);
            uint64_t len = cigar_length(cigar);
            uint64_t source_nid = nidx.get_id(e.source_name);
            uint64_t sink_nid = nidx.get_id(e.sink_name);
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
                e.alignment = align_ends(source_seq, sink_seq, len);
                std::cout << e.to_string_1() << std::endl;
            }
        });

    return(0);
}
