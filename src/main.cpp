#include <iostream>
#include <cstdint>
#include <cstdio>
#include <string>
#include "exists.hpp"
#include "threads.hpp"
#include "args.hxx"
#include "gssw.h"

using namespace gimbricate;

int main(int argc, char** argv) {
    args::ArgumentParser parser("gimbricate: recompute GFA link overlaps");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "use this GFA FILE as input", {'g', "gfa"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
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

    return(0);
}
