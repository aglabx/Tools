//
// Created by Aleksey Komissarov on 22/01/16.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <thread>
#include "emphf/common.hpp"
#include <cstdint>
#include <cstring>
#include <limits.h>
#include <math.h>
#include <mutex>
#include <algorithm>
#include <assert.h>
#include "kmers.hpp"
#include <atomic>
#include "hash.hpp"
#include <cassert>
#include "read.hpp"
#include "emphf/common.hpp"
#include <cassert>

int main(int argc, char** argv) {

    if (argc < 8) {
        std::cerr << "Compute AIndex index for genome with pf." << std::endl;
        std::cerr << "Expected arguments: " << argv[0]
        << " <reads_file> <pf_file> <index_prefix> <output_prefix> <p> <k> <tf_file>" << std::endl;
        std::terminate();
    }

    std::vector<READS::READ *> reads;
    PHASH_MAP hash_map;

    std::string read_file = argv[1];
    std::string hash_filename = argv[2]; // precomputed emphf hash
    std::string index_prefix = argv[3];
    std::string output_prefix = argv[4];
    static const uint num_threads = atoi(argv[5]);
    Settings::K = atoi(argv[6]);
    std::string line;

    std::string tf_file = argv[7];

    emphf::logger() << "Loading hash..." << std::endl;

    if (Settings::K == 13) {
        load_hash_full_tf(hash_map, index_prefix, tf_file, hash_filename);
    } else {
        load_hash(hash_map, index_prefix, tf_file, hash_filename);
    }

    emphf::logger() << "\tDone. Kmers: " << hash_map.n << std::endl;

    emphf::logger() << "Load and reads and build docid index..." << std::endl;
    emphf::logger() << "Opening read_file: " << read_file << std::endl;

    std::vector<size_t> start_positions;
    std::map<size_t, uint32_t> start2rid;

    std::ifstream infile(read_file);
    if (!infile) {
        emphf::logger() << "Failed open read_file: " << read_file << std::endl;
        exit(10);
    }

    infile.seekg(0, std::ios::end);
    size_t length = infile.tellg();
    char *contents = new char[length + 1];
    if (contents == nullptr) {
        emphf::logger() << "Failed to allocate for reads: " << length + 1 << std::endl;
        exit(10);
    }

    infile.seekg(0, std::ios::beg);
    infile.read(contents, length);
    infile.close();
    contents[length] = 0;
    std::string line_seq = "";
    uint32_t rid = 0;
    size_t pos = 0;
    start_positions.push_back(pos);
    start2rid[pos] = rid;

    std::cout << "Init aindex..." << std::endl;
    AIndexCompressed aindex(hash_map);
    std::cout << "Done." << std::endl;

    size_t current_position = 0;
    size_t nreads = 0;
    for (size_t i = 0; i < length; i++) {
        if (contents[i] == '\n') {
            start_positions.push_back(i + 1);
            rid += 1;
            start2rid[i + 1] = rid;
            nreads += 1;
            if (rid % 1000000 == 0) {
                emphf::logger() << "Loaded read: " << rid << std::endl;
            }
            current_position = i + 1;
        }
    }

    emphf::logger() << "\tLoaded: " << nreads << " nreads and " << current_position << " symbols" << std::endl;
    start_positions.push_back(length);
    emphf::logger() << "\tDone." << std::endl;

    aindex.fill_index_from_reads(contents, length, num_threads, hash_map);

    aindex.save(output_prefix, start_positions, hash_map);

    emphf::logger() << "\tDone." << std::endl;
    delete[] contents;


    return 0;
}
