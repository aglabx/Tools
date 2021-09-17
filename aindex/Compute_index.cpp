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
#include "debrujin.hpp"


#include <cstdint>
#include <cstring>
#include <limits.h>
#include <math.h>
#include <mutex>

#include "kmers.hpp"

#include "hash.hpp"

#include "read.hpp"
#include "emphf/common.hpp"
#include <cassert>

static std::mutex barrier;
emphf::stl_string_adaptor str_adapter;


int main(int argc, char** argv) {

    if (argc < 6) {
        std::cerr << "Compute LU index for reads with pf." << std::endl;
        std::cerr << "Expected arguments: " << argv[0]
        << " <dat_file> <pf_file> <output_prefix> <nthreads> <mock_flag if no dat file>" << std::endl;
        std::terminate();
    }

    PHASH_MAP hash_map;

    std::string dat_filename = argv[1];
    std::string hash_filename = argv[2]; // precomputed emphf hash
    std::string output_prefix = argv[3];
    int n_threads = atoi(argv[4]);
    int MOCK_DAT_FILE = atoi(argv[5]);


    emphf::logger() << "Loading hash..." << std::endl;
    index_hash_pp(hash_map, dat_filename, hash_filename, n_threads, MOCK_DAT_FILE);
    emphf::logger() << "\tDone." << std::endl;

    emphf::logger() << "Save raw data." << std::endl;
    emphf::logger() << "N=" << hash_map.n << std::endl;

    std::ofstream fout3(output_prefix+".kmers.bin", std::ios::out | std::ios::binary);
    emphf::logger() << "Kmer array size: " << sizeof(uint64_t) * hash_map.n <<  std::endl;
    fout3.write(reinterpret_cast<const char*> (hash_map.checker), sizeof(uint64_t) * hash_map.n);
    fout3.close();

    std::ofstream fout4(output_prefix+".tf.bin", std::ios::out | std::ios::binary);
    emphf::logger() << "TF array size: " << sizeof(unsigned int) * hash_map.n <<  std::endl;
    fout4.write(reinterpret_cast<const char*> (hash_map.tf_values), sizeof(unsigned int) * hash_map.n);
    fout4.close();

//    std::cout << "Dumping kmers" << std::endl;

//    std::ofstream fout5(output_prefix+".pf.kmers", std::ios::out);
//    for (size_t pos=0; pos < hash_map.n; ++pos) {
//        if (pos && pos % 1000000 == 0) {
//            std::cout << pos << std::endl;
//        }
//        fout5 << get_bitset_dna23(hash_map.checker[pos]) << "\n";
//    }
//    fout5.close();

    emphf::logger() << "Done." << std::endl;

    return 0;
}
