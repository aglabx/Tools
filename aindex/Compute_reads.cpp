//
// Created by Aleksey Komissarov on 22/01/16.
//

//
// Created by Aleksey Komissarov on 28/12/15.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <thread>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits.h>
#include <math.h>
#include <mutex>
#include "emphf/common.hpp"
#include "read.hpp"

static std::mutex barrier;
emphf::stl_string_adaptor str_adapter;

int main(int argc, char** argv) {

    if (argc < 5) {
        std::cerr << "Convert fasta or fastq reads to simple reads." << std::endl;
        std::cerr << "Expected arguments: " << argv[0]
        << " <fastq_file1|fasta_file1> <fastq_file2|fasta_file2> <fastq|fasta> <output_file>" << std::endl;
        std::terminate();
    }

    std::string file_name1 = argv[1];
    std::string file_name2 = argv[2];
    std::string read_type = argv[3];
    std::string output_file = argv[4];


    emphf::logger() << "Starting..." << std::endl;

    emphf::logger() << "Converting reads..." << std::endl;
    size_t n_reads = 0;

    std::string line1;
    std::string line2;


    std::ifstream fin1(file_name1, std::ios::in);
    std::ifstream fin2(file_name2, std::ios::in);
    std::ofstream fout(output_file, std::ios::out);


    if (read_type == "fastq") {
        while (std::getline(fin1, line1)) {
            std::getline(fin1, line1);
            std::getline(fin2, line2);
            std::getline(fin2, line2);
            std::string rline2 = line2;
            get_revcomp(line2, rline2);

            fout << line1;
            fout << "~";
            fout << rline2;
            fout << "\n";

            std::getline(fin1, line1);
            std::getline(fin1, line1);
            std::getline(fin2, line2);
            std::getline(fin2, line2);
            n_reads += 1;

            if (n_reads % 1000000 == 0) {
                emphf::logger() << "Completed: " << n_reads << std::endl;
            }
        }
    } else {
        emphf::logger() << "Unknown format." << std::endl;
        exit(2);
    }

    fout.close();
    fin1.close();
    fin2.close();

    return 0;
}
