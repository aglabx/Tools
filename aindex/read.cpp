//
// Created by Aleksey Komissarov on 29/08/15.
//

#include "read.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <mutex>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include "emphf/common.hpp"

namespace READS {

    std::mutex barrier;

    void read_reads(std::string file_name, std::vector<READ *> &reads, size_t n) {
        // Read fastq file into memory.

        barrier.lock();
        emphf::logger() << "Starting load reads..." << std::endl;
        barrier.unlock();

        std::ifstream infile(file_name);
        infile.seekg(0, std::ios::end);
        size_t length = infile.tellg();
        char *contents = new char[length + 1];

        infile.seekg(0, std::ios::beg);
        infile.read(contents, length);
        infile.close();

        contents[length] = 0;
        std::stringstream ss_contents;
        ss_contents << contents;

        std::string line_head = "";
        std::string line_seq = "";
        std::string line_strand = "";
        std::string line_Q = "";

        size_t loaded = 1;
        while (std::getline(ss_contents, line_head)) {
            std::getline(ss_contents, line_seq);
            std::getline(ss_contents, line_strand);
            std::getline(ss_contents, line_Q);
            READ *read = new READ(line_seq, line_Q);
            reads.push_back(read);
            if (loaded % 100000 == 0){
                barrier.lock();
                emphf::logger() << "Loaded " << loaded << "/" << n << std::endl;
                barrier.unlock();
            }
            loaded += 1;

        }
        barrier.lock();
        emphf::logger() << "Loaded reads: " << reads.size() << std::endl;
        barrier.unlock();

        delete[] contents;
    }

    void read_sreads(std::string file_name, std::vector<READ *> &reads) {
        // Read springs file into memory.

        barrier.lock();
        emphf::logger() << "Starting load sreads..." << std::endl;
        barrier.unlock();

        std::ifstream infile(file_name);
        infile.seekg(0, std::ios::end);
        size_t length = infile.tellg();
        char *contents = new char[length + 1];

        infile.seekg(0, std::ios::beg);
        infile.read(contents, length);
        infile.close();

        contents[length] = 0;
        std::stringstream ss_contents;
        ss_contents << contents;

        std::string line_head = "";
        std::string line_seq = "";
        std::string line_strand = "";
        std::string line_Q = "";

        while (std::getline(ss_contents, line_seq)) {
            READ *read = new READ(line_head, line_seq, line_strand, line_Q);
            reads.push_back(read);
        }

        barrier.lock();
        emphf::logger() << "Read: " << reads.size() << std::endl;
        barrier.unlock();

        delete[] contents;

    }

    void read_springs(std::string file_name, std::vector<SPRING *> &reads) {
        // Read springs file into memory.

        barrier.lock();
        emphf::logger() << "Starting load springs..." << std::endl;
        barrier.unlock();

        std::ifstream infile(file_name);
        infile.seekg(0, std::ios::end);
        size_t length = infile.tellg();
        char *contents = new char[length + 1];

        infile.seekg(0, std::ios::beg);
        infile.read(contents, length);
        infile.close();

        contents[length] = 0;
        std::stringstream ss_contents;
        ss_contents << contents;

        std::string line_seq = "";

        size_t loaded = 1;
        while (std::getline(ss_contents, line_seq)) {
            SPRING *read = new SPRING(line_seq);
            loaded += 1;
            if (loaded % 1000000 == 0 ) {
                barrier.lock();
                emphf::logger() << "Loaded " << loaded << " reads." << std::endl;
                barrier.unlock();
            }
            reads.push_back(read);
        }
        barrier.lock();
        emphf::logger() << "Read: " << reads.size() << std::endl;
        barrier.unlock();

        delete[] contents;
    }




    char* sgets(char *s, int n, char **strp){
        if (**strp == '\0') return NULL;
        int i;
        for(i=0; i<n-1; ++i, ++(*strp)){
            s[i] = **strp;
            if(**strp == '\0')
                break;
            if(**strp == '\n') {
                s[i]='\0';
                ++(*strp);
                break;
            }
        }
        if(i==n-1)
            s[i] = '\0';
        return s;
    }

    void read_spring_pairs(std::string file_name, std::vector<SPRING_PAIR *> &reads, size_t n_reads) {
        // Read springs file into memory.

        barrier.lock();
        emphf::logger() << "Starting load springs..." << std::endl;
        barrier.unlock();

        std::ifstream infile(file_name, std::ios_base::binary);
        infile.seekg(0, std::ios::end);
        size_t length = infile.tellg();

        emphf::logger() << "Trying allocate (GB): " << (sizeof(char) * (length+1)) /(1024ll*1024*1024*8)  << std::endl;
        char *contents = new char[length + 1], *contents0 = contents;
        if (contents == nullptr) {
            emphf::logger() << "Failed to allocate memory for file content (char): " << length+1 << std::endl;
            exit(10);
        }
        emphf::logger() << "Done." << std::endl;

        infile.seekg(0, std::ios::beg);
        infile.read(contents, length);
        infile.close();

        contents[length] = '\n';
        emphf::logger() << "Done read." << std::endl;

        std::string line_seq = "";

        char buff[1000];
        char **p = &contents;

        size_t loaded = 0;
        while(NULL != sgets(buff, sizeof(buff), p)) {
            line_seq = std::string(buff);
            SPRING_PAIR *read = new SPRING_PAIR(line_seq);
            loaded += 1;
            if (loaded % 1000000 == 0 ) {
                barrier.lock();
                emphf::logger() << "Loaded " << loaded << " reads (" << n_reads << ")" << std::endl;
                barrier.unlock();
            }
            reads.push_back(read);
            if (loaded == n_reads) {
                break;
            }
        }
        barrier.lock();
        emphf::logger() << "Read: " << reads.size() << std::endl;
        barrier.unlock();

        delete[] contents0;
    }

    void read_simple_spring_pairs(std::string file_name, std::vector<SPRING_PAIR *> &reads, size_t n_reads) {


        barrier.lock();
        emphf::logger() << "Starting load springs..." << std::endl;
        barrier.unlock();

        std::ifstream infile(file_name, std::ios_base::binary);
        infile.seekg(0, std::ios::end);
        size_t length = infile.tellg();

        emphf::logger() << "Trying allocate (GB): " << (sizeof(char) * (length+1ll)) /(1024ll*1024*1024*8)  << std::endl;

        char *contents = new char[length + 1], *contents0 = contents;
        if (contents == nullptr) {
            emphf::logger() << "Failed to allocate memory for file content (char): " << length+1 << std::endl;
            exit(10);
        }
        emphf::logger() << "Done." << std::endl;

        infile.seekg(0, std::ios::beg);
        infile.read(contents, length);
        infile.close();

        contents[length] = '\n';
        emphf::logger() << "Done read." << std::endl;

        std::string line_seq = "";

        char buff[1000];
        char **p = &contents;

        size_t loaded = 0;
        while(NULL != sgets(buff, sizeof(buff), p)) {
            line_seq = std::string(buff);
            SPRING_PAIR *read = new SPRING_PAIR(line_seq);
            loaded += 1;
            if (loaded && loaded % 1000000 == 0 ) {
                barrier.lock();
                emphf::logger() << "Loaded " << loaded << " reads (" << n_reads << ")" << std::endl;
                barrier.unlock();
            }
            reads.push_back(read);
            if (loaded == n_reads) {
                break;
            }
        }
        barrier.lock();
        emphf::logger() << "Read: " << reads.size() << std::endl;
        barrier.unlock();
        delete[] contents0;

    }

    void read_pair_reads(std::string file_name1, std::string file_name2, std::vector<READ_PAIR *> &read_pairs) {

        barrier.lock();
        emphf::logger() << "Starting load read pairs..." << std::endl;
        barrier.unlock();




        std::string line_head1 = "";
        std::string line_seq1 = "";
        std::string line_strand1 = "";
        std::string line_Q1 = "";

        std::string line_head2 = "";
        std::string line_seq2 = "";
        std::string line_strand2 = "";
        std::string line_Q2 = "";

        std::ifstream fh1(file_name1);
        std::ifstream fh2(file_name2);
        size_t n_reads = 0;


        while (std::getline(fh1, line_head1)) {
            std::getline(fh1, line_seq1);
            std::getline(fh1, line_strand1);
            std::getline(fh1, line_Q1);

            std::getline(fh2, line_head2);
            std::getline(fh2, line_seq2);
            std::getline(fh2, line_strand2);
            std::getline(fh2, line_Q2);

            READ *read1 = new READ(line_head1, line_seq1, line_strand1, line_Q1);
            READ *read2 = new READ(line_head2, line_seq2, line_strand2, line_Q2);

            READ_PAIR *read_pair = new READ_PAIR(read1, read2);
            read_pairs.push_back(read_pair);

            n_reads += 1;

            if (n_reads && n_reads % 1000000 == 0) {
                emphf::logger() << "\tloaded " << n_reads << std::endl;
            }
        }

        fh1.close();
        fh2.close();
    }

    void print_read(READ &read, int coverage) {
        std::cout << "ID: " << read.rid << std::endl;
        std::cout << "Header: " << read.head << std::endl;
        std::cout << "Sequence: " << read.seq << std::endl;
        std::cout << "Length: " << read.seq.length() << std::endl;
        std::cout << "Strand: " << read.strand << std::endl;
        std::cout << "Quality : " << read.Q << std::endl;
        std::cout << "Status: " << read.status << std::endl;
        std::cout << "Position: " << read.position << std::endl;
        std::cout << "Left Status: " << read.left_status << std::endl;
        std::cout << "Left Position: " << read.left_position << std::endl;
        int i;
        char c;
        std::cout << "FM: ";
        for (size_t p = 0; p < read.seq.length() - Settings::K + 1; p++) {
            i = read.fm[p];
            std::cout << i << " ";
        }
        std::cout << std::endl;
        read.set_am(coverage);
        std::cout << "AM: ";
        for (size_t p = 0; p < read.seq.length() - Settings::K + 1; p++) {
            c = read.am[p];
            std::cout << c << " ";
        }
        std::cout << std::endl;

        std::cout << "Adapters (" << read.adapters.size() << "): ";
        for (auto it = read.adapters.begin(); it != read.adapters.end(); it++) {
            std::cout << *it << " ";
        }
        std::cout << std::endl;

        char cmeanq = (char) read.meanq;
        std::cout << "MeanQ : " << cmeanq << " (" << read.meanq << ")" << std::endl;
        std::cout << "Fixed: ";
        for (size_t p = 0; p < read.fixed.size(); p++) {
            i = read.fixed[p];
            std::cout << i << " ";
        }
        std::cout << std::endl;

        std::cout << "Corrections (" << read.corrections.size() << "): " << std::endl;
        for (auto it = read.corrections.begin(); it != read.corrections.end(); it++) {
            std::cout << " " << it->prev_type << " ";
            std::cout << " " << it->type << " ";
            std::cout << " " << it->start_position << " ";
            std::cout << " " << it->stop_point << " ";
            std::cout << " " << it->end_postition << " ";
            std::cout << " " << it->status << std::endl;
        }
    }

    void print_read_adapter_view(READ &read, std::string &message) {
        std::cout << "From print_read_adapter_view: " << message << std::endl;
        std::cout << read.seq << std::endl;
        std::cout << read.Q << std::endl;
        for (auto it = read.adapters.begin(); it != read.adapters.end(); it++) {
            std::cout << *it << " ";
        }
        std::cout << std::endl;
        for (size_t j = 0; j < read.seq.length(); j++) {
            std::cout << (char) read.am[j];
        }
        std::cout << std::endl;
        for (size_t j = 0; j < read.seq.length(); j++) {
            std::cout << read.fm[j] << " ";
        }
        std::cout << std::endl;
    }

    void print_read_fragment_view(READ &read, std::string &message, size_t start, size_t end) {
        std::cout << "From print_read_adapter_view: " << message << std::endl;
        std::cout << read.seq.substr(start, end - start) << std::endl;
        std::cout << read.Q.substr(start, end - start) << std::endl;
        std::cout << std::endl;
        for (size_t j = start; j < end; j++) {
            std::cout << (char) read.am[j];
        }
        std::cout << std::endl;
        for (size_t j = start; j < end; j++) {
            std::cout << read.fm[j] << " ";
        }
        std::cout << std::endl;
    }

    int get_fm_mode(READ &read) {
        int size = read.seq.length() - Settings::K + 1;
        int *_reps = new int[size];
        for (int i = 0; i < size; ++i) {
            _reps[i] = 0;
            int j = 0;
            while ((j < i) && (read.fm[i] != read.fm[j])) {
                if (read.fm[i] != read.fm[j]) {
                    ++j;
                }
            }
            ++(_reps[j]);
        }
        int _max = 0;
        for (int i = 1; i < size; ++i) {
            if (_reps[i] > _reps[_max]) {
                _max = i;
            }
        }
        delete[] _reps;
        return read.fm[_max];
    }
}