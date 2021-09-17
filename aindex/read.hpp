//
// Created by Aleksey Komissarov on 29/08/15.
// akomissarov@dell:~/Dropbox/Ariadna/Stirka$ g++ -std=c++11 -pthread -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive Compute_forks.cpp kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp correction.cpp correction.hpp read.cpp read.hpp statistics.hpp settings.hpp settings.cpp -o bin/V2_forks.exe -O3 -rdynamic



#ifndef STIRKA_READ_H
#define STIRKA_READ_H

#include <string>
#include <vector>
#include "settings.hpp"
#include "hash.hpp"
#include <math.h>
#include <bitset>
#include "dna_bitseq.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h> // memchr
#include <fcntl.h>

const int MINIMAL_PHRED = 33;
const int MAXIMAL_PHRED = 125;


//bool has_ending(std::string const &fullString, std::string const &ending) {
//    if (fullString.length() >= ending.length()) {
//        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
//    } else {
//        return false;
//    }
//}

namespace READS {

    struct Correction {

        std::string type = "unknown";
        std::string prev_type = "unknown";
        unsigned int start_position = 0;
        unsigned int end_postition = 0;
        unsigned int stop_point = 0;
        std::vector<int> fixed;
        int status = 0;

        Correction();

        Correction(std::string _type, std::string _prev_type, unsigned int _start_position) {
            type = _type;
            prev_type = _prev_type;
            start_position = _start_position;
        }
    };

    static int total_rid = 0;

    struct READ {

        int rid = 0;
        std::string head;
        std::string seq;
        std::string strand;
        std::string Q;
        int *fm;
        int *am;
        std::string cov = "";
        int status = 0;
        int position = 0;
        std::vector<int> fixed;
        std::vector<int> adapters;
        int left_status = 0;
        int right_status = 0;
        int left_position = 0;
        size_t n = 0;
        int meanq = 0;
        int solution = 0;
        int mi = 0;
        int ma = 0;
        int lmi = 0;
        int lma = 0;
        int rmi = 0;
        int rma = 0;
        bool is_string = false;

        std::vector<Correction> corrections;

        READ() {
            head = "";
            seq = "";
            strand = "";
            Q = "";

            n = 0;
            fm = nullptr;
            am = nullptr;

        }

        READ(std::string &line_head, std::string &line_seq, std::string &line_strand, std::string &line_Q) {

            head = line_head;
            seq = line_seq;
            strand = line_strand;
            Q = line_Q;

            n = (size_t) line_seq.length();
            fm = new int[n];
            am = new int[n];

            for (size_t i = 0; i < n; i++) {
                fm[i] = 0;
                am[i] = MAXIMAL_PHRED;
            }

            rid = total_rid;
            total_rid += 1;

        }

        size_t vlength() {
            return right_status - left_status;
        }

        void reverse() {
            std::string rev = seq;
            get_revcomp(seq, rev);
            seq = rev;
            size_t temp = left_status;
            left_status = seq.length() - right_status;
            right_status = seq.length() - temp;
        }


        READ(std::string &line_seq, std::string &line_Q) {

            seq = line_seq;
            Q = line_Q;

            n = (size_t) line_seq.length();

            rid = total_rid;
            total_rid += 1;

            fm = nullptr;
            am = nullptr;
        }

        READ(std::string &line_seq) {

            seq = line_seq;
            n = (size_t) line_seq.length();

            rid = total_rid;
            total_rid += 1;

            fm = nullptr;
            am = nullptr;
        }

        size_t length() {
            return n;
        }

        ~READ() {
            if (fm != nullptr) delete[] fm;
            if (am != nullptr) delete[] am;
        }

        void save_as_fastq(std::ofstream &fh) {
            fh << head << "\n";
            fh << seq << "\n";
            fh << strand << " " << status << "\n";
            fh << Q << "\n";
        }

        void save_as_read(std::ofstream &fh) {
            fh << seq << "\n";
        }

        void save_as_stats(std::ofstream &fh, int coverage, PHASH_MAP &kmers) {

            fh << "0 0 0 0 0 0 S" << status << " " << position << " " << left_status << " " << left_position << " " << meanq;
            if (adapters.size() > 0) {
                fh << " A ";
                for (size_t p = 0; p < adapters.size(); p++) {
                    fh << adapters[p] << " ";
                }
            }
            fh << "\n";
            fh << seq << "\n";
            for (size_t p = 0; p < seq.length() - Settings::K+1; p++) {
                fh << (char) am[p];
            }
            fh << "\n";
            for (size_t p = 0; p < seq.length()-Settings::K+1; p++) {
                fh << 1;
                if (p != seq.length()-Settings::K) {
                    fh << ',';
                }
            }
            fh << "\n";
            for (size_t p = 0; p < seq.length()-Settings::K+1; p++) {
                fh << fm[p];
                if (p != seq.length()-Settings::K) {
                    fh << ',';
                }
            }
            fh << "\n";
        }

        bool set_fm(PHASH_MAP &kmers) {

            bool correct_one = true;

            if (fm == nullptr) {
                fm = new int[n];
            }
            for (size_t i = 0; i < n - Settings::K + 1; i++) {
                std::string kmer = seq.substr(i, Settings::K);
                fm[i] = kmers.get_freq(kmer);
                if (fm[i] <= Settings::TRUE_ERRORS) {
                    correct_one = false;
                }

            }
            for (size_t i = n - Settings::K + 1; i < n; i++) {
                fm[i] = 0;
            }

            return correct_one;

        }

        void set_am(int coverage) {
            /*
             *   Convert Q to freqQ.
             */
            if (am == nullptr) {
                am = new int[n];
            }
            for (size_t i = 0; i < n - Settings::K + 1; i++) {
                am[i] = std::min( (int)std::round(fm[i]/coverage), 9);
            }
            for (size_t i = n - Settings::K + 1; i < n; i++) {
                am[i] = 0;
            }
        }

        void cut_end_from(size_t position) {
            if (position == 0) {
                Q = "";
                seq = "";
                n = 0;
                return;
            }
            seq = seq.substr(0, position);
            if (Q.length()) Q = Q.substr(0, position);
            n = (size_t) seq.length();
        }

        void cut_start_to(size_t position) {
            if (position+1 >= n) {
                Q = "";
                seq = "";
                n = 0;
                return;
            }
            seq = seq.substr(position+1);
            if (Q.length()) Q = Q.substr(position+1);
            n = (size_t) seq.length();
        }

        char at(size_t pos) {
            return seq[pos];
        }

        char atq(size_t pos) {
            return (int)Q[pos]-33;
        }
    };

    struct SIMPLE_READ {

        int rid = 0;
        dna_bitset* bitdna;
        size_t n = 0;

        SIMPLE_READ(std::string &line_seq) {

            bitdna = new dna_bitset(line_seq.c_str(), line_seq.length());
            n = (size_t) line_seq.length();
        }

        size_t length() {
            return n;
        }

        ~SIMPLE_READ() {
            delete bitdna;
        }

        std::string kmer(size_t pos, int k=23) {
            char* kmer = new char[k+1];
            bitdna->kmer(20, k, kmer);
            return std::string(kmer);
        }

        char at(size_t pos) {
            return bitdna->at(pos);
        }

        void save_as_read(std::ofstream &fh) {
            fh << bitdna->to_string() << "\n";
        }

        void cut_end_from(size_t position) {

            if (position == 0) {
                delete bitdna;
                bitdna = new dna_bitset("\0", 0);
                n = 0;
                return;
            }
            std::string read = std::string(bitdna->to_string());
            read = read.substr(0, position);
            n = (size_t) read.length();
            bitdna = new dna_bitset(read.c_str(), n);
        }

        void cut_start_to(size_t position) {
            if (position+1 >= n) {
                delete bitdna;
                bitdna = new dna_bitset("\0", 0);
                n = 0;
                return;
            }
            std::string read = std::string(bitdna->to_string());
            read = read.substr(position+1);
            n = (size_t) read.length();
            bitdna = new dna_bitset(read.c_str(), n);
        }

        std::string seq() {
            return std::string(bitdna->to_string());
        }

    };

    struct SPRING {

        int rid = 0;
        std::string seq;
        size_t n = 0;
        int mi = 0;
        int ma = 0;
        int lmi = 0;
        int lma = 0;
        int rmi = 0;
        int rma = 0;

        SPRING() {
            seq = "";
            n = 0;
        }

        SPRING(std::string &line_seq) {

            seq = line_seq;

            n = (size_t) line_seq.length();
            rid = total_rid;
            total_rid += 1;

        }

        size_t length() {
            return n;
        }

        ~SPRING() {

        }

    };

    struct SPRING_PAIR {

        READ *read1;
        READ *read2;


        int rid = 0;
        std::string seq;
        size_t n = 0;
        int mi = 0;
        int ma = 0;
        int lmi = 0;
        int lma = 0;
        int rmi = 0;
        int rma = 0;

        SPRING_PAIR() {
            seq = "";
            n = 0;
        }

        SPRING_PAIR(std::string &line_seq) {

            seq = line_seq;

            n = (size_t) line_seq.length();
            rid = total_rid;
            total_rid += 1;

            int pos = seq.find('~', 0);
            std::string seq1 = seq.substr(0, pos);
            std::string seq2 = seq.substr(pos+1);
            read1 = new READ(seq1);
            read2 = new READ(seq2);
        }

        SPRING_PAIR(READS::READ* read1, READS::READ* read2) {
            n = read1->length() + read2->length();
            rid = total_rid;
            total_rid += 1;
//            read1 = new READ(read1.seq);
//            read2 = new READ(read2.seq);
        }

        size_t length() {
            return n;
        }

        ~SPRING_PAIR() {

        }

        void save_as_read(std::ofstream &fh) {
            if (read1->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << read1->seq;
            }
            if (read1->n >= Settings::MINIMAL_READ_LENGTH && read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << "~";
            }
            if (read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << read2->seq;
            }
            fh << "\n";

        }

        void save_for_jellyfish(std::ofstream &fh) {

            if (read1->n < Settings::MINIMAL_READ_LENGTH && read2->n < Settings::MINIMAL_READ_LENGTH) {
                return;
            }

            fh << ">" << rid << "\n";

            if (read1->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << read1->seq;
            }
            if (read1->n >= Settings::MINIMAL_READ_LENGTH && read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << "N";
            }
            if (read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << read2->seq;
            }
            fh << "\n";

        }

        void save_as_am(std::ofstream &fh) {

            if (read1->n >= Settings::MINIMAL_READ_LENGTH) {
                for (size_t i = 0; i < read1->n; i++) {
                    fh << read1->am[i];
                }
            }
            if (read1->n >= Settings::MINIMAL_READ_LENGTH && read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << "~";
            }
            if (read2->n >= Settings::MINIMAL_READ_LENGTH) {
                for (size_t i = 0; i < read2->n; i++) {
                    fh << read2->am[i];
                }
            }
            fh << "\n";
        }

        void save_as_fm(std::ofstream &fh) {

            int tf = 0;
            if (read1->n >= Settings::MINIMAL_READ_LENGTH) {
                for (size_t i = 0; i < read1->n; i++) {
                    tf = read1->fm[i];
                    fh << tf;
                    if (i < read1->n-1) {
                        fh << " ";
                    }
                }
            }


            if (read1->n >= Settings::MINIMAL_READ_LENGTH && read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << "~";
            }
            if (read2->n >= Settings::MINIMAL_READ_LENGTH) {
                for (size_t i = 0; i < read2->n; i++) {
                    tf = read2->fm[i];
                    fh << tf;
                    if (i < read2->n-1) {
                        fh << " ";
                    }
                }
            }
            fh << "\n";

        }

    };

    struct SIMPLE_SPRING_PAIR {

        SIMPLE_READ *read1;
        SIMPLE_READ *read2;

        int rid = 0;
        std::string seq;
        size_t n = 0;

        SIMPLE_SPRING_PAIR() {
            seq = "";
            n = 0;
        }

        SIMPLE_SPRING_PAIR(std::string &line_seq) {

            seq = line_seq;
            n = (size_t) line_seq.length();
            rid = total_rid;
            total_rid += 1;

            int pos = seq.find('~', 0);
            std::string seq1 = seq.substr(0, pos);
            std::string seq2 = seq.substr(pos+1);
            read1 = new SIMPLE_READ(seq1);
            read2 = new SIMPLE_READ(seq2);

        }


        size_t length() {
            return n;
        }

        ~SIMPLE_SPRING_PAIR() {

        }

        void save_as_read(std::ofstream &fh) {
            if (read1->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << read1->seq();
            }
            if (read1->n >= Settings::MINIMAL_READ_LENGTH && read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << "~";
            }
            if (read2->n >= Settings::MINIMAL_READ_LENGTH) {
                fh << read2->seq();
            }
            fh << "\n";

        }

    };

    struct INDEXER {

        size_t n_reads = 0;
        size_t batch_size = 0;
        size_t *index = nullptr;

        const size_t BUFSIZE = 1024 * 1024;

        int * result = nullptr;

        void init_int_result(size_t size) {
            result = new int[n_reads*size]; // status1 start1 end1 status2 start2 end2
            for (size_t i=0; i<n_reads*size; ++i) {
                result[i] = 0;
            }
        }

        void set_result(size_t pos, int val) {

            result[pos] = val;
        }

        int get_val(size_t pos) {
            return result[pos];
        }

        INDEXER(const INDEXER& that) = delete;


        void save_index(std::string index_file) {
            emphf::logger() << "Saving index array..." << std::endl;
            std::ofstream fout(index_file, std::ios::out | std::ios::binary);
            emphf::logger() << "Index array size: " << sizeof(size_t) * n_reads << std::endl;
            fout.write(reinterpret_cast<const char *> (index), sizeof(size_t) * n_reads);
            fout.close();
            emphf::logger() << "Done." << std::endl;
        }

        void load_index(std::string index_file) {

            std::ifstream fout(index_file, std::ios::in | std::ios::binary);
            fout.seekg(0, std::ios::end);
            size_t length = fout.tellg();
            fout.close();

            FILE* in = std::fopen(index_file.c_str(), "rb");
            index = (size_t*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fileno(in), 0);
            if (index == nullptr) {
                std::cerr << "Failed index loading" << std::endl;
                exit(10);
            }
            fclose(in);
        }

        void load_index_raw(std::string index_file) {
            size_t f = 0;
            size_t pos = 0;
            index = new size_t[n_reads];
            std::ifstream fout(index_file, std::ios::in | std::ios::binary);
            emphf::logger() << "Loading index array..." << std::endl;
            while(fout.read(reinterpret_cast<char *>(&f), sizeof(f))) {
                index[pos] = f;
                pos += 1;
                if (pos && pos % 1000000 == 0) {
                    emphf::logger() << "\tcomputed: " << pos << std::endl;
                }
            }
            fout.close();
            emphf::logger() << "Done." << std::endl;
        }

        size_t build_fastq_index(std::string fastq_file1) {

            emphf::logger() << "\tComputing lines new way... " << std::endl;
            n_reads = 0;

//            if (has_ending(fastq_file1, ".gz") {
//
//            } else {
//                fp = open(fastq_file1.c_str(), O_RDONLY);
//            }

            int fp;
            fp = open(fastq_file1.c_str(), O_RDONLY);
            if (!fp) {
                exit(12);
            }

            std::ifstream fout(fastq_file1, std::ios::in | std::ios::binary);
            fout.seekg(0, std::ios::end);
            size_t length = fout.tellg();
            fout.close();

            char buf[BUFSIZE+1];
            int bytes_read;
            size_t readed = 0;

            while((bytes_read = read(fp, buf, BUFSIZE)) > 0){
                char* p = buf;
                while ((p = static_cast<char*>(memchr (p, '\n', (buf + bytes_read) - p)))) {
                    ++p;
                    ++n_reads;

                    if (n_reads && n_reads % 4000000 == 0) {
                        emphf::logger() << "\tcomputed: " << n_reads << " or " << 100*BUFSIZE*readed/length << "%" << std::endl;
                    }
                }

                readed += 1;
            }
            close(fp);
            n_reads /= 4;

            emphf::logger() << n_reads << " reads" << std::endl;

            index = new size_t[n_reads];
            fp = open(fastq_file1.c_str(), O_RDONLY);

            size_t read_n = 0;
            size_t prev_pos = 0;
            readed = 0;


            while((bytes_read = read(fp, buf, BUFSIZE)) > 0){
                char* p = buf;
                while ((p = static_cast<char*>(memchr (p, '\n', (buf + bytes_read) - p)))) {
                    ++p;
                    if (read_n % 4 == 0) {
                        index[read_n/4] = prev_pos;
                    }
                    if (read_n && read_n % 4000000 == 0) {
                        emphf::logger() << "\tcomputed: " << read_n << " from " << n_reads * 4 << " or " <<
                        100 * (read_n) / (n_reads * 4) << "%" << std::endl;
                    }
                    ++read_n;
                    prev_pos = BUFSIZE*readed + p - buf;
                }

                readed += 1;
            }
            close(fp);
            return n_reads;
        }

        size_t build_reads_index(std::string reads_file) {

            emphf::logger() << "\tComputing lines new way... " << std::endl;
            n_reads = 0;
            int fp = open(reads_file.c_str(), O_RDONLY);

            if (!fp) {
                exit(12);
            }

            std::ifstream fout(reads_file, std::ios::in | std::ios::binary);
            fout.seekg(0, std::ios::end);
            size_t length = fout.tellg();
            fout.close();

            char buf[BUFSIZE+1];
            int bytes_read;
            size_t readed = 0;

            while((bytes_read = read(fp, buf, BUFSIZE)) > 0){
                char* p = buf;
                while ((p = static_cast<char*>(memchr (p, '\n', (buf + bytes_read) - p)))) {
                    ++p;
                    ++n_reads;

                    if (n_reads && n_reads % 10000000 == 0) {
                        emphf::logger() << "\tcomputed: " << n_reads << " or " << 100*BUFSIZE*readed/length << "%" << std::endl;
                    }
                }
                readed += 1;
            }
            close(fp);
            emphf::logger() << n_reads << " reads" << std::endl;

            index = new size_t[n_reads];
            fp = open(reads_file.c_str(), O_RDONLY);

            size_t read_n = 0;
            size_t prev_pos = 0;
            readed = 0;

            while((bytes_read = read(fp, buf, BUFSIZE)) > 0){
                char* p = buf;
                while ((p = static_cast<char*>(memchr (p, '\n', (buf + bytes_read) - p)))) {
                    ++p;
                    index[read_n] = prev_pos;
                    if (read_n && read_n % 10000000 == 0) {
                        emphf::logger() << "\tcomputed: " << read_n << " from " << n_reads << " or " <<
                        100 * (read_n) / (n_reads) << "%" << std::endl;
                    }
                    ++read_n;
                    prev_pos = BUFSIZE*readed + p - buf;
                }

                readed += 1;
            }
            close(fp);
            return n_reads;
        }

        READS::READ * get_fastq_reads(std::ifstream &r_file, size_t rid) {
            std::string line_head = "";
            std::string line_seq = "";
            std::string line_strand = "";
            std::string line_Q = "";

            r_file.seekg(index[rid]);

            std::getline(r_file, line_head);
            std::getline(r_file, line_seq);
            std::getline(r_file, line_strand);
            std::getline(r_file, line_Q);

            READS::READ * read = new READS::READ(line_head, line_seq, line_strand, line_Q);

            return read;

        }

        READS::SPRING_PAIR * get_spring_read(std::ifstream &r_file, size_t rid) {
            std::string line_seq = "";
            r_file.seekg(index[rid]);
            std::getline(r_file, line_seq);
            READS::SPRING_PAIR *read = new SPRING_PAIR(line_seq);
            return read;
        }


        void worker_count_lines(std::string read_file, std::vector<size_t> *r, size_t start, size_t end, int i) {
//            barrier.lock();
//            emphf::logger() << "Started from " << start << " to " << end << std::endl;
//            barrier.unlock();
        }

//        size_t build_reads_index_pp(std::string reads_file, int num_threads) {
//
//            /// compute file size
//
//            emphf::logger() << "\tComputing lines new way... " << std::endl;
//            n_reads = 0;
//            int fp = open(reads_file.c_str(), O_RDONLY);
//            if (!fp) {
//                exit(12);
//            }
//            std::ifstream fout(reads_file, std::ios::in | std::ios::binary);
//            fout.seekg(0, std::ios::end);
//            size_t length = fout.tellg();
//            fout.close();
//
//            /// split it in batches
//
//            size_t batch_size;
//            batch_size = (length / num_threads) + 1;
//            std::vector<std::thread> t;
//            void* results = operator new(num_threads);
//
//            for (size_t i = 0; i < num_threads; ++i) {
//                std::vector<size_t> t;
//                results[i] = &t;
//            }
//
//            for (size_t i = 0; i < num_threads; ++i) {
//                size_t start = i * batch_size;
//                size_t end = (i + 1) * batch_size;
//
//                if (end > n_reads) {
//                    end = n_reads;
//                }
//
//                t.push_back(std::thread(worker_count_lines,
//                                        std::ref(read_file),
//                                        std::ref(results[i]),
//                                        start,
//                                        end,
//                                        i
//                ));
//            }
//            for (size_t i = 0; i < num_threads; ++i) {
//                t[i].join();
//            }
//            emphf::logger() << "\tDone." << std::endl;
//
//            operator delete(results);
//
//            return 0;
//
//
//            char buf[BUFSIZE+1];
//            int bytes_read;
//            size_t readed = 0;
//
//            while((bytes_read = read(fp, buf, BUFSIZE)) > 0){
//                char* p = buf;
//                while ((p = static_cast<char*>(memchr (p, '\n', (buf + bytes_read) - p)))) {
//                    ++p;
//                    ++n_reads;
//
//                    if (n_reads && n_reads % 1000000 == 0) {
//                        emphf::logger() << "\tcomputed: " << n_reads << " or " << 100*BUFSIZE*readed/length << "%" << std::endl;
//                    }
//                }
//                readed += 1;
//            }
//            close(fp);
//            emphf::logger() << n_reads << " reads" << std::endl;
//
//            index = new size_t[n_reads];
//            fp = open(reads_file.c_str(), O_RDONLY);
//
//            size_t read_n = 0;
//            size_t prev_pos = 0;
//            readed = 0;
//
//            while((bytes_read = read(fp, buf, BUFSIZE)) > 0){
//                char* p = buf;
//                while ((p = static_cast<char*>(memchr (p, '\n', (buf + bytes_read) - p)))) {
//                    ++p;
//                    index[read_n] = prev_pos;
//                    if (read_n && read_n % 1000000 == 0) {
//                        emphf::logger() << "\tcomputed: " << read_n << " from " << n_reads * 4 << " or " <<
//                        100 * (read_n) / (n_reads * 4) << "%" << std::endl;
//                    }
//                    ++read_n;
//                    prev_pos = BUFSIZE*readed + p - buf;
//                }
//
//                readed += 1;
//            }
//            close(fp);
//            return n_reads;
//        }

        INDEXER() {

        }

        ~INDEXER() {
            if (index != nullptr) {
                delete [] index;
            }
            if (result != nullptr) {
                delete [] result;
            }
        }
    };

    struct READ_PAIR {

        READ *read1;
        READ *read2;
        bool spring;

        READ_PAIR(READ *_read1, READ *_read2) {
            read1 = _read1;
            read2 = _read2;
            spring = false;
        }

        ~READ_PAIR() {
            delete read1;
            delete read2;
        }

        int save_as_fastq(std::ofstream &fh1, std::ofstream &fh2) {
            if (read1->n >= Settings::MINIMAL_READ_LENGTH && read2->n >= Settings::MINIMAL_READ_LENGTH) {
                read1->save_as_fastq(fh1);
                read2->save_as_fastq(fh2);
                return 1;
            }
            return 0;
        }
    };



    void read_sreads(std::string file_name, std::vector<READ *> &reads);

    void read_springs(std::string file_name, std::vector<SPRING *> &reads);

    void read_reads(std::string file_name, std::vector<READ *> &reads, size_t n);

    void read_pair_reads(std::string file_name1, std::string file_name2, std::vector<READ_PAIR *> &read_pairs);

    void print_read(READ &read, int coverage);

    int get_fm_mode(READ &read);

    void read_spring_pairs(std::string file_name, std::vector<SPRING_PAIR *> &reads, size_t n_reads);

    void read_simple_spring_pairs(std::string file_name, std::vector<SPRING_PAIR *> &reads, size_t nreads);
}

#endif //STIRKA_READ_H
