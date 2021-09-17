#include <iostream>
#include <fstream>
#include <string>
#include "hash.hpp"
#include <vector>
#include "emphf/common.hpp"
#include <algorithm>
#include <sys/mman.h>


emphf::stl_string_adaptor str_adapter;

struct Hit {
    size_t rid;
    size_t start;
    std::string read;
    size_t pos;
    int ori;
    bool rev;
};

class AindexWrapper{

    size_t *positions = nullptr;
    size_t *indices = nullptr;

    size_t n = 0;
    uint32_t max_tf = 0;
    size_t indices_length = 0;

public:

    PHASH_MAP *hash_map;
    size_t n_reads = 0;
    std::map<size_t, uint32_t> start2rid;
    std::vector<size_t> start_positions;

    size_t *start_postitions_raw = nullptr;

    size_t reads_length = 0;
    char *reads = nullptr;

    AindexWrapper() {

    }

    ~AindexWrapper() {
        emphf::logger() << "NOTE: Calling aindex deconstructor..." << std::endl;
        if (positions != nullptr) munmap(positions, n*sizeof(size_t));
        if (indices != nullptr) munmap(indices, indices_length);
        if (reads != nullptr) munmap(reads, reads_length);

        if (start_postitions_raw != nullptr) delete[] start_postitions_raw;

        delete hash_map;

        reads = nullptr;
        indices = nullptr;
        positions = nullptr;

    }

    void load(std::string index_prefix){

        hash_map = new PHASH_MAP();
        // Load perfect hash into hash_map into memory
        emphf::logger() << "Reading index and hash..." << std::endl;
        std::string tf_file = index_prefix + ".tf.bin";
        std::string hash_filename = index_prefix + ".pf";
        emphf::logger() << "...files: " << index_prefix << std::endl;
        emphf::logger() << "...files: " << tf_file << std::endl;
        emphf::logger() << "...files: " << hash_filename << std::endl;
        load_hash(*hash_map, index_prefix, tf_file, hash_filename);
        emphf::logger() << "\tDone" << std::endl;
    }

    void load_hash_file(std::string hash_filename) {
        std::cout << "Loading only hash..." << std::endl;
        load_only_hash(*hash_map, hash_filename);
    }

    void load_reads(std::string reads_file) {
        // Memory map reads
        emphf::logger() << "Memory mapping reads file..." << std::endl;
        std::ifstream fout(reads_file, std::ios::in | std::ios::binary);
        fout.seekg(0, std::ios::end);
        size_t length = fout.tellg();
        fout.close();

        FILE* in = std::fopen(reads_file.c_str(), "rb");
        reads = (char*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in), 0);
        if (reads == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in);

        n_reads = 0;
        reads_length = length;

        for (size_t i=0; i < length; ++i) {
            if (reads[i] == '\n') n_reads += 1;
        }
        emphf::logger() << "\tloaded reads: " << n_reads << std::endl;

        emphf::logger() << "\tbuilding start pos index over reads: " << std::endl;
        start_postitions_raw = new size_t[n_reads+1];

        size_t rid = 0;
        for (size_t i=0; i < length; ++i) {
            if (reads[i] == '\n') {
                start_postitions_raw[rid] = i;
//                start2rid[i] = rid;
                rid += 1;
            }
        }
        emphf::logger() << "\tDone" << std::endl;

    }

    void load_index(std::string aindex_prefix, uint32_t _max_tf) {
        // Load aindex.

        // std::cout << "Inside load index: " << _max_tf << std::endl;

        n = hash_map->n;
        max_tf = _max_tf;

        std::string pos_file = aindex_prefix + ".pos.bin";
        std::string index_file = aindex_prefix + ".index.bin";
        std::string indices_file = aindex_prefix + ".indices.bin";

        // std::cout << "END" << std::endl;

//        emphf::logger() << "Reading aindex.pos.bin array..." << std::endl;
//
//        size_t f = 0;
//        size_t rid = 0;
//        std::ifstream fout2(pos_file, std::ios::in | std::ios::binary);
//        while(fout2.read(reinterpret_cast<char *>(&f), sizeof(f))) {
//            start_positions.push_back(f);
//            start2rid[f] = rid;
//            rid += 1;
//
//        }
//
//        start_positions.pop_back();
//        start_positions.pop_back();
//        fout2.close();
//        emphf::logger() << "\tDone" << std::endl;

        emphf::logger() << "Reading aindex.indices.bin array..." << std::endl;

        size_t pos = 0;
        std::ifstream fout5(indices_file, std::ios::in | std::ios::binary);
        fout5.seekg(0, std::ios::end);
        size_t length = fout5.tellg();
        fout5.close();

        FILE* in1 = std::fopen(indices_file.c_str(), "rb");
        indices = (size_t*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in1), 0);
        if (indices == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in1);
        indices_length = length;
        emphf::logger() << "\tDone" << std::endl;

        emphf::logger() << "Reading aindex.index.bin array..." << std::endl;

        pos = 0;
        std::ifstream fout6(index_file, std::ios::in | std::ios::binary);
        fout6.seekg(0, std::ios::end);
        length = fout6.tellg();
        fout6.close();

        FILE* in = std::fopen(index_file.c_str(), "rb");
        positions = (size_t*)mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fileno(in), 0);
        if (positions == nullptr) {
            std::cerr << "Failed position loading" << std::endl;
            exit(10);
        }
        fclose(in);
        emphf::logger() << "\tDone" << std::endl;

    }



    std::string get_read_by_rid(uint32_t rid, int ori) {

        size_t start = start_positions[rid];
        size_t end = start;
        size_t spring_pos = 0;

        while (true) {
            if (reads[end] == '\n') {
                char rkmer[end-spring_pos];
                std::memcpy(rkmer, &reads[spring_pos+1], end-spring_pos-1);
                rkmer[end-spring_pos-1] = '\0';
                return std::string(rkmer);
            } else if (reads[end] == '~') {
                if (ori == 0) {
                    char lkmer[end-start+1];
                    std::memcpy(lkmer, &reads[start], end-start);
                    lkmer[end-start] = '\0';
                    return std::string(lkmer);
                } else {
                    spring_pos = end;
                }
            }
            end += 1;
        }
    }

    size_t get_n() {
        return hash_map->n;
    }

    size_t get_rid(size_t pos) {
        // Get rid by position.
        while (true){
            if (reads[pos] == '\n') {
                return pos+1;
            }
            if (pos == 0) {
                return pos;
            }
            pos -= 1;
        }
    }

    size_t get_kid_by_kmer(std::string _kmer) {
        uint64_t kmer = get_dna23_bitset(_kmer);
        return hash_map->get_pfid_by_umer_safe(kmer);

    }

    

    void get_positions(size_t* r, std::string kmer) {
        // Get read positions and save them to given r
        auto h1 = hash_map->get_pfid(kmer);
        size_t j = 0;
        for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {
            if (j == max_tf - 1) {
                break;
            }
            r[j] = positions[i];
            j += 1;
        }
        r[j] = 0;
    }

    size_t get(char* ckmer) {
        // Return tf for given char * kmer
        std::string kmer = std::string(ckmer);
        return get(kmer);
    }

    size_t get(std::string kmer) {
        // Return tf for given kmer
        uint64_t ukmer = get_dna23_bitset(kmer);
        auto h1 = hash_map->hasher.lookup(kmer, str_adapter);
        if (h1 >= hash_map->n || hash_map->checker[h1] != ukmer) {
            std::string rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
            uint64_t urev_kmer = reverseDNA(ukmer);
            get_bitset_dna23(urev_kmer, rev_kmer);
            auto h2 = hash_map->hasher.lookup(rev_kmer, str_adapter);
            if (h2 >= hash_map->n || hash_map->checker[h2] != urev_kmer) {
                return 0;
            } else {
                return hash_map->tf_values[h2];
            }
        } else {
            return hash_map->tf_values[h1];
        }
        return 0;
    }


    void get_kmer_by_kid(size_t r, char* kmer) {
            // if (r >= hash_map->n) {
            //     return;
            // }
            uint64_t ukmer = hash_map->checker[r];
            get_bitset_dna23_c(ukmer, kmer, 23);            
    }



    size_t get_kmer(size_t p, char* kmer, char* rkmer) {
        // Get tf, kmer and rev_kmer stored in given arrays.
        uint64_t ukmer = hash_map->checker[p];
        uint64_t urev_kmer = reverseDNA(ukmer);
        get_bitset_dna23_c(ukmer, kmer, 23);
        get_bitset_dna23_c(urev_kmer, rkmer, 23);
        return hash_map->tf_values[p];
    }

    size_t get_hash_size() {
        return n;
    }

    void increase(char* ckmer) {
        std::string kmer = std::string(ckmer);
        hash_map->increase(kmer);
    }

    void decrease(char* ckmer) {
        std::string kmer = std::string(ckmer);
        hash_map->decrease(kmer);
    }

    void set_positions(size_t* r, std::string kmer) {
        // Set read positions
        auto h1 = hash_map->get_pfid(kmer);
        size_t j = 0;
        for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {
            positions[i] = r[j];
            j += 1;
        }
    }

    void check_aindex() {

        for (size_t h1=0; h1<hash_map->n; ++h1) {
            size_t tf = hash_map->tf_values[h1];
            size_t xtf = 0;

            if (h1 && h1 % 1000000 == 0) {
                std::cout << "Completed: " << h1 << "/" << hash_map->n << std::endl;
            }

            for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {
                if (positions[i] == 0) {
                    break;
                }

                xtf += 1;

                size_t pos = positions[i]-1;

                char ckmer[Settings::K];

                std::memcpy(ckmer, &reads[pos], Settings::K);
                ckmer[Settings::K] = '\0';
                std::string data_kmer = std::string(ckmer);

                uint64_t h1_kmer = hash_map->checker[h1];
                std::string kmer = get_bitset_dna23(h1_kmer);
                if (data_kmer != kmer) {
                    size_t rh1 = reverseDNA(h1_kmer);
                    std::string rkmer = get_bitset_dna23(rh1);
                    if (data_kmer != rkmer) {
                        std::cout << h1 << " " << i << " " << tf << " " << xtf << " " <<  data_kmer << " " << kmer << " " << rkmer << std::endl;
                    }
                }
            }
            if (tf != xtf) {
                std::cout << tf << " " << xtf << std::endl;

            }
        }
    }

    void check_aindex_reads() {

        bool* used_reads = new bool[1];
        std::vector<Hit> hits;

        for (size_t h1=0; h1<hash_map->n; ++h1) {

            if (h1 && h1 % 1000000 == 0) {
                std::cout << "Completed: " << h1 << "/" << hash_map->n << std::endl;
            }

            uint64_t h1_kmer = hash_map->checker[h1];
            std::string kmer = get_bitset_dna23(h1_kmer);
            hits.clear();
            get_reads_se_by_kmer(kmer, h1, used_reads, hits);

            size_t max_pos = 0;

            for (auto hit: hits) {
                std::max(max_pos, hit.pos);
                std::string subkmer = hit.read.substr(hit.pos, Settings::K);
                assert(subkmer == kmer);
                std::cout << kmer << " " << subkmer << " " << h1 << " " << hash_map->tf_values[h1] << std::endl;
            }
        }
    }


    void get_reads_se_by_kmer(std::string const kmer, size_t h1, bool* used_reads, std::vector<Hit> &hits) {


        for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {

            if (positions[i] == 0) {
                break;
            }

            size_t position = positions[i] - 1;
            size_t start = get_rid(position);

            size_t end = start;
            size_t spring_pos = 0;

            size_t pos = position - start;
            std::string left_read;
            std::string right_read;

            while (true) {
                if (reads[end] == '\n') {
                    if (spring_pos > 0) {
                        char rkmer[end-spring_pos];
                        std::memcpy(rkmer, &reads[spring_pos+1], end-spring_pos-1);
                        rkmer[end-spring_pos-1] = '\0';
                        right_read = std::string(rkmer);
                    }
                    break;
                } else if (reads[end] == '~') {
                    char lkmer[end-start+1];
                    std::memcpy(lkmer, &reads[start], end-start);
                    lkmer[end-start] = '\0';
                    left_read = std::string(lkmer);
                    spring_pos = end;
                }
                end += 1;
            }

            size_t real_rid = start2rid[start];

            Hit hit;
            hit.rid = real_rid;
            hit.start = start;
            hit.pos = pos;
            hit.rev = 0;
            spring_pos = spring_pos - start;

            if (pos < spring_pos) {
                hit.read = left_read;
                hit.ori = 0;

                if (hit.read.substr(hit.pos, Settings::K) != kmer) {
                    std::string rleft_read = hit.read;
                    get_revcomp(hit.read, rleft_read);
                    hit.pos = hit.read.length() - pos - Settings::K;
                    if (rleft_read.substr(hit.pos, Settings::K) != kmer) {
                        std::cout << rleft_read << std::endl;
                        std::cout << left_read << std::endl;
                        std::cout << right_read << std::endl;
                        std::cout << kmer << " " << pos <<  std::endl;
                        continue;
                    }
                    hit.read = rleft_read;
                    hit.rev = 1;
//                    hit.ori = 1;
                }
            } else {

                if (hit.pos == spring_pos) {
                    hit.pos = hit.pos - spring_pos;
                    std::cout <<  left_read << std::endl;
                    std::cout <<  right_read << std::endl;
                    std::cout << kmer << std::endl;

                } else {
                    hit.pos = hit.pos - spring_pos - 1;
                }

                hit.read = right_read;
                hit.ori = 1;

                if (hit.read.substr(hit.pos, Settings::K) != kmer) {
                    std::string rright_read = hit.read;
                    get_revcomp(hit.read, rright_read);
                    hit.pos = hit.read.length() - hit.pos - Settings::K;
                    if (rright_read.substr(hit.pos, Settings::K) != kmer) {
                        continue;
                    }
                    hit.read = rright_read;
                    hit.rev = 1;
//                    hit.ori = 0;
                }

            }

//            if (reversed_reads[real_rid]) {
//                if (hit.ori) {
//                    hit.ori = 0;
//                } else {
//                    hit.ori = 1;
//                }
//            }

            if (used_reads[2*hit.rid+hit.ori]) {
                continue;
            }
            hits.push_back(hit);

        }
    }


    void freeme(char* ptr)
    {
        std::cout << "freeing address: " << ptr << std::endl;
        free(ptr);
    }

};

extern "C" {

    AindexWrapper* AindexWrapper_new(){ return new AindexWrapper(); }
    void AindexWrapper_load(AindexWrapper* foo, char* index_prefix){ foo->load(index_prefix); }


    void AindexWrapper_freeme(AindexWrapper* foo, char* ptr){ foo->freeme(ptr); }


    void AindexWrapper_load_hash_file(AindexWrapper* foo, char* hash_filename){ foo->load(hash_filename); }

    void AindexWrapper_load_reads(AindexWrapper* foo, char* reads_file){ foo->load_reads(reads_file); }

    void AindexWrapper_load_index(AindexWrapper* foo, char* index_prefix, uint32_t max_tf){ foo->load_index(index_prefix, max_tf); }
    
    void AindexWrapper_increase(AindexWrapper* foo, char* kmer){ foo->increase(kmer); }
    void AindexWrapper_decrease(AindexWrapper* foo, char* kmer){ foo->decrease(kmer); }

    size_t AindexWrapper_get_kid_by_kmer(AindexWrapper* foo, char* kmer){ return foo->get_kid_by_kmer(kmer); }

    void AindexWrapper_get_kmer_by_kid(AindexWrapper* foo, size_t kid, char* kmer){ foo->get_kmer_by_kid(kid, kmer); }

    size_t AindexWrapper_get(AindexWrapper* foo, char* kmer){ return foo->get(kmer); }

    size_t AindexWrapper_get_n(AindexWrapper* foo){ return foo->get_n(); }

    size_t AindexWrapper_get_rid(AindexWrapper* foo, size_t pos){ return foo->get_rid(pos); }

//    char* AindexWrapper_get_read(AindexWrapper* foo, size_t start, int ori){ return foo->get_read(pos, ori); }

    void AindexWrapper_get_positions(AindexWrapper* foo, size_t* r, char* kmer){ foo->get_positions(r, kmer); }

    void AindexWrapper_set_positions(AindexWrapper* foo, size_t* r, char* kmer){ foo->set_positions(r, kmer); }

    size_t AindexWrapper_get_kmer(AindexWrapper* foo, size_t p, char* kmer, char* rkmer){ return foo->get_kmer(p, kmer, rkmer); }

    size_t AindexWrapper_get_hash_size(AindexWrapper* foo){ return foo->get_hash_size(); }
}