//
// Created by Aleksey Komissarov on 30/08/15.
//

#ifndef STIRKA_HASH_H
#define STIRKA_HASH_H

#include <stdint.h>
#include <unordered_map>
#include <string>
#include "emphf/common.hpp"
#include "emphf/perfutils.hpp"
#include "emphf/mphf.hpp"
#include "emphf/base_hash.hpp"
#include <atomic>
#include "kmers.hpp"
#include <stdint.h>
#include "settings.hpp"
#include <mutex>
#include <thread>
#include <functional>
#include <vector>
#include <algorithm>

typedef emphf::mphf<emphf::jenkins64_hasher> HASHER;
typedef int *VAULT;

typedef std::atomic<unsigned int> ATOMIC;
typedef std::atomic<unsigned long> ATOMIC_LONG;

typedef std::atomic<uint8_t> ATOMIC8;
typedef std::atomic<size_t> ATOMIC64;

typedef std::unordered_map<uint64_t, int> HASH_MAP;
typedef std::unordered_map<std::string, int> HASH_MAP13;


struct Stats {

    size_t zero = 0;
    size_t unique = 0;
    size_t distinct = 0;
    size_t total = 0;
    size_t max_count = 0;
    size_t coverage = 0;

    int *profile;

    Stats() {
        profile = nullptr;
    }

    void init(int coverage) {
        zero = 0;
        unique = 0;
        distinct = 0;
        total = 0;
        max_count = 0;
        if (profile != nullptr) {
            delete [] profile;
            profile = nullptr;
        }
        profile = new int[coverage+coverage/2];
        for (size_t i=0; i<coverage+coverage/2; i++) {
            profile[i] = 0;
        }

//
//        for (size_t i=0; i < coverage+coverage/2; ++i) {
//            std::cout << i << " " << profile[i] << std::endl;
//        }
    }

    ~Stats() {
        if (profile != nullptr) {
            delete [] profile;
        }

    }
};



struct PHASH_MAP {

    HASHER hasher;
    ATOMIC *tf_values;
    ATOMIC_LONG *left_qtf_values;
    ATOMIC_LONG *right_qtf_values;
    uint64_t *checker;
    size_t n = 0;

    Stats stats;

    emphf::stl_string_adaptor str_adapter;

    PHASH_MAP() {
        tf_values = nullptr;
        left_qtf_values = nullptr;
        right_qtf_values = nullptr;
        checker = nullptr;
        n = 0;
    }

    size_t get_n() {
        return n;
    }

    ~PHASH_MAP() {
        if (tf_values != nullptr) {
            delete [] tf_values;
            tf_values = nullptr;
        }
        if (left_qtf_values != nullptr) delete [] left_qtf_values;
        if (right_qtf_values != nullptr) delete [] right_qtf_values;
        if (checker != nullptr) delete [] checker;

    }

    inline unsigned int get_freq(uint64_t kmer) {
        std::string _kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        std::string _rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        get_bitset_dna23(kmer, _kmer);
        auto h1 = hasher.lookup(_kmer, str_adapter);
        // try to find forward
        if (h1 < n && checker[h1] == kmer) {
            return tf_values[h1].load();
        }
        // else try to find reverse
        uint64_t rev_kmer = reverseDNA(kmer);
        get_bitset_dna23(rev_kmer, _rev_kmer);
        auto h2 = hasher.lookup(_rev_kmer, str_adapter);
        if (h2 < n && checker[h2] == rev_kmer) {
            return tf_values[h2].load();
        }
        return 0;
    }


    inline size_t get_index_unsafe(std::string kmer) {
        return hasher.lookup(kmer, str_adapter);
    }

    inline size_t get_pfid(std::string &_kmer) {
        uint64_t kmer = get_dna23_bitset(_kmer);
        std::string _rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        uint64_t rev_kmer = reverseDNA(kmer);
        get_bitset_dna23(rev_kmer, _rev_kmer);
        if (_kmer.compare(_rev_kmer) <= 0) {
            size_t h1 = hasher.lookup(_kmer, str_adapter);
            if (h1 < n && checker[h1] == kmer) {
                return h1;
            } else {
                return n;
            }
        } else {
            size_t h1 = hasher.lookup(_rev_kmer, str_adapter);
            if (h1 < n && checker[h1] == rev_kmer) {
                return h1;
            } else {
                return n;
            }
        }
    }

    inline size_t get_pfid_by_umer_safe(uint64_t kmer) {

        std::string _kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        get_bitset_dna23(kmer, _kmer, Settings::K);

        std::string _rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        uint64_t rev_kmer = reverseDNA(kmer);
        get_bitset_dna23(rev_kmer, _rev_kmer);
        if (_kmer.compare(_rev_kmer) <= 0) {
            size_t h1 = hasher.lookup(_kmer, str_adapter);
            if (h1 < n && checker[h1] == kmer) {
                return h1;
            } else {
                return n;
            }
        } else {
            size_t h1 = hasher.lookup(_rev_kmer, str_adapter);
            if (h1 < n && checker[h1] == rev_kmer) {
                return h1;
            } else {
                return n;
            }
        }
    }

    inline size_t get_pfid_by_umer_unsafe(uint64_t kmer) {
        std::string _kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        get_bitset_dna23(kmer, _kmer);
        return hasher.lookup(_kmer, str_adapter);
    }

    inline unsigned int get_freq(std::string &kmer) {
        uint64_t _kmer = get_dna23_bitset(kmer);
        return get_freq(_kmer);
    }

    inline uint64_t get_kmer(size_t p) {
        return checker[p];
    }

    inline ATOMIC& get_atomic(uint64_t kmer) {
        std::string _kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
        get_bitset_dna23(kmer, _kmer);
        return tf_values[hasher.lookup(_kmer, str_adapter)];
    }

    inline void increase(std::string &kmer) {
        auto h1 = hasher.lookup(kmer, str_adapter);
        uint64_t _kmer = get_dna23_bitset(kmer);
        if (h1 < n && checker[h1] == _kmer) {
            tf_values[h1]++;
        } else {
            std::string _rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
            uint64_t rev_kmer = reverseDNA(_kmer);
            get_bitset_dna23(rev_kmer, _rev_kmer);
            auto h2 = hasher.lookup(_rev_kmer, str_adapter);
            if (h2 < n && checker[h2] == rev_kmer) {
                tf_values[h2]++;
            }
        }
    }

    inline void increase_raw(std::string &kmer) {
        auto h1 = hasher.lookup(kmer, str_adapter);
        tf_values[h1]++;
    }

    inline void decrease(std::string &kmer) {
        auto h1 = hasher.lookup(kmer, str_adapter);
        uint64_t _kmer = get_dna23_bitset(kmer);
        if (h1 < n && checker[h1] == _kmer && tf_values[h1].load()>0) {
            if (tf_values[h1] > 0) tf_values[h1]--;
        } else {
            std::string _rev_kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
            uint64_t rev_kmer = reverseDNA(_kmer);
            get_bitset_dna23(rev_kmer, _rev_kmer);
            auto h2 = hasher.lookup(_rev_kmer, str_adapter);
            if (h2 < n && checker[h2] == rev_kmer) {
                if (tf_values[h2] > 0) tf_values[h2]--;
            }
        }
    }

    void save_values(std::string &file_name) {
        std::ofstream fh(file_name);

        if (!fh) {
            std::cerr << "Cannot open file for writting.";
            exit(12);
        }

        for (size_t i=0; i < n; i++) {
            std::string kmer = "NNNNNNNNNNNNNNNNNNNNNNN";
            get_bitset_dna23(checker[i], kmer);
            fh << kmer << "\t" << tf_values[i].load() << "\n";
        }
        fh.close();

    }

    void reset_tf() {
        for (size_t i=0; i < n; i++) {
            tf_values[i] = 0;
        }
    }

    void set_stats(int coverage) {

        stats.init(coverage);

        int max_coverage = coverage + coverage/2;

        for (size_t i=0; i < n; i++) {
            stats.total += tf_values[i];
            if (tf_values[i] == 0) {
                stats.zero += 1;
            }
            if (tf_values[i] == 1) {
                stats.unique += 1;
            }
            if (tf_values[i] > 0) {
                stats.distinct += 1;
            }
            if (tf_values[i] < max_coverage) {
                stats.profile[tf_values[i]] += 1;
            } else {
                stats.profile[max_coverage-1] += 1;
            }
            if (tf_values[i] > stats.max_count) {
                stats.max_count = tf_values[i];
            }
        }
    }

    void print_stats() {

        std::cout << "Zero: " << stats.zero << std::endl;
        std::cout << "Unique: " << stats.unique << std::endl;
        std::cout << "Distinct: " << stats.distinct << std::endl;
        std::cout << "Total: " << stats.total << std::endl;
        std::cout << "Coverage: " << stats.coverage << std::endl;
        std::cout << "Max value: " << stats.max_count << std::endl;

    }

    void print_stats_profile(int coverage) {
        for (size_t i=0; i<coverage+coverage/2;i++) {
            std::cout << i << ":" << stats.profile[i] << " ";
        }
        std::cout << std::endl;
    }

    std::string print_and_set_coverage(int coverage) {
        set_stats(coverage);
        print_stats_profile(coverage);
        std::string res = "Z: " + std::to_string(stats.zero) + " U: " + std::to_string(stats.unique) + " D: " + std::to_string(stats.distinct) + " T: " + std::to_string(stats.total) + " C: " + std::to_string(stats.coverage) + " M: " + std::to_string(stats.max_count);
        std::cout << res << std::endl;
        return res;
    }

private:

};


void lu_compressed_worker(int worker_id, size_t start, size_t end, char *contents,  ATOMIC64 *positions, ATOMIC64 *ppositions, size_t* indices, PHASH_MAP &hash_map);


struct AIndexCompressed {

    size_t* indices; // position indices
    ATOMIC64* ppositions; // position completness
    ATOMIC64* positions; // position itself
    size_t total_size = 0;
    size_t max_tf = 0;

    AIndexCompressed(PHASH_MAP &hash_map) {

        std::cout << "Allocate indices" << std::endl;
        indices = new size_t[hash_map.n+1];
        if (indices == nullptr) {
            emphf::logger() << "Failed to allocate memory for positions: " << hash_map.n+1 << std::endl;
            exit(10);
        }
        indices[0] = 0;
        for (size_t i=1; i<hash_map.n+1; ++i) {
//            if (hash_map.tf_values[i-1] > 1)
//                std::cout << indices[i-1] + hash_map.tf_values[i-1] << " " << hash_map.tf_values[i-1] << std::endl;
            indices[i] = indices[i-1] + hash_map.tf_values[i-1];
            total_size += hash_map.tf_values[i-1];
            max_tf = std::max(max_tf, (size_t)hash_map.tf_values[i-1]);
        }

        std::cout << "\tmax_tf: " << max_tf << std::endl;
        std::cout << "\ttotal_size: " << total_size << std::endl;

        std::cout << "Allocate ppositions" << std::endl;
        ppositions = new ATOMIC64[hash_map.n];
        if (ppositions == nullptr) {
            emphf::logger() << "Failed to allocate memory for positions: " << hash_map.n << std::endl;
            exit(10);
        }
        memset(ppositions, 0, hash_map.n * sizeof(ATOMIC64));

        std::cout << "Allocate positions" << std::endl;
        positions = new ATOMIC64[total_size];
        if (positions == nullptr) {
            emphf::logger() << "Failed to allocate memory for positions: " << total_size << std::endl;
            exit(10);
        }
        memset(positions, 0, total_size * sizeof(ATOMIC64));

    }

    ~AIndexCompressed() {
        if (indices != nullptr) delete [] indices;
        if (ppositions != nullptr) delete [] ppositions;
        if (positions != nullptr) delete [] positions;
    }

    void fill_index_from_reads(char *contents, size_t length, uint num_threads, PHASH_MAP &hash_map) {

        emphf::logger() << "Building index..." << std::endl;

        size_t batch_size = (length / num_threads) + 1;
        std::vector<std::thread> t;

        for (size_t worker_id = 0; worker_id < num_threads; ++worker_id) {
            size_t start = worker_id * batch_size;
            size_t end = (worker_id + 1) * batch_size;
            if (end > length) {
                end = length;
            }
            // inner worker takes kmers from start..end-k+1
            if (start > Settings::K) {
                start -= (Settings::K-1);
            }

            t.push_back(std::thread(
                    lu_compressed_worker,
                    worker_id,
                    start,
                    end,
                    contents,
                    std::ref(positions),
                    std::ref(ppositions),
                    std::ref(indices),
                    std::ref(hash_map)
            )
            );
        }

        for (size_t worker_id = 0; worker_id < num_threads; ++worker_id) {
            t[worker_id].join();
        }

        emphf::logger() << "\tDone." << std::endl;
    }

    void get_positions(std::string kmer, unsigned int* r, PHASH_MAP &hash_map) {

        memset(r, 0, max_tf * sizeof(unsigned int));
        auto h1 = hash_map.get_pfid(kmer);
        size_t j = 0;
        for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {
            r[j] = positions[i];
            j += 1;
        }
    }

    void set_positions(std::string kmer, unsigned int* r, PHASH_MAP &hash_map) {
        auto h1 = hash_map.get_pfid(kmer);
        size_t j = 0;
        for (size_t i=indices[h1]; i < indices[h1+1]; ++i) {
            positions[i] = r[j];
            j += 1;
        }
    }

    void save(std::string output_prefix, std::vector<size_t> start_positions, PHASH_MAP &hash_map) {
        //
        emphf::logger() << "Saving pos.bin array..." << std::endl;
        std::ofstream fout2(output_prefix + ".pos.bin", std::ios::out | std::ios::binary);
        fout2.write((char *) &start_positions[0], start_positions.size() * sizeof(size_t));
        fout2.close();

        emphf::logger() << "Saving index.bin array..." << std::endl;
        std::ofstream fout3(output_prefix+ ".index.bin", std::ios::out | std::ios::binary);
        emphf::logger() << "Positions array size: " << sizeof(size_t) * total_size << std::endl;
        fout3.write(reinterpret_cast<const char *> (positions), sizeof(size_t) * total_size);
        fout3.close();

        emphf::logger() << "Saving indices array..." << std::endl;
        std::ofstream fout4(output_prefix+ ".indices.bin", std::ios::out | std::ios::binary);
        emphf::logger() << "Indices array size: " << sizeof(size_t) * total_size << std::endl;
        fout4.write(reinterpret_cast<const char *> (indices), sizeof(size_t) * (hash_map.n+1));
        fout4.close();

        emphf::logger() << "\tDone." << std::endl;
    }

private:

};


struct AtomicCounter {

    std::atomic<unsigned int> value;

    void increment(){
        ++value;
    }

    void decrement(){
        --value;
    }

    int get(){
        return value.load();
    }
};

extern void load_hash(PHASH_MAP &hash_map, std::string &output_prefix, std::string &tf_file, std::string &hash_filename);
extern void load_only_hash(PHASH_MAP &hash_map, std::string &hash_filename);
void construct_hash_unordered_hash_illumina(std::string data_file, HASH_MAP13 &kmers);
void load_hash_for_qkmer(PHASH_MAP &hash_map, size_t n, std::string &data_filename, std::string &hash_filename);
void index_hash(PHASH_MAP &hash_map, std::string &dat_filename, std::string &hash_filename);
void index_hash_pp(PHASH_MAP &hash_map, std::string &dat_filename, std::string &hash_filename, int num_threads, int mock_dat=0);
void load_full_hash(PHASH_MAP &hash_map, std::string &hash_filename, int k, size_t n);
void load_hash_full_tf(PHASH_MAP &hash_map, std::string &output_prefix, std::string &tf_file, std::string &hash_filename);
#endif //STIRKA_HASH_H
