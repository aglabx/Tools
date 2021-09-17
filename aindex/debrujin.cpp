//
// Created by Aleksey Komissarov on 30/08/15.
//

#include "debrujin.hpp"
#include "kmers.hpp"
#include <stdint.h>
#include "hash.hpp"
#include "read.hpp"

namespace DEBRUJIN {

    int get_freq(uint64_t kmer, PHASH_MAP &kmers) {

        return kmers.get_freq(kmer);
    }


    void print_next(uint64_t kmer, PHASH_MAP &kmers, CONT &cont, unsigned int cutoff = 0) {
        /*
        */
        cont.n = 0;
        uint64_t kmer_A = ((kmer << 2) | ((uint64_t) 0 & 3)) & 0x00003fffffffffff;
        uint64_t kmer_C = ((kmer << 2) | ((uint64_t) 1 & 3)) & 0x00003fffffffffff;
        uint64_t kmer_G = ((kmer << 2) | ((uint64_t) 2 & 3)) & 0x00003fffffffffff;
        uint64_t kmer_T = ((kmer << 2) | ((uint64_t) 3 & 3)) & 0x00003fffffffffff;

        cont.A = kmers.get_freq(kmer_A);
        cont.C = kmers.get_freq(kmer_C);
        cont.G = kmers.get_freq(kmer_G);
        cont.T = kmers.get_freq(kmer_T);

        if (cutoff > 0) {
            if (cont.A <= cutoff) cont.A = 0;
            if (cont.C <= cutoff) cont.C = 0;
            if (cont.G <= cutoff) cont.G = 0;
            if (cont.T <= cutoff) cont.T = 0;
        }

        cont.sum = cont.A + cont.C + cont.G + cont.T;
        cont.n = (unsigned int) (bool) cont.A + (unsigned int) (bool) cont.C + (unsigned int) (bool) cont.G +
                 (unsigned int) (bool) cont.T;

        if (cont.A >= cont.C && cont.A >= cont.G && cont.A >= cont.T) {
            cont.best_hit = 'A';
            cont.best_ukmer = kmer_A;
        }
        if (cont.C >= cont.A && cont.C >= cont.G && cont.C >= cont.T) {
            cont.best_hit = 'C';
            cont.best_ukmer = kmer_C;
        }
        if (cont.G >= cont.C && cont.G >= cont.A && cont.G >= cont.T)  {
            cont.best_hit = 'G';
            cont.best_ukmer = kmer_G;
        }
        if (cont.T >= cont.C && cont.T >= cont.G && cont.T >= cont.A)  {
            cont.best_hit = 'T';
            cont.best_ukmer = kmer_T;
        }
    }


    void print_prev(uint64_t kmer, PHASH_MAP &kmers, CONT &cont, unsigned int cutoff = 0) {
        /*
        */
        cont.n = 0;
        uint64_t kmer_A = (kmer >> 2) | (((uint64_t) 0 & 3) << 44);
        uint64_t kmer_C = (kmer >> 2) | (((uint64_t) 1 & 3) << 44);
        uint64_t kmer_G = (kmer >> 2) | (((uint64_t) 2 & 3) << 44);
        uint64_t kmer_T = (kmer >> 2) | (((uint64_t) 3 & 3) << 44);

        cont.A = kmers.get_freq(kmer_A);
        cont.C = kmers.get_freq(kmer_C);
        cont.G = kmers.get_freq(kmer_G);
        cont.T = kmers.get_freq(kmer_T);

        if (cutoff > 0) {
            if (cont.A <= cutoff) cont.A = 0;
            if (cont.C <= cutoff) cont.C = 0;
            if (cont.G <= cutoff) cont.G = 0;
            if (cont.T <= cutoff) cont.T = 0;
        }

        cont.sum = cont.A + cont.C + cont.G + cont.T;

        if (cont.A >= cont.C && cont.A >= cont.G && cont.A >= cont.T) {
            cont.best_hit = 'A';
            cont.best_ukmer = kmer_A;
        }
        if (cont.C >= cont.A && cont.C >= cont.G && cont.C >= cont.T) {
            cont.best_hit = 'C';
            cont.best_ukmer = kmer_C;
        }
        if (cont.G >= cont.C && cont.G >= cont.A && cont.G >= cont.T)  {
            cont.best_hit = 'G';
            cont.best_ukmer = kmer_G;
        }
        if (cont.T >= cont.C && cont.T >= cont.G && cont.T >= cont.A)  {
            cont.best_hit = 'T';
            cont.best_ukmer = kmer_T;
        }

        cont.n = (unsigned int) (bool) cont.A + (unsigned int) (bool) cont.C + (unsigned int) (bool) cont.G +
                 (unsigned int) (bool) cont.T;
    }


    void set_fm_for_read(READS::READ &read, PHASH_MAP &kmers) {
        std::string kmer;
        for (size_t i = 0; i < read.seq.length() - Settings::K + 1; i++) {
            kmer = read.seq.substr(i, 23);
            read.fm[i] = kmers.get_freq(kmer);
        }
    }

    void set_fm_for_read(READS::READ &read, PHASH_MAP &kmers, int from, size_t to) {
        std::string kmer;
        if (from < 0) {
            from = 0;
        }
        if (to > read.seq.length()) {
            to = read.seq.length();
        }
        for (int i = from; i < to; i++) {
            kmer = read.seq.substr(i, Settings::K);
            read.fm[i] = kmers.get_freq(kmer);
        }
    }
}