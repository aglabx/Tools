//
// Created by Aleksey Komissarov on 30/08/15.
//
#include <stdint.h>
#include <string>
#include <limits.h>
#include <iostream>
#include "settings.hpp"

/// CONVERTERS to uint 23-mers and 13-mers from strings and char*

uint64_t get_dna23_bitset(std::string dna_str) {
    /*
     * Convert 23-mer to bit 23-mer.
     */
    uint64_t num = 0;
    for (int8_t n=0; n<Settings::K; n++) {
        num = num << 2;
        if (dna_str[n] == 'A') num += 0;
        if (dna_str[n] == 'C') num += 1;
        if (dna_str[n] == 'G') num += 2;
        if (dna_str[n] == 'T') num += 3;
    }
    return num;
}

uint32_t get_dna13_bitset(std::string dna_str) {
    /*
     * Convert 13-mer to bit 13-mer.
     */
    uint32_t num = 0;
    for (int8_t n=0; n<13; n++) {
        num = num << 2;
        if (dna_str[n] == 'A') num += 0;
        if (dna_str[n] == 'C') num += 1;
        if (dna_str[n] == 'G') num += 2;
        if (dna_str[n] == 'T') num += 3;
    }
    return num;
}

uint64_t get_dna23_bitset(char* dna_str) {
    /*
     * Convert 23-mer to bit 23-mer.
     */
    uint64_t num = 0;
    for (int8_t n=0; n<Settings::K; n++) {
        num = num << 2;
        if (dna_str[n] == 'A') num += 0;
        if (dna_str[n] == 'C') num += 1;
        if (dna_str[n] == 'G') num += 2;
        if (dna_str[n] == 'T') num += 3;
    }
    return num;
}

uint32_t get_dna13_bitset(char* dna_str) {
    /*
     * Convert 13-mer to bit 13-mer.
     */
    uint32_t num = 0;
    for (int8_t n=0; n<13; n++) {
        num = num << 2;
        if (dna_str[n] == 'A') num += 0;
        if (dna_str[n] == 'C') num += 1;
        if (dna_str[n] == 'G') num += 2;
        if (dna_str[n] == 'T') num += 3;
    }
    return num;
}

/// CONVERTERS from uint to string kmers

void get_bitset_dna23(uint64_t x, std::string &res, int k) {
    /*
     * Convert bit 23-mer to string.
     */
    for (int8_t i = k-1; i+1; i--) {
        switch(x & 3) {
            case 0:
                res[i] = 'A';
                break;

            case 1:
                res[i] = 'C';
                break;

            case 2:
                res[i] = 'G';
                break;

            case 3:
                res[i] = 'T';
                break;
        }

        x >>= 2;
    }
}

void get_bitset_dna23_c(uint64_t x, char *res, int k) {
    /*
     * Convert bit 23-mer to string.
     */
    for (int8_t i = k-1; i+1; i--) {
        switch(x & 3) {
            case 0:
                res[i] = 'A';
                break;

            case 1:
                res[i] = 'C';
                break;

            case 2:
                res[i] = 'G';
                break;

            case 3:
                res[i] = 'T';
                break;
        }

        x >>= 2;
    }
}

std::string get_bitset_dna23(uint64_t x) {
    /*
     * Convert bit 23-mer to string.
     */
    std::string res = "NNNNNNNNNNNNNNNNNNNNNNN";
    for (int8_t i = 22; i+1; i--) {
        switch(x & 3) {
            case 0:
                res[i] = 'A';
                break;

            case 1:
                res[i] = 'C';
                break;

            case 2:
                res[i] = 'G';
                break;

            case 3:
                res[i] = 'T';
                break;
        }

        x >>= 2;
    }

    return res;

}

void get_bitset_dna13(uint32_t x, std::string &res, int k) {
    /*
     * Convert bit 13-mer to string.
     */
    for (int8_t i = k-1; i+1; i--) {
        switch(x & 3) {
            case 0:
                res[i] = 'A';
                break;

            case 1:
                res[i] = 'C';
                break;

            case 2:
                res[i] = 'G';
                break;

            case 3:
                res[i] = 'T';
                break;
        }

        x >>= 2;
    }
}

void get_bitset_dna13_c(uint32_t x, char *res, int k) {
    /*
     * Convert bit 13-mer to string.
     */
    for (int8_t i = k-1; i+1; i--) {
        switch(x & 3) {
            case 0:
                res[i] = 'A';
                break;

            case 1:
                res[i] = 'C';
                break;

            case 2:
                res[i] = 'G';
                break;

            case 3:
                res[i] = 'T';
                break;
        }

        x >>= 2;
    }
}

std::string get_bitset_dna13(uint32_t x) {
    /*
     * Convert bit 13-mer to string.
     */
    std::string res = "NNNNNNNNNNNNN";
    for (int8_t i = 12; i+1; i--) {
        switch(x & 3) {
            case 0:
                res[i] = 'A';
                break;

            case 1:
                res[i] = 'C';
                break;

            case 2:
                res[i] = 'G';
                break;

            case 3:
                res[i] = 'T';
                break;
        }

        x >>= 2;
    }

    return res;

}

/// REVERSE COMPLEMENT SECTION

//unsigned int reverse(uint64_t t) {
//    /*
//     * Reverse complement for bit 23-mer.
//     */
//    size_t n = t; //store in 64 bit number for call to BSWAP
//    __asm__("BSWAP %0" : "=r"(n) : "0"(n));
//    n >>= ((sizeof(size_t) - sizeof(uint64_t)) * 8);
//    n = ((n & 0xaaaaaaaa) >> 1) | ((n & 0x55555555) << 1);
//    n = ((n & 0xcccccccc) >> 2) | ((n & 0x33333333) << 2);
//    n = ((n & 0xf0f0f0f0) >> 4) | ((n & 0x0f0f0f0f) << 4);
//    return n;
//}
//
//uint32_t reverse(uint32_t t) {
//    /*
//     * Reverse complement for bit 13-mer.
//     */
//    size_t n = t; //store in 64 bit number for call to BSWAP
//    __asm__("BSWAP %0" : "=r"(n) : "0"(n));
//    n >>= ((sizeof(size_t) - sizeof(uint32_t)) * 8);
//    n = ((n & 0xaaaaaaaa) >> 1) | ((n & 0x55555555) << 1);
//    n = ((n & 0xcccccccc) >> 2) | ((n & 0x33333333) << 2);
//    n = ((n & 0xf0f0f0f0) >> 4) | ((n & 0x0f0f0f0f) << 4);
//    return (uint32_t)n;
//}


void get_revcomp(std::string &input, std::string &output) {
    /*
     * Revcomp for string.
     */
    size_t n = input.length();
    for(int y = n-1; y >= 0; y--) {
        if (input[y] == 'A') {
            output[n-1-y] = 'T';
        } else if (input[y] == 'C') {
            output[n-1-y] = 'G';
        } else if (input[y] == 'G') {
            output[n-1-y] = 'C';
        } else if (input[y] == 'T') {
            output[n-1-y] = 'A';
        } else if (input[y] == '~') {
            output[n-1-y] = '~';
        } else {
            output[n-1-y] = 'N';
        }
    }
}

std::string get_revcomp(std::string &input) {
    /*
     * Revcomp for  string.
     */
    size_t n = input.length();
    std::string output(n, 'N');
    for(int y = n-1; y >= 0; y--) {
        if (input[y] == 'A') {
            output[n-1-y] = 'T';
        } else if (input[y] == 'C') {
            output[n-1-y] = 'G';
        } else if (input[y] == 'G') {
            output[n-1-y] = 'C';
        } else if (input[y] == 'T') {
            output[n-1-y] = 'A';
        } else {
            output[n-1-y] = 'N';
        }
    }
    return output;
}

uint64_t _reversePairs(uint64_t num) {
    /*
     * Reverse bit 23-mer helper.
     */
    uint64_t count = CHAR_BIT * sizeof num - 2;
    uint64_t reverse_num = num;
    for (num >>= 2; num; reverse_num <<= 2, reverse_num |= num & 3, num >>= 2, count -= 2);
    return reverse_num << count;
}

uint32_t _reversePairs(uint32_t num) {
    /*
     * Reverse bit 13-mer helper.
     */
    uint32_t count = CHAR_BIT * sizeof num - 2;
    uint32_t reverse_num = num;
    for (num >>= 2; num; reverse_num <<= 2, reverse_num |= num & 3, num >>= 2, count -= 2);
    return reverse_num << count;
}


uint64_t reverseDNA(uint64_t num) {
    /*
     * Reverse bit 23-mer.
     */
    return (~_reversePairs(num) >> 18);
}

uint32_t reverseDNA(uint32_t num) {
    /*
     * Reverse bit 13-mer.
     */
    return (~_reversePairs(num) >> 6);
}
