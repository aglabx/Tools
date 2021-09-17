//
// Created by Aleksey Komissarov on 30/08/15.
//

#ifndef STIRKA_KMERS_H
#define STIRKA_KMERS_H

#include <stdint.h>
#include <string>


#define BASE_MASK 0x3 /* binary: 11 */

/* useful constants */
enum {
    BASE_A = 0x0, /* binary: 00 */
    BASE_C = 0x1, /*'binary: 01 */
    BASE_G = 0x2, /* binary: 10 */
    BASE_T = 0x3, /* binary: 11 */
};

void get_bitset_dna23(uint64_t x, std::string &res, int k=23);
void get_bitset_dna23_c(uint64_t x, char *res, int k);
std::string get_bitset_dna23(uint64_t x);

void get_bitset_dna13(uint32_t x, std::string &res, int k=13);
void get_bitset_dna13_c(uint32_t x, char *res, int k);
std::string get_bitset_dna13(uint32_t x);

uint64_t get_dna23_bitset(std::string dna_str);
uint64_t get_dna23_bitset(char* dna_str);

uint32_t get_dna13_bitset(std::string dna_str);
uint32_t get_dna13_bitset(char* dna_str);

void get_revcomp(std::string &input, std::string &output);
std::string get_revcomp(std::string &input);

uint64_t reverseDNA(uint64_t num);
uint32_t reverseDNA(uint32_t num);

#endif //STIRKA_KMERS_H
