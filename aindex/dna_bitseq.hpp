//
// Created by Aleksey Komissarov on 07/07/16.
//

#include <string>
#include <iostream>
#include "kmers.hpp"

#ifndef STIRKA_DNA_BITSEQ_HPP
#define STIRKA_DNA_BITSEQ_HPP

class dna_bitset {
public:
    /**
     * @brief constructor
     * @param dna_str a string containing a DNA sequence
     * @param dna_len length of the DNA sequence
     */
    dna_bitset (const char* dna_str, size_t dna_len)
    {
        m_len = dna_len;

        /* bytes necessary to store dna_str as a bitset */
        size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);

        m_data = new uint8_t[dna_bytes];

        std::memset(m_data, 0, dna_bytes);

        /* for each base of the DNA sequence */
        for (size_t i = 0; i < dna_len; i++)
        {
            uint8_t shift = 6 - 2*(i % 4);

            switch (dna_str[i])
            {
                case 'A':
                    m_data[i/4] |= BASE_A << shift;
                    break;
                case 'C':
                    m_data[i/4] |= BASE_C << shift;
                    break;
                case 'G':
                    m_data[i/4] |= BASE_G << shift;
                    break;
                case 'T':
                    m_data[i/4] |= BASE_T << shift;
                    break;
                default:
                    m_data[i/4] |= BASE_A << shift;
                    break;
            }

            shift = (shift == 0) ? 6 : shift-2;
        }
    }

    void kmer(size_t pos, int k, char* dna_str) {

        for (int i=0; i < k; ++i) {
            uint8_t shift = 6 - 2 * ((pos + i) % 4);
            uint8_t mask = BASE_MASK << shift;
            uint8_t base = (m_data[(i + pos) / 4] & mask) >> shift;
            switch (base) {
                case BASE_A:
                    dna_str[i] = 'A';
                    break;
                case BASE_C:
                    dna_str[i] = 'C';
                    break;
                case BASE_G:
                    dna_str[i] = 'G';
                    break;
                case BASE_T:
                    dna_str[i] = 'T';
                    break;
                default:
                    throw std::runtime_error("invalid DNA base");
            }
        }
    }

    char at(size_t pos) {

        uint8_t shift = 6 - 2 * ((pos) % 4);
        uint8_t mask = BASE_MASK << shift;
        uint8_t base = (m_data[(pos) / 4] & mask) >> shift;
        switch (base) {
            case BASE_A:
                return 'A';
                break;
            case BASE_C:
                return 'C';
                break;
            case BASE_G:
                return 'G';
                break;
            case BASE_T:
                return 'T';
                break;
            default:
                throw std::runtime_error("invalid DNA base");
        }

    }

    void cut_start(size_t pos) {

    }

    void cut_end(size_t pos) {

    }

    uint64_t ukmer(size_t pos, int k) {

        uint64_t num = 0;
        for (int i=0; i < k; ++i) {
            uint8_t shift = 6 - 2 * ((pos + i) % 4);
            uint8_t mask = BASE_MASK << shift;
            uint8_t base = (m_data[(i + pos) / 4] & mask) >> shift;
            num = num << 2;
            switch (base) {
                case BASE_A:
                    num += 0;
                    break;
                case BASE_C:
                    num += 1;
                    break;
                case BASE_G:
                    num += 2;
                    break;
                case BASE_T:
                    num += 3;
                    break;
                default:
                    throw std::runtime_error("invalid DNA base");
            }
        }
        return num;

    }

    /**
     * @brief destructor
     */
    ~dna_bitset ()
    {
        delete[] m_data;
    }

    /**
     * @brief returns the stored DNA sequence as a string
     */
    char* to_string() const
    {
        char* dna_str = new char[m_len+1];

        /* for each base of the DNA sequence */
        for (size_t i = 0; i < m_len; i++)
        {
            uint8_t shift = 6 - 2*(i % 4);

            uint8_t mask = BASE_MASK << shift;

            /* get the i-th DNA base */
            uint8_t base = (m_data[i/4] & mask) >> shift;

            switch (base)
            {
                case BASE_A:
                    dna_str[i] = 'A';
                    break;
                case BASE_C:
                    dna_str[i] = 'C';
                    break;
                case BASE_G:
                    dna_str[i] = 'G';
                    break;
                case BASE_T:
                    dna_str[i] = 'T';
                    break;
                default:
                    throw std::runtime_error("invalid DNA base");
            }
        }

        dna_str[m_len] = '\0';
        return dna_str;
    }

private:
    uint8_t* m_data;
    size_t   m_len;
};

#endif //STIRKA_DNA_BITSEQ_HPP