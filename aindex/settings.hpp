//
// Created by Aleksey Komissarov on 02/09/15.
//

#ifndef STIRKA_SETTINGS_H
#define STIRKA_SETTINGS_H

#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <memory>
#include <ctime>
#include <cstring>

namespace Settings {

    extern unsigned int MINIMAL_READ_LENGTH;
    extern unsigned int TRUE_ERRORS;
    extern unsigned int PRE_ERRORS;
    extern unsigned int K;
    extern unsigned int MIN_Q;
    extern unsigned int MIN_MI;
    extern unsigned int MIN_MI_N;
    extern unsigned int MAX_MA;
    extern unsigned int TRIM_LOW_COV_TOLERANCE;
    extern unsigned int STIRKA_FILL_OVERLAP;
    extern unsigned int STIRKA_RESOLVE_SNP;

}
#endif //STIRKA_SETTINGS_H
