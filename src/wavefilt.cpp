/*
* Copyright (c) 2014, Rafat Hussain
* Copyright (C) 2016  Holger Nahrstaedt
* Daubechies wavelets coefficents Computed by Kazuo Hatano, Aichi Institute of Technology. http://phase.hpcc.jp/phase/wavelet/
*/
#include "wavefilt.h"

#include "filters/coif.hpp"
#include "filters/db.hpp"
#include "filters/h.hpp"
#include "filters/meyer.hpp"
#include "filters/sym.hpp"

#include <memory>
#include <string_view>

using namespace std::string_view_literals;

auto filtlength(char const* name) -> int
{
    int len = strlen(name);
    int i = 0;
    int N = 0;
    if (name == "haar"sv || name == "db1"sv) {
        return 2;
    }
    if (len > 2 && strstr(name, "db") != nullptr) {
        auto new_str = std::make_unique<char[]>((len - 2 + 1));
        for (i = 2; i < len + 1; i++) {
            new_str[i - 2] = name[i];
        }

        N = atoi(new_str.get());
        if (N > 38) {
            printf("\n Filter Not in Database \n");
            return -1;
        }

        return N * 2;
    }
    if (name == "bior1.1"sv) {
        return 2;
    }

    if (name == "bior1.3"sv) {
        return 6;
    }

    if (name == "bior1.5"sv) {
        return 10;
    }

    if (name == "bior2.2"sv) {
        return 6;
    }

    if (name == "bior2.4"sv) {
        return 10;
    }

    if (name == "bior2.6"sv) {
        return 14;
    }
    if (name == "bior2.8"sv) {
        return 18;
    }

    if (name == "bior3.1"sv) {
        return 4;
    }
    if (name == "bior3.3"sv) {
        return 8;
    }
    if (name == "bior3.5"sv) {
        return 12;
    }

    if (name == "bior3.7"sv) {
        return 16;
    }
    if (name == "bior3.9"sv) {
        return 20;
    }
    if (name == "bior4.4"sv) {
        return 10;
    }
    if (name == "bior5.5"sv) {
        return 12;
    }
    if (name == "bior6.8"sv) {
        return 18;
    }
    if (name == "rbior1.1"sv) {
        return 2;
    }

    if (name == "rbior1.3"sv) {
        return 6;
    }

    if (name == "rbior1.5"sv) {
        return 10;
    }

    if (name == "rbior2.2"sv) {
        return 6;
    }

    if (name == "rbior2.4"sv) {
        return 10;
    }

    if (name == "rbior2.6"sv) {
        return 14;
    }
    if (name == "rbior2.8"sv) {
        return 18;
    }

    if (name == "rbior3.1"sv) {
        return 4;
    }
    if (name == "rbior3.3"sv) {
        return 8;
    }
    if (name == "rbior3.5"sv) {
        return 12;
    }

    if (name == "rbior3.7"sv) {
        return 16;
    }
    if (name == "rbior3.9"sv) {
        return 20;
    }
    if (name == "rbior4.4"sv) {
        return 10;
    }
    if (name == "rbior5.5"sv) {
        return 12;
    }
    if (name == "rbior6.8"sv) {
        return 18;
    }
    if (len > 4 && strstr(name, "coif") != nullptr) {
        auto new_str = std::make_unique<char[]>((len - 4 + 1));
        for (i = 4; i < len + 1; i++) {
            new_str[i - 4] = name[i];
        }

        N = atoi(new_str.get());
        if (N > 17) {
            printf("\n Filter Not in Database \n");
            return -1;
        }

        return N * 6;
    }
    if (len > 3 && strstr(name, "sym") != nullptr) {
        auto new_str = std::make_unique<char[]>((len - 3 + 1));
        for (i = 3; i < len + 1; i++) {
            new_str[i - 3] = name[i];
        }

        N = atoi(new_str.get());
        if (N > 20 || N < 2) {
            printf("\n Filter Not in Database \n");
            return -1;
        }

        return N * 2;
    }
    if (strcmp(name, "meyer") == 0) {
        return 102;
    }
    printf("\n Filter Not in Database \n");
    return -1;
}

void copy_reverse(double const* in, int N, double* out)
{
    int count = 0;
    for (count = 0; count < N; count++) {
        out[count] = in[N - count - 1];
    }
}

void qmf_wrev(double const* in, int N, double* out)
{
    auto sigOutTemp = std::make_unique<double[]>(N);

    qmf_even(in, N, sigOutTemp.get());
    copy_reverse(sigOutTemp.get(), N, out);
}

void qmf_even(double const* in, int N, double* out)
{
    int count = 0;
    for (count = 0; count < N; count++) {
        out[count] = in[N - count - 1];
        if (count % 2 != 0) {
            out[count] = -1 * out[count];
        }
    }
}
void copy(double const* in, int N, double* out)
{
    int count = 0;
    for (count = 0; count < N; count++) {
        out[count] = in[count];
    }
}

auto filtcoef(char const* name, double* lp1, double* hp1, double* lp2, double* hp2) -> int
{
    int i = 0;
    int N = filtlength(name);
    if ((strcmp(name, "haar") == 0) || (strcmp(name, "db1") == 0)) {
        copy_reverse(db1, N, lp1);
        qmf_wrev(db1, N, hp1);
        copy(db1, N, lp2);
        qmf_even(db1, N, hp2);

        return N;
    }
    if (name == "db2"sv) {
        copy_reverse(db2, N, lp1);
        qmf_wrev(db2, N, hp1);
        copy(db2, N, lp2);
        qmf_even(db2, N, hp2);

        return N;
    }
    if (name == "db3"sv) {
        copy_reverse(db3, N, lp1);
        qmf_wrev(db3, N, hp1);
        copy(db3, N, lp2);
        qmf_even(db3, N, hp2);

        return N;
    }
    if (name == "db4"sv) {
        copy_reverse(db4, N, lp1);
        qmf_wrev(db4, N, hp1);
        copy(db4, N, lp2);
        qmf_even(db4, N, hp2);

        return N;
    }
    if (name == "db5"sv) {
        copy_reverse(db5, N, lp1);
        qmf_wrev(db5, N, hp1);
        copy(db5, N, lp2);
        qmf_even(db5, N, hp2);

        return N;
    }
    if (name == "db6"sv) {
        copy_reverse(db6, N, lp1);
        qmf_wrev(db6, N, hp1);
        copy(db6, N, lp2);
        qmf_even(db6, N, hp2);

        return N;
    }
    if (name == "db7"sv) {
        copy_reverse(db7, N, lp1);
        qmf_wrev(db7, N, hp1);
        copy(db7, N, lp2);
        qmf_even(db7, N, hp2);

        return N;
    }
    if (name == "db8"sv) {
        copy_reverse(db8, N, lp1);
        qmf_wrev(db8, N, hp1);
        copy(db8, N, lp2);
        qmf_even(db8, N, hp2);

        return N;
    }
    if (name == "db9"sv) {
        copy_reverse(db9, N, lp1);
        qmf_wrev(db9, N, hp1);
        copy(db9, N, lp2);
        qmf_even(db9, N, hp2);

        return N;
    }

    if (name == "db10"sv) {
        copy_reverse(db10, N, lp1);
        qmf_wrev(db10, N, hp1);
        copy(db10, N, lp2);
        qmf_even(db10, N, hp2);

        return N;
    }

    if (name == "db11"sv) {
        copy_reverse(db11, N, lp1);
        qmf_wrev(db11, N, hp1);
        copy(db11, N, lp2);
        qmf_even(db11, N, hp2);

        return N;
    }
    if (name == "db12"sv) {
        copy_reverse(db12, N, lp1);
        qmf_wrev(db12, N, hp1);
        copy(db12, N, lp2);
        qmf_even(db12, N, hp2);

        return N;
    }

    if (name == "db13"sv) {
        copy_reverse(db13, N, lp1);
        qmf_wrev(db13, N, hp1);
        copy(db13, N, lp2);
        qmf_even(db13, N, hp2);

        return N;
    }

    if (name == "db14"sv) {
        copy_reverse(db14, N, lp1);
        qmf_wrev(db14, N, hp1);
        copy(db14, N, lp2);
        qmf_even(db14, N, hp2);

        return N;
    }

    if (name == "db15"sv) {
        copy_reverse(db15, N, lp1);
        qmf_wrev(db15, N, hp1);
        copy(db15, N, lp2);
        qmf_even(db15, N, hp2);

        return N;
    }

    if (name == "db16"sv) {
        copy_reverse(db16, N, lp1);
        qmf_wrev(db16, N, hp1);
        copy(db16, N, lp2);
        qmf_even(db16, N, hp2);

        return N;
    }

    if (name == "db17"sv) {
        copy_reverse(db17, N, lp1);
        qmf_wrev(db17, N, hp1);
        copy(db17, N, lp2);
        qmf_even(db17, N, hp2);

        return N;
    }

    if (name == "db18"sv) {
        copy_reverse(db18, N, lp1);
        qmf_wrev(db18, N, hp1);
        copy(db18, N, lp2);
        qmf_even(db18, N, hp2);

        return N;
    }

    if (name == "db19"sv) {
        copy_reverse(db19, N, lp1);
        qmf_wrev(db19, N, hp1);
        copy(db19, N, lp2);
        qmf_even(db19, N, hp2);

        return N;
    }

    if (name == "db20"sv) {
        copy_reverse(db20, N, lp1);
        qmf_wrev(db20, N, hp1);
        copy(db20, N, lp2);
        qmf_even(db20, N, hp2);

        return N;
    }

    if (name == "db21"sv) {
        copy_reverse(db21, N, lp1);
        qmf_wrev(db21, N, hp1);
        copy(db21, N, lp2);
        qmf_even(db21, N, hp2);

        return N;
    }

    if (name == "db22"sv) {
        copy_reverse(db22, N, lp1);
        qmf_wrev(db22, N, hp1);
        copy(db22, N, lp2);
        qmf_even(db22, N, hp2);

        return N;
    }

    if (name == "db23"sv) {
        copy_reverse(db23, N, lp1);
        qmf_wrev(db23, N, hp1);
        copy(db23, N, lp2);
        qmf_even(db23, N, hp2);

        return N;
    }

    if (name == "db24"sv) {
        copy_reverse(db24, N, lp1);
        qmf_wrev(db24, N, hp1);
        copy(db24, N, lp2);
        qmf_even(db24, N, hp2);

        return N;
    }

    if (name == "db25"sv) {
        copy_reverse(db25, N, lp1);
        qmf_wrev(db25, N, hp1);
        copy(db25, N, lp2);
        qmf_even(db25, N, hp2);

        return N;
    }

    if (name == "db26"sv) {
        copy_reverse(db26, N, lp1);
        qmf_wrev(db26, N, hp1);
        copy(db26, N, lp2);
        qmf_even(db26, N, hp2);
        return N;
    }

    if (name == "db27"sv) {
        copy_reverse(db27, N, lp1);
        qmf_wrev(db27, N, hp1);
        copy(db27, N, lp2);
        qmf_even(db27, N, hp2);

        return N;
    }

    if (name == "db28"sv) {
        copy_reverse(db28, N, lp1);
        qmf_wrev(db28, N, hp1);
        copy(db28, N, lp2);
        qmf_even(db28, N, hp2);

        return N;
    }

    if (name == "db29"sv) {
        copy_reverse(db29, N, lp1);
        qmf_wrev(db29, N, hp1);
        copy(db29, N, lp2);
        qmf_even(db29, N, hp2);

        return N;
    }

    if (name == "db30"sv) {
        copy_reverse(db30, N, lp1);
        qmf_wrev(db30, N, hp1);
        copy(db30, N, lp2);
        qmf_even(db30, N, hp2);

        return N;
    }

    if (name == "db31"sv) {
        copy_reverse(db31, N, lp1);
        qmf_wrev(db31, N, hp1);
        copy(db31, N, lp2);
        qmf_even(db31, N, hp2);

        return N;
    }

    if (name == "db32"sv) {
        copy_reverse(db32, N, lp1);
        qmf_wrev(db32, N, hp1);
        copy(db32, N, lp2);
        qmf_even(db32, N, hp2);

        return N;
    }

    if (name == "db33"sv) {
        copy_reverse(db33, N, lp1);
        qmf_wrev(db33, N, hp1);
        copy(db33, N, lp2);
        qmf_even(db33, N, hp2);

        return N;
    }

    if (name == "db34"sv) {
        copy_reverse(db34, N, lp1);
        qmf_wrev(db34, N, hp1);
        copy(db34, N, lp2);
        qmf_even(db34, N, hp2);

        return N;
    }

    if (name == "db35"sv) {
        copy_reverse(db35, N, lp1);
        qmf_wrev(db35, N, hp1);
        copy(db35, N, lp2);
        qmf_even(db35, N, hp2);

        return N;
    }

    if (name == "db36"sv) {
        copy_reverse(db36, N, lp1);
        qmf_wrev(db36, N, hp1);
        copy(db36, N, lp2);
        qmf_even(db36, N, hp2);

        return N;
    }

    if (name == "db37"sv) {
        copy_reverse(db37, N, lp1);
        qmf_wrev(db37, N, hp1);
        copy(db37, N, lp2);
        qmf_even(db37, N, hp2);

        return N;
    }

    if (name == "db38"sv) {
        copy_reverse(db38, N, lp1);
        qmf_wrev(db38, N, hp1);
        copy(db38, N, lp2);
        qmf_even(db38, N, hp2);

        return N;
    }

    if (name == "bior1.1"sv) {
        copy_reverse(hm1_11, N, lp1);
        qmf_wrev(h1 + 4, N, hp1);
        copy(h1 + 4, N, lp2);
        qmf_even(hm1_11, N, hp2);
        return N;
    }

    if (name == "bior1.3"sv) {
        copy_reverse(hm1_13, N, lp1);
        qmf_wrev(h1 + 2, N, hp1);
        copy(h1 + 2, N, lp2);
        qmf_even(hm1_13, N, hp2);
        return N;
    }

    if (name == "bior1.5"sv) {
        copy_reverse(hm1_15, N, lp1);
        qmf_wrev(h1, N, hp1);
        copy(h1, N, lp2);
        qmf_even(hm1_15, N, hp2);
        return N;
    }

    if (name == "bior2.2"sv) {
        copy_reverse(hm2_22, N, lp1);
        qmf_wrev(h2 + 6, N, hp1);
        copy(h2 + 6, N, lp2);
        qmf_even(hm2_22, N, hp2);
        return N;
    }

    if (name == "bior2.4"sv) {
        copy_reverse(hm2_24, N, lp1);
        qmf_wrev(h2 + 4, N, hp1);
        copy(h2 + 4, N, lp2);
        qmf_even(hm2_24, N, hp2);
        return N;
    }

    if (name == "bior2.6"sv) {
        copy_reverse(hm2_26, N, lp1);
        qmf_wrev(h2 + 2, N, hp1);
        copy(h2 + 2, N, lp2);
        qmf_even(hm2_26, N, hp2);
        return N;
    }

    if (name == "bior2.8"sv) {
        copy_reverse(hm2_28, N, lp1);
        qmf_wrev(h2, N, hp1);
        copy(h2, N, lp2);
        qmf_even(hm2_28, N, hp2);
        return N;
    }

    if (name == "bior3.1"sv) {
        copy_reverse(hm3_31, N, lp1);
        qmf_wrev(h3 + 8, N, hp1);
        copy(h3 + 8, N, lp2);
        qmf_even(hm3_31, N, hp2);
        return N;
    }

    if (name == "bior3.3"sv) {
        copy_reverse(hm3_33, N, lp1);
        qmf_wrev(h3 + 6, N, hp1);
        copy(h3 + 6, N, lp2);
        qmf_even(hm3_33, N, hp2);
        return N;
    }

    if (name == "bior3.5"sv) {
        copy_reverse(hm3_35, N, lp1);
        qmf_wrev(h3 + 4, N, hp1);
        copy(h3 + 4, N, lp2);
        qmf_even(hm3_35, N, hp2);
        return N;
    }

    if (name == "bior3.7"sv) {
        copy_reverse(hm3_37, N, lp1);
        qmf_wrev(h3 + 2, N, hp1);
        copy(h3 + 2, N, lp2);
        qmf_even(hm3_37, N, hp2);
        return N;
    }

    if (name == "bior3.9"sv) {
        copy_reverse(hm3_39, N, lp1);
        qmf_wrev(h3, N, hp1);
        copy(h3, N, lp2);
        qmf_even(hm3_39, N, hp2);
        return N;
    }

    if (name == "bior4.4"sv) {
        copy_reverse(hm4_44, N, lp1);
        qmf_wrev(h4, N, hp1);
        copy(h4, N, lp2);
        qmf_even(hm4_44, N, hp2);
        return N;
    }

    if (name == "bior5.5"sv) {
        copy_reverse(hm5_55, N, lp1);
        qmf_wrev(h5, N, hp1);
        copy(h5, N, lp2);
        qmf_even(hm5_55, N, hp2);
        return N;
    }

    if (name == "bior6.8"sv) {
        copy_reverse(hm6_68, N, lp1);
        qmf_wrev(h6, N, hp1);
        copy(h6, N, lp2);
        qmf_even(hm6_68, N, hp2);
        return N;
    }

    if (name == "rbior1.1"sv) {
        copy_reverse(h1 + 4, N, lp1);
        qmf_wrev(hm1_11, N, hp1);
        copy(hm1_11, N, lp2);
        qmf_even(h1 + 4, N, hp2);
        return N;
    }

    if (name == "rbior1.3"sv) {
        copy_reverse(h1 + 2, N, lp1);
        qmf_wrev(hm1_13, N, hp1);
        copy(hm1_13, N, lp2);
        qmf_even(h1 + 2, N, hp2);
        return N;
    }

    if (name == "rbior1.5"sv) {
        copy_reverse(h1, N, lp1);
        qmf_wrev(hm1_15, N, hp1);
        copy(hm1_15, N, lp2);
        qmf_even(h1, N, hp2);
        return N;
    }

    if (name == "rbior2.2"sv) {
        copy_reverse(h2 + 6, N, lp1);
        qmf_wrev(hm2_22, N, hp1);
        copy(hm2_22, N, lp2);
        qmf_even(h2 + 6, N, hp2);
        return N;
    }

    if (name == "rbior2.4"sv) {
        copy_reverse(h2 + 4, N, lp1);
        qmf_wrev(hm2_24, N, hp1);
        copy(hm2_24, N, lp2);
        qmf_even(h2 + 4, N, hp2);
        return N;
    }

    if (name == "rbior2.6"sv) {
        copy_reverse(h2 + 2, N, lp1);
        qmf_wrev(hm2_26, N, hp1);
        copy(hm2_26, N, lp2);
        qmf_even(h2 + 2, N, hp2);
        return N;
    }

    if (name == "rbior2.8"sv) {
        copy_reverse(h2, N, lp1);
        qmf_wrev(hm2_28, N, hp1);
        copy(hm2_28, N, lp2);
        qmf_even(h2, N, hp2);
        return N;
    }

    if (name == "rbior3.1"sv) {
        copy_reverse(h3 + 8, N, lp1);
        qmf_wrev(hm3_31, N, hp1);
        copy(hm3_31, N, lp2);
        qmf_even(h3 + 8, N, hp2);
        return N;
    }

    if (name == "rbior3.3"sv) {
        copy_reverse(h3 + 6, N, lp1);
        qmf_wrev(hm3_33, N, hp1);
        copy(hm3_33, N, lp2);
        qmf_even(h3 + 6, N, hp2);
        return N;
    }

    if (name == "rbior3.5"sv) {
        copy_reverse(h3 + 4, N, lp1);
        qmf_wrev(hm3_35, N, hp1);
        copy(hm3_35, N, lp2);
        qmf_even(h3 + 4, N, hp2);
        return N;
    }

    if (name == "rbior3.7"sv) {
        copy_reverse(h3 + 2, N, lp1);
        qmf_wrev(hm3_37, N, hp1);
        copy(hm3_37, N, lp2);
        qmf_even(h3 + 2, N, hp2);
        return N;
    }

    if (name == "rbior3.9"sv) {
        copy_reverse(h3, N, lp1);
        qmf_wrev(hm3_39, N, hp1);
        copy(hm3_39, N, lp2);
        qmf_even(h3, N, hp2);
        return N;
    }

    if (name == "rbior4.4"sv) {
        copy_reverse(h4, N, lp1);
        qmf_wrev(hm4_44, N, hp1);
        copy(hm4_44, N, lp2);
        qmf_even(h4, N, hp2);
        return N;
    }

    if (name == "rbior5.5"sv) {
        copy_reverse(h5, N, lp1);
        qmf_wrev(hm5_55, N, hp1);
        copy(hm5_55, N, lp2);
        qmf_even(h5, N, hp2);
        return N;
    }

    if (name == "rbior6.8"sv) {
        copy_reverse(h6, N, lp1);
        qmf_wrev(hm6_68, N, hp1);
        copy(hm6_68, N, lp2);
        qmf_even(h6, N, hp2);
        return N;
    }

    if (name == "coif1"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif1, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif2"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif2, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif3"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif3, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif4"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif4, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif5"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif5, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif6"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif6, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif7"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif7, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif8"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif8, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif9"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif9, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }

    if (name == "coif10"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif10, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "coif11"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif11, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "coif12"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif12, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "coif13"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif13, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "coif14"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif14, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "coif15"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif15, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "coif16"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif16, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "coif17"sv) {
        auto coeffTemp = std::make_unique<double[]>(N);

        copy(coif17, N, coeffTemp.get());
        for (i = 0; i < N; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copy_reverse(coeffTemp.get(), N, lp1);
        qmf_wrev(coeffTemp.get(), N, hp1);
        copy(coeffTemp.get(), N, lp2);
        qmf_even(coeffTemp.get(), N, hp2);

        return N;
    }
    if (name == "sym2"sv) {
        copy_reverse(sym2, N, lp1);
        qmf_wrev(sym2, N, hp1);
        copy(sym2, N, lp2);
        qmf_even(sym2, N, hp2);
        return N;
    }

    if (name == "sym3"sv) {
        copy_reverse(sym3, N, lp1);
        qmf_wrev(sym3, N, hp1);
        copy(sym3, N, lp2);
        qmf_even(sym3, N, hp2);
        return N;
    }

    if (name == "sym4"sv) {
        copy_reverse(sym4, N, lp1);
        qmf_wrev(sym4, N, hp1);
        copy(sym4, N, lp2);
        qmf_even(sym4, N, hp2);
        return N;
    }

    if (name == "sym5"sv) {
        copy_reverse(sym5, N, lp1);
        qmf_wrev(sym5, N, hp1);
        copy(sym5, N, lp2);
        qmf_even(sym5, N, hp2);
        return N;
    }

    if (name == "sym6"sv) {
        copy_reverse(sym6, N, lp1);
        qmf_wrev(sym6, N, hp1);
        copy(sym6, N, lp2);
        qmf_even(sym6, N, hp2);
        return N;
    }

    if (name == "sym7"sv) {
        copy_reverse(sym7, N, lp1);
        qmf_wrev(sym7, N, hp1);
        copy(sym7, N, lp2);
        qmf_even(sym7, N, hp2);
        return N;
    }

    if (name == "sym8"sv) {
        copy_reverse(sym8, N, lp1);
        qmf_wrev(sym8, N, hp1);
        copy(sym8, N, lp2);
        qmf_even(sym8, N, hp2);
        return N;
    }

    if (name == "sym9"sv) {
        copy_reverse(sym9, N, lp1);
        qmf_wrev(sym9, N, hp1);
        copy(sym9, N, lp2);
        qmf_even(sym9, N, hp2);
        return N;
    }

    if (name == "sym10"sv) {
        copy_reverse(sym10, N, lp1);
        qmf_wrev(sym10, N, hp1);
        copy(sym10, N, lp2);
        qmf_even(sym10, N, hp2);
        return N;
    }
    if (name == "sym11"sv) {
        copy_reverse(sym11, N, lp1);
        qmf_wrev(sym11, N, hp1);
        copy(sym11, N, lp2);
        qmf_even(sym11, N, hp2);
        return N;
    }
    if (name == "sym12"sv) {
        copy_reverse(sym12, N, lp1);
        qmf_wrev(sym12, N, hp1);
        copy(sym12, N, lp2);
        qmf_even(sym12, N, hp2);
        return N;
    }
    if (name == "sym13"sv) {
        copy_reverse(sym13, N, lp1);
        qmf_wrev(sym13, N, hp1);
        copy(sym13, N, lp2);
        qmf_even(sym13, N, hp2);
        return N;
    }
    if (name == "sym14"sv) {
        copy_reverse(sym14, N, lp1);
        qmf_wrev(sym14, N, hp1);
        copy(sym14, N, lp2);
        qmf_even(sym14, N, hp2);
        return N;
    }
    if (name == "sym15"sv) {
        copy_reverse(sym15, N, lp1);
        qmf_wrev(sym15, N, hp1);
        copy(sym15, N, lp2);
        qmf_even(sym15, N, hp2);
        return N;
    }
    if (name == "sym16"sv) {
        copy_reverse(sym16, N, lp1);
        qmf_wrev(sym16, N, hp1);
        copy(sym16, N, lp2);
        qmf_even(sym16, N, hp2);
        return N;
    }
    if (name == "sym17"sv) {
        copy_reverse(sym17, N, lp1);
        qmf_wrev(sym17, N, hp1);
        copy(sym17, N, lp2);
        qmf_even(sym17, N, hp2);
        return N;
    }
    if (name == "sym18"sv) {
        copy_reverse(sym18, N, lp1);
        qmf_wrev(sym18, N, hp1);
        copy(sym18, N, lp2);
        qmf_even(sym18, N, hp2);
        return N;
    }
    if (name == "sym19"sv) {
        copy_reverse(sym19, N, lp1);
        qmf_wrev(sym19, N, hp1);
        copy(sym19, N, lp2);
        qmf_even(sym19, N, hp2);
        return N;
    }
    if (name == "sym20"sv) {
        copy_reverse(sym20, N, lp1);
        qmf_wrev(sym20, N, hp1);
        copy(sym20, N, lp2);
        qmf_even(sym20, N, hp2);
        return N;
    }
    if (name == "meyer"sv) {
        copy_reverse(meyer, N, lp1);
        qmf_wrev(meyer, N, hp1);
        copy(meyer, N, lp2);
        qmf_even(meyer, N, hp2);
        return N;
    }

    printf("\n Filter Not in Database \n");
    return -1;
}
