#pragma once

#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <vector>

inline auto readFileToVector(char const* filePath) -> std::vector<double>
{
    auto* ifp = std::fopen(filePath, "r");
    if (ifp == nullptr) {
        std::printf("Cannot Open File: %s\n", filePath);
        std::exit(100);
    }

    double temp[1200] {};
    auto i = std::size_t { 0 };
    while (std::feof(ifp) == 0) {
        std::fscanf(ifp, "%lf \n", &temp[i]);
        i++;
    }

    std::fclose(ifp);
    return std::vector<double>(std::begin(temp), std::begin(temp) + i);
}