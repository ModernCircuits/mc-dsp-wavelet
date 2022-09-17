#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/random.hpp>
#include <mc/testing/test.hpp>

namespace mc {

auto generateRnd() -> float
{
    std::random_device rd{};
    auto gen = std::mt19937{rd()};
    auto dis = std::uniform_real_distribution<float>{1.0F, 100.0F};
    return dis(gen);
}

auto split(String const& s, char delim) -> Vector<String>
{
    auto result = Vector<String>{};
    auto ss     = std::stringstream(s);
    auto item   = String{};

    while (std::getline(ss, item, delim)) { result.push_back(item); }
    return result;
}

auto loadTestData(char const* filePath) -> TestData<float>
{
    auto parseLine = [](auto const& line) {
        auto splits = split(line, ' ');
        auto values = Vector<float>{};
        for (auto const& s : splits) { values.push_back(std::stof(s)); }
        return values;
    };

    auto file   = std::fstream{filePath, std::ios::in};
    auto tmp    = String{};
    auto result = TestData<float>{};

    if (file.is_open()) {
        while (std::getline(file, tmp)) { result.push_back(parseLine(tmp)); }
        file.close();
    }

    return result;
}

auto toFloat(TestData<float> const& d) -> TestData<float>
{
    auto result = TestData<float>{};
    for (auto const& line : d) {
        result.emplace_back();
        auto& lineFloat = result.back();
        for (auto value : line) { lineFloat.push_back(static_cast<float>(value)); }
    }

    return result;
}

auto generateRandomTestData(std::size_t n) -> Vector<float>
{
    Vector<float> data(n);
    auto rd  = std::random_device{};
    auto gen = std::mt19937{rd()};
    auto dis = std::uniform_real_distribution<float>{-1.0F, 1.0F};
    std::generate(data.begin(), data.end(), [&]() { return dis(gen); });
    return data;
}

}  // namespace mc
