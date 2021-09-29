#include "lt/cmath.hpp"
#include "lt/random.hpp"

#include "lt/testing/test.hpp"

#include "lt/cstdlib.hpp"

auto absmax(float* array, std::size_t n) -> float
{
    auto max = 0.0;
    for (auto i = std::size_t { 0 }; i < n; ++i) {
        if (fabs(array[i]) >= max) {
            max = fabs(array[i]);
        }
    }
    return max;
}

auto sum1(float const* array, std::size_t n) -> float
{
    auto sum = 0.0;
    for (auto i = std::size_t { 0 }; i < n; ++i) {
        sum += array[i];
    }
    return sum;
}

auto sum2(float const* array, std::size_t n) -> float
{
    auto sum = 0.0;
    for (std::size_t i = 0; i < n; i += 2) {
        sum += array[i];
    }
    return sum;
}
auto sum3(float const* array, std::size_t n) -> float
{
    auto sum = 0.0;
    for (std::size_t i = 1; i < n; i += 2) {
        sum += array[i];
    }
    return sum;
}
// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
auto sum4(float const* array, std::size_t n) -> float
{
    auto sum = 0.0;
    for (std::size_t i = 0; i < n; i += 1) {
        sum += array[i] * array[i];
    }
    return sum;
}
// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
auto sum5(float const* array, std::size_t n, std::size_t m) -> float
{
    auto sum = 0.0;
    for (std::size_t i = 2 * m; i < n; i += 1) {
        sum += array[i] * array[i - 2 * m];
    }
    return sum;
}

auto rmsError(float const* data, float const* rec, std::size_t n) -> float
{
    float sum = 0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += (data[i] - rec[i]) * (data[i] - rec[i]);
    }
    return std::sqrt(sum / ((float)n - 1));
}

auto relError(float const* data, float const* rec, std::size_t n) -> float
{
    float sum1 = 0;
    float sum2 = 0;
    for (std::size_t i = 0; i < n; ++i) {
        sum1 += (data[i] - rec[i]) * (data[i] - rec[i]);
        sum2 += data[i] * data[i];
    }
    return std::sqrt(sum1) / std::sqrt(sum2);
}

auto generateRnd() -> float
{
    std::random_device rd {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<float> { 1.0, 100.0 };
    return dis(gen);
}

auto split(std::string const& s, char delim) -> std::vector<std::string>
{
    auto result = std::vector<std::string> {};
    auto ss = std::stringstream(s);
    auto item = std::string {};

    while (std::getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}

auto loadTestData(char const* filePath) -> TestData<float>
{
    auto parseLine = [](auto const& line) {
        auto splits = split(line, ' ');
        auto values = std::vector<float> {};
        for (auto const& s : splits) {
            values.push_back(std::stod(s));
        }
        return values;
    };

    auto file = std::fstream { filePath, std::ios::in };
    auto tmp = std::string {};
    auto result = TestData<float> {};

    if (file.is_open()) {
        while (std::getline(file, tmp)) {
            result.push_back(parseLine(tmp));
        }
        file.close();
    }

    return result;
}

auto toFloat(TestData<float> const& d) -> TestData<float>
{
    auto result = TestData<float> {};
    for (auto const& line : d) {
        result.emplace_back();
        auto& lineFloat = result.back();
        for (auto value : line) {
            lineFloat.push_back(static_cast<float>(value));
        }
    }

    return result;
}