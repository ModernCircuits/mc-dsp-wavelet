#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/random.hpp>
#include <mc/testing/test.hpp>

namespace mc {

auto absmax(float* array, std::size_t n) -> float
{
    auto max = 0.0F;
    for (auto i = std::size_t{0}; i < n; ++i) {
        if (std::abs(array[i]) >= max) { max = std::abs(array[i]); }
    }
    return max;
}

auto sum1(float const* array, std::size_t n) -> float
{
    auto sum = 0.0F;
    for (auto i = std::size_t{0}; i < n; ++i) { sum += array[i]; }
    return sum;
}

auto sum2(float const* array, std::size_t n) -> float
{
    auto sum = 0.0F;
    for (std::size_t i = 0; i < n; i += 2) { sum += array[i]; }
    return sum;
}

auto sum3(float const* array, std::size_t n) -> float
{
    auto sum = 0.0F;
    for (std::size_t i = 1; i < n; i += 2) { sum += array[i]; }
    return sum;
}

// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
auto sum4(float const* array, std::size_t n) -> float
{
    auto sum = 0.0F;
    for (std::size_t i = 0; i < n; i += 1) { sum += array[i] * array[i]; }
    return sum;
}

// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
auto sum5(float const* array, std::size_t n, std::size_t m) -> float
{
    auto sum = 0.0F;
    for (std::size_t i = 2 * m; i < n; i += 1) { sum += array[i] * array[i - 2 * m]; }
    return sum;
}

auto rmsError(float const* data, float const* rec, std::size_t n) -> float
{
    float sum = 0;
    for (std::size_t i = 0; i < n; ++i) { sum += (data[i] - rec[i]) * (data[i] - rec[i]); }
    return sqrt(sum / ((float)n - 1));
}

auto relError(float const* data, float const* rec, std::size_t n) -> float
{
    float sum1 = 0;
    float sum2 = 0;
    for (std::size_t i = 0; i < n; ++i) {
        sum1 += (data[i] - rec[i]) * (data[i] - rec[i]);
        sum2 += data[i] * data[i];
    }
    return sqrt(sum1) / sqrt(sum2);
}

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

auto corrcoef(int n, float const* x, float const* y) -> float
{
    float cc   = NAN;
    float xm   = NAN;
    float ym   = NAN;
    float tx   = NAN;
    float ty   = NAN;
    float num  = NAN;
    float den1 = NAN;
    float den2 = NAN;
    int i      = 0;
    xm = ym = 0.0F;
    for (i = 0; i < n; ++i) {
        xm += x[i];
        ym += y[i];
    }

    xm  = xm / n;
    ym  = ym / n;
    num = den1 = den2 = 0.0F;

    for (i = 0; i < n; ++i) {
        tx = x[i] - xm;
        ty = y[i] - ym;
        num += (tx * ty);
        den1 += (tx * tx);
        den2 += (ty * ty);
    }

    cc = num / sqrt(den1 * den2);

    return cc;
}

}  // namespace mc
