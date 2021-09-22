#ifndef WAVELIB_H_
#define WAVELIB_H_

#include "tcb/span.hpp"

#include <algorithm>
#include <memory>
#include <string>

#if defined(_MSC_VER)
#pragma warning(disable : 4200)
#pragma warning(disable : 4996)
#endif

#ifndef fft_type
#define fft_type double
#endif

#ifndef cplx_type
#define cplx_type double
#endif

template <typename T>
auto makeZeros(std::size_t length) -> std::unique_ptr<T[]>
{
    auto ptr = std::make_unique<T[]>(length);
    std::fill(ptr.get(), ptr.get() + length, T {});
    return ptr;
}

enum struct SignalExtension {
    periodic,
    symmetric,
};

[[nodiscard]] inline auto toString(SignalExtension ext) -> std::string
{
    if (ext == SignalExtension::periodic) {
        return "periodic";
    }
    return "symmetric";
}

enum struct ConvolutionMethod {
    direct,
    fft,
};

[[nodiscard]] inline auto toString(ConvolutionMethod method) -> std::string
{
    if (method == ConvolutionMethod::direct) {
        return "direct";
    }
    return "fft";
}

struct CplxData {
    cplx_type re;
    cplx_type im;
};

struct Wavelet {
    explicit Wavelet(char const* wname);

    [[nodiscard]] auto size() const noexcept -> int { return static_cast<int>(size_); }
    [[nodiscard]] auto name() const noexcept -> std::string const& { return name_; }

    [[nodiscard]] auto lpd() const noexcept -> double const* { return lpd_.data(); }
    [[nodiscard]] auto hpd() const noexcept -> double const* { return hpd_.data(); }
    [[nodiscard]] auto lpr() const noexcept -> double const* { return lpr_.data(); }
    [[nodiscard]] auto hpr() const noexcept -> double const* { return hpr_.data(); }

    [[nodiscard]] auto lpdLen() const noexcept -> int { return static_cast<int>(lpd_.size()); }
    [[nodiscard]] auto hpdLen() const noexcept -> int { return static_cast<int>(hpd_.size()); }
    [[nodiscard]] auto lprLen() const noexcept -> int { return static_cast<int>(lpr_.size()); }
    [[nodiscard]] auto hprLen() const noexcept -> int { return static_cast<int>(hpr_.size()); }

private:
    std::string name_;

    // When all filters are of the same length.
    // [Matlab uses zero-padding to make all filters of the same length]
    std::size_t size_;
    std::unique_ptr<double[]> params_;

    lt::span<double> lpd_;
    lt::span<double> hpd_;
    lt::span<double> lpr_;
    lt::span<double> hpr_;
};

struct FftData {
    fft_type re;
    fft_type im;
};

struct FftSet {
    int N;
    int sgn;
    int factors[64];
    int lf;
    int lt;
    std::unique_ptr<FftData[]> data;
};

auto fftInit(int n, int sgn) -> std::unique_ptr<FftSet>;

struct FftRealSet {
    std::unique_ptr<FftSet> cobj;
    std::unique_ptr<FftData[]> data;
};

auto fftRealInit(int n, int sgn) -> std::unique_ptr<FftRealSet>;

struct ConvSet {
    std::unique_ptr<FftRealSet> fobj;
    std::unique_ptr<FftRealSet> iobj;
    int ilen1;
    int ilen2;
    int clen;
};

auto convInit(int n, int l) -> std::unique_ptr<ConvSet>;

struct WaveletTransform {
    WaveletTransform(Wavelet& wave, char const* method, int siglength, int j);

    [[nodiscard]] auto wave() const noexcept -> Wavelet const& { return *wave_; }
    [[nodiscard]] auto levels() const noexcept -> int { return levels_; }
    [[nodiscard]] auto method() const noexcept -> std::string const& { return method_; }

    auto extension(SignalExtension ext) -> void;
    [[nodiscard]] auto extension() const noexcept -> SignalExtension { return ext_; }

    auto convMethod(ConvolutionMethod method) -> void;
    [[nodiscard]] auto convMethod() const noexcept -> ConvolutionMethod { return cmethod_; }

    [[nodiscard]] auto output() const noexcept -> lt::span<double>;
    [[nodiscard]] auto approx() const noexcept -> lt::span<double>;
    [[nodiscard]] auto detail(std::size_t level) const noexcept -> lt::span<double>;

private:
    Wavelet* wave_;
    int levels_;
    std::string method_;
    SignalExtension ext_;
    ConvolutionMethod cmethod_;

    double* output_;

public:
    std::unique_ptr<ConvSet> cobj;
    int siglength; // Length of the original signal.
    int modwtsiglength; // Modified signal length for MODWT
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise

    int N; //
    int cfftset;
    int zpad;
    int length[102];
    std::unique_ptr<double[]> params;
};

struct WtreeSet {
    Wavelet* wave;
    ConvSet* cobj;
    std::string method;
    int siglength; // Length of the original signal.
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise
    std::string ext; // Type of Extension used - "per" or "sym"

    int N; //
    int nodes;
    int cfftset;
    int zpad;
    int length[102];
    double* output;
    int* nodelength;
    int* coeflength;
    std::unique_ptr<double[]> params;
};

auto wtreeInit(Wavelet* wave, int siglength, int j) -> WtreeSet*;

struct WptSet {
    Wavelet* wave;
    ConvSet* cobj;
    int siglength; // Length of the original signal.
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise
    std::string ext; // Type of Extension used - "per" or "sym"
    std::string entropy;
    double eparam;

    int N; //
    int nodes;
    int length[102];
    double* output;
    double* costvalues;
    double* basisvector;
    int* nodeindex;
    int* numnodeslevel;
    int* coeflength;
    std::unique_ptr<double[]> params;
};

auto wptInit(Wavelet* wave, int siglength, int j) -> WptSet*;

struct CwaveletTransform {
    std::string wave; // Wavelet - morl/morlet,paul,dog/dgauss
    int siglength; // Length of Input Data
    int J; // Total Number of Scales
    double s0; // Smallest scale. It depends on the sampling rate. s0 <= 2 * dt for most wavelets
    double dt; // Sampling Rate
    double dj; // Separation between scales. eg., scale = s0 * 2 ^ ( [0:N-1] *dj ) or scale = s0 *[0:N-1] * dj
    std::string type; // Scale Type - Power or Linear
    int pow; // Base of Power in case type = pow. Typical value is pow = 2
    int sflag;
    int pflag;
    int npad;
    int mother;
    double m; // Wavelet parameter param
    double smean; // Input Signal mean

    CplxData* output;
    double* scale;
    double* period;
    double* coi;
    std::unique_ptr<double[]> params;
};

auto cwtInit(char const* wave, double param, int siglength, double dt, int j) -> CwaveletTransform*;

struct Wt2Set {
    Wavelet* wave;
    std::string method;
    int rows; // Matrix Number of rows
    int cols; // Matrix Number of columns
    int outlength; // Length of the output DWT vector
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    std::string ext; // Type of Extension used - "per" or "sym"
    int coeffaccesslength;

    int N; //
    int* dimensions;
    int* coeffaccess;
    std::unique_ptr<int[]> params;
};

auto wt2Init(Wavelet& wave, char const* method, int rows, int cols, int j) -> Wt2Set*;

auto dwt(WaveletTransform& wt, double const* inp) -> void;

auto idwt(WaveletTransform& wt, double* dwtop) -> void;

auto wtree(WtreeSet* wt, double const* inp) -> void;

auto dwpt(WptSet* wt, double const* inp) -> void;

auto idwpt(WptSet* wt, double* dwtop) -> void;

auto swt(WaveletTransform& wt, double const* inp) -> void;

auto iswt(WaveletTransform& wt, double* swtop) -> void;

auto modwt(WaveletTransform& wt, double const* inp) -> void;

auto imodwt(WaveletTransform& wt, double* oup) -> void;

auto setWTREEExtension(WtreeSet* wt, char const* extension) -> void;

auto setDWPTExtension(WptSet* wt, char const* extension) -> void;

auto setDWT2Extension(Wt2Set* wt, char const* extension) -> void;

auto setDWPTEntropy(WptSet* wt, char const* entropy, double eparam) -> void;

auto getWTREENodelength(WtreeSet* wt, int x) -> int;

auto getWTREECoeffs(WtreeSet* wt, int x, int y, double* coeffs, int n) -> void;

auto getDWPTNodelength(WptSet* wt, int x) -> int;

auto setCWTScales(CwaveletTransform* wt, double s0, double dj, char const* type, int power) -> void;

auto cwt(CwaveletTransform* wt, double const* inp) -> void;

auto icwt(CwaveletTransform* wt, double* cwtop) -> void;

auto dwt2(Wt2Set* wt, double* inp) -> std::unique_ptr<double[]>;

auto idwt2(Wt2Set* wt, double* wavecoeff, double* oup) -> void;

auto swt2(Wt2Set* wt, double* inp) -> std::unique_ptr<double[]>;

auto iswt2(Wt2Set* wt, double const* wavecoeffs, double* oup) -> void;

auto modwt2(Wt2Set* wt, double* inp) -> std::unique_ptr<double[]>;

auto imodwt2(Wt2Set* wt, double* wavecoeff, double* oup) -> void;

auto getWT2Coeffs(Wt2Set* wt, double* wcoeffs, int level, char const* type, int* rows, int* cols) -> double*;

auto dispWT2Coeffs(double* a, int row, int col) -> void;

auto waveSummary(Wavelet const& obj) -> void;

auto wtSummary(WaveletTransform const& wt) -> void;

auto wtreeSummary(WtreeSet* wt) -> void;

auto wptSummary(WptSet* wt) -> void;

auto cwtSummary(CwaveletTransform* wt) -> void;

auto wt2Summary(Wt2Set* wt) -> void;

auto wtreeFree(WtreeSet* object) -> void;

auto wptFree(WptSet* object) -> void;

auto cwtFree(CwaveletTransform* object) -> void;

auto wt2Free(Wt2Set* wt) -> void;

#endif /* WAVELIB_H_ */
