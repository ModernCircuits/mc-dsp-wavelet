#ifndef WAVELIB_H_
#define WAVELIB_H_

#include "tcb/span.hpp"

#include "ConvolutionMethod.hpp"
#include "SignalExtension.hpp"
#include "Wavelet.hpp"

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

struct CplxData {
    cplx_type re;
    cplx_type im;
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

struct Convolution {
    std::unique_ptr<FftRealSet> fobj;
    std::unique_ptr<FftRealSet> iobj;
    int ilen1;
    int ilen2;
    int clen;
};

auto convInit(int n, int l) -> std::unique_ptr<Convolution>;

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
    std::unique_ptr<Convolution> cobj;
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

struct WaveletTree {
    Wavelet* wave;
    Convolution* cobj;
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

auto wtreeInit(Wavelet* wave, int siglength, int j) -> WaveletTree*;

struct WaveletPacketTransform {
    Wavelet* wave;
    Convolution* cobj;
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

auto wptInit(Wavelet* wave, int siglength, int j) -> WaveletPacketTransform*;

struct ComplexWaveletTransform {
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

auto cwtInit(char const* wave, double param, int siglength, double dt, int j) -> ComplexWaveletTransform*;

struct WaveletTransform2D {
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

auto wt2Init(Wavelet& wave, char const* method, int rows, int cols, int j) -> WaveletTransform2D*;

auto dwt(WaveletTransform& wt, double const* inp) -> void;

auto idwt(WaveletTransform& wt, double* dwtop) -> void;

auto wtree(WaveletTree* wt, double const* inp) -> void;

auto dwpt(WaveletPacketTransform* wt, double const* inp) -> void;

auto idwpt(WaveletPacketTransform* wt, double* dwtop) -> void;

auto swt(WaveletTransform& wt, double const* inp) -> void;

auto iswt(WaveletTransform& wt, double* swtop) -> void;

auto modwt(WaveletTransform& wt, double const* inp) -> void;

auto imodwt(WaveletTransform& wt, double* oup) -> void;

auto setWTREEExtension(WaveletTree* wt, char const* extension) -> void;

auto setDWPTExtension(WaveletPacketTransform* wt, char const* extension) -> void;

auto setDWT2Extension(WaveletTransform2D* wt, char const* extension) -> void;

auto setDWPTEntropy(WaveletPacketTransform* wt, char const* entropy, double eparam) -> void;

auto getWTREENodelength(WaveletTree* wt, int x) -> int;

auto getWTREECoeffs(WaveletTree* wt, int x, int y, double* coeffs, int n) -> void;

auto getDWPTNodelength(WaveletPacketTransform* wt, int x) -> int;

auto setCWTScales(ComplexWaveletTransform* wt, double s0, double dj, char const* type, int power) -> void;

auto cwt(ComplexWaveletTransform* wt, double const* inp) -> void;

auto icwt(ComplexWaveletTransform* wt, double* cwtop) -> void;

auto dwt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>;

auto idwt2(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void;

auto swt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>;

auto iswt2(WaveletTransform2D* wt, double const* wavecoeffs, double* oup) -> void;

auto modwt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>;

auto imodwt2(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void;

auto getWT2Coeffs(WaveletTransform2D* wt, double* wcoeffs, int level, char const* type, int* rows, int* cols) -> double*;

auto dispWT2Coeffs(double* a, int row, int col) -> void;

auto waveSummary(Wavelet const& obj) -> void;

auto wtSummary(WaveletTransform const& wt) -> void;

auto wtreeSummary(WaveletTree* wt) -> void;

auto wptSummary(WaveletPacketTransform* wt) -> void;

auto cwtSummary(ComplexWaveletTransform* wt) -> void;

auto wt2Summary(WaveletTransform2D* wt) -> void;

auto wtreeFree(WaveletTree* object) -> void;

auto wptFree(WaveletPacketTransform* object) -> void;

auto cwtFree(ComplexWaveletTransform* object) -> void;

auto wt2Free(WaveletTransform2D* wt) -> void;

#endif /* WAVELIB_H_ */
