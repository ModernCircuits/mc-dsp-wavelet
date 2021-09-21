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

struct cplx_data {
    cplx_type re;
    cplx_type im;
};

struct wavelet {
    explicit wavelet(char const* wname);

    [[nodiscard]] auto size() const noexcept -> int { return static_cast<int>(size_); }
    [[nodiscard]] auto name() const noexcept -> std::string const& { return name_; }

    [[nodiscard]] auto lpd() const noexcept -> double const* { return lpd_.data(); }
    [[nodiscard]] auto hpd() const noexcept -> double const* { return hpd_.data(); }
    [[nodiscard]] auto lpr() const noexcept -> double const* { return lpr_.data(); }
    [[nodiscard]] auto hpr() const noexcept -> double const* { return hpr_.data(); }

    [[nodiscard]] auto lpd_len() const noexcept -> int { return static_cast<int>(lpd_.size()); }
    [[nodiscard]] auto hpd_len() const noexcept -> int { return static_cast<int>(hpd_.size()); }
    [[nodiscard]] auto lpr_len() const noexcept -> int { return static_cast<int>(lpr_.size()); }
    [[nodiscard]] auto hpr_len() const noexcept -> int { return static_cast<int>(hpr_.size()); }

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

struct fft_data {
    fft_type re;
    fft_type im;
};

struct fft_set {
    int N;
    int sgn;
    int factors[64];
    int lf;
    int lt;
    std::unique_ptr<fft_data[]> data;
};

auto fft_init(int N, int sgn) -> std::unique_ptr<fft_set>;

struct fft_real_set {
    std::unique_ptr<fft_set> cobj;
    std::unique_ptr<fft_data[]> data;
};

auto fft_real_init(int N, int sgn) -> std::unique_ptr<fft_real_set>;

struct conv_set {
    std::unique_ptr<fft_real_set> fobj;
    std::unique_ptr<fft_real_set> iobj;
    int ilen1;
    int ilen2;
    int clen;
};

auto conv_init(int N, int L) -> std::unique_ptr<conv_set>;

struct wt_set {
    wavelet* wave;
    std::unique_ptr<conv_set> cobj;
    std::string method;
    int siglength; // Length of the original signal.
    int modwtsiglength; // Modified signal length for MODWT
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise
    std::string ext; // Type of Extension used - "per" or "sym"
    std::string cmethod; // Convolution Method - "direct" or "FFT"

    int N; //
    int cfftset;
    int zpad;
    int length[102];
    double* output;
    std::unique_ptr<double[]> params;
};

auto wt_init(wavelet& wave, char const* method, int siglength, int J) -> wt_set*;

struct wtree_set {
    wavelet* wave;
    conv_set* cobj;
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

auto wtree_init(wavelet* wave, int siglength, int J) -> wtree_set*;

struct wpt_set {
    wavelet* wave;
    conv_set* cobj;
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

auto wpt_init(wavelet* wave, int siglength, int J) -> wpt_set*;

struct cwt_set {
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

    cplx_data* output;
    double* scale;
    double* period;
    double* coi;
    std::unique_ptr<double[]> params;
};

auto cwt_init(char const* wave, double param, int siglength, double dt, int J) -> cwt_set*;

struct wt2_set {
    wavelet* wave;
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

auto wt2_init(wavelet& wave, char const* method, int rows, int cols, int J) -> wt2_set*;

void dwt(wt_set* wt, double const* inp);

void idwt(wt_set* wt, double* dwtop);

void wtree(wtree_set* wt, double const* inp);

void dwpt(wpt_set* wt, double const* inp);

void idwpt(wpt_set* wt, double* dwtop);

void swt(wt_set* wt, double const* inp);

void iswt(wt_set* wt, double* swtop);

void modwt(wt_set* wt, double const* inp);

void imodwt(wt_set* wt, double* oup);

void setDWTExtension(wt_set* wt, char const* extension);

void setWTREEExtension(wtree_set* wt, char const* extension);

void setDWPTExtension(wpt_set* wt, char const* extension);

void setDWT2Extension(wt2_set* wt, char const* extension);

void setDWPTEntropy(wpt_set* wt, char const* entropy, double eparam);

void setWTConv(wt_set* wt, char const* cmethod);

auto getWTREENodelength(wtree_set* wt, int X) -> int;

void getWTREECoeffs(wtree_set* wt, int X, int Y, double* coeffs, int N);

auto getDWPTNodelength(wpt_set* wt, int X) -> int;

void setCWTScales(cwt_set* wt, double s0, double dj, char const* type, int power);

void cwt(cwt_set* wt, double const* inp);

void icwt(cwt_set* wt, double* cwtop);

auto dwt2(wt2_set* wt, double* inp) -> std::unique_ptr<double[]>;

void idwt2(wt2_set* wt, double* wavecoeff, double* oup);

auto swt2(wt2_set* wt, double* inp) -> std::unique_ptr<double[]>;

void iswt2(wt2_set* wt, double const* wavecoeffs, double* oup);

auto modwt2(wt2_set* wt, double* inp) -> std::unique_ptr<double[]>;

void imodwt2(wt2_set* wt, double* wavecoeff, double* oup);

auto getWT2Coeffs(wt2_set* wt, double* wcoeffs, int level, char const* type, int* rows, int* cols) -> double*;

void dispWT2Coeffs(double* A, int row, int col);

void wave_summary(wavelet const& obj);

void wt_summary(wt_set* wt);

void wtree_summary(wtree_set* wt);

void wpt_summary(wpt_set* wt);

void cwt_summary(cwt_set* wt);

void wt2_summary(wt2_set* wt);

void wt_free(wt_set* object);

void wtree_free(wtree_set* object);

void wpt_free(wpt_set* object);

void cwt_free(cwt_set* object);

void wt2_free(wt2_set* wt);

#endif /* WAVELIB_H_ */
