#include "OverlapSaveConvolver.hpp"

#define REAL 0
#define IMAG 1

// Raises an exception if complex_size!=(real_size/2+1), being "/" an integer division.
auto checkRealComplexRatio(std::size_t const realSize, std::size_t const complexSize, std::string const& funcName) -> void
{
    if (complexSize != (realSize / 2 + 1)) {
        throw std::runtime_error(std::string("[ERROR] ") + funcName
            + ": size of ComplexSignal must equal size(DoubleSignal)/2+1. "
            + " Sizes were (double, complex): " + iterableToString({ realSize, complexSize }));
    }
}

// Raises an exception with the given message if a>b.
auto checkALessEqualB(std::size_t const a, std::size_t const b, std::string const& message) -> void
{
    checkTwoElements(
        a, b, [](auto ax, auto bx) { return ax > bx; }, message);
}

/// HELPERS

auto pow2Ceil(std::size_t x) -> size_t { return std::pow(2, std::ceil(std::log2(x))); }

// the basic constructor allocates an aligned, double array, which is zeroed by the superclass
DoubleSignal::DoubleSignal(std::size_t size)
    : Signal(fftw_alloc_real(size), size)
{
}

DoubleSignal::DoubleSignal(double* data, size_t size)
    : DoubleSignal(size)
{
    std::memcpy(data_, data, sizeof(double) * size);
}
DoubleSignal::DoubleSignal(double* data, size_t size, size_t padBef, size_t padAft)
    : DoubleSignal(size + padBef + padAft)
{
    std::memcpy(data_ + padBef, data, sizeof(double) * size);
}
// the destructor frees the only resource allocated
DoubleSignal::~DoubleSignal() { fftw_free(data_); }
auto DoubleSignal::operator+=(double const x) -> void
{
    for (std::size_t i = 0; i < size_; ++i) {
        data_[i] += x;
    }
}
auto DoubleSignal::operator*=(double const x) -> void
{
    for (std::size_t i = 0; i < size_; ++i) {
        data_[i] *= x;
    }
}
auto DoubleSignal::operator/=(double const x) -> void
{
    for (std::size_t i = 0; i < size_; ++i) {
        data_[i] /= x;
    }
}

ComplexSignal::ComplexSignal(std::size_t size)
    : Signal(fftw_alloc_complex(size), size)
{
}

ComplexSignal::~ComplexSignal() { fftw_free(data_); }

auto ComplexSignal::operator*=(double const x) -> void
{
    for (std::size_t i = 0; i < size_; ++i) {
        data_[i][REAL] *= x;
        data_[i][IMAG] *= x;
    }
}

auto ComplexSignal::operator+=(double const x) -> void
{
    for (std::size_t i = 0; i < size_; ++i) {
        data_[i][REAL] += x;
    }
}

auto ComplexSignal::operator+=(const fftw_complex x) -> void
{
    for (std::size_t i = 0; i < size_; ++i) {
        data_[i][REAL] += x[REAL];
        data_[i][IMAG] += x[IMAG];
    }
}

auto ComplexSignal::print(std::string const& /*name*/) -> void
{
    for (std::size_t i = 0; i < size_; ++i) {
        // printf("%s[%zu]\t=\t(%f, i%f)\n", name.c_str(), i, data_[i][REAL], data_[i][IMAG]);
    }
}

auto spectralConvolution(ComplexSignal const& a, ComplexSignal const& b, ComplexSignal& result) -> void
{
    std::size_t const kSizeA = a.size();
    std::size_t const kSizeB = b.size();
    std::size_t const kSizeResult = result.size();
    checkAllEqual({ kSizeA, kSizeB, kSizeResult }, "SpectralConvolution: all sizes must be equal and are");
    for (std::size_t i = 0; i < kSizeA; ++i) {
        // a+ib * c+id = ac+iad+ibc-bd = ac-bd + i(ad+bc)
        result[i][REAL] = a[i][REAL] * b[i][REAL] - a[i][IMAG] * b[i][IMAG];
        result[i][IMAG] = a[i][IMAG] * b[i][REAL] + a[i][REAL] * b[i][IMAG];
    }
}

auto spectralCorrelation(ComplexSignal const& a, ComplexSignal const& b, ComplexSignal& result) -> void
{
    std::size_t const kSizeA = a.size();
    std::size_t const kSizeB = b.size();
    std::size_t const kSizeResult = result.size();
    checkAllEqual({ kSizeA, kSizeB, kSizeResult }, "SpectralCorrelation: all sizes must be equal and are");
    for (std::size_t i = 0; i < kSizeA; ++i) {
        // a * conj(b) = a+ib * c-id = ac-iad+ibc+bd = ac+bd + i(bc-ad)
        result[i][REAL] = a[i][REAL] * b[i][REAL] + a[i][IMAG] * b[i][IMAG];
        result[i][IMAG] = a[i][IMAG] * b[i][REAL] - a[i][REAL] * b[i][IMAG];
    }
}

FftForwardPlan::FftForwardPlan(DoubleSignal& fs, ComplexSignal& cs)
    : FftPlan(fftw_plan_dft_r2c_1d(fs.size(), fs.data(), cs.data(), FFTW_ESTIMATE))
{
    checkRealComplexRatio(fs.size(), cs.size(), "FftForwardPlan");
}

FftBackwardPlan::FftBackwardPlan(ComplexSignal& cs, DoubleSignal& fs)
    : FftPlan(fftw_plan_dft_c2r_1d(fs.size(), cs.data(), fs.data(), FFTW_ESTIMATE))
{
    checkRealComplexRatio(fs.size(), cs.size(), "FftBackwardPlan");
}

OverlapSaveConvolver::OverlapSaveConvolver(DoubleSignal& signal, DoubleSignal& patch, std::string const& /*wisdomPath*/)
    : signalSize_(signal.size())
    , patchSize_(patch.size())
    , resultSize_(signalSize_ + patchSize_ - 1)
    , paddedPatch_(patch.data(), patchSize_, 0, 2 * pow2Ceil(patchSize_) - patchSize_)
    , resultChunksize_(paddedPatch_.size())
    , resultChunksizeComplex_(resultChunksize_ / 2 + 1)
    , result_stride_(resultChunksize_ - patchSize_ + 1)
    , paddedPatchComplex_(resultChunksizeComplex_)
    , paddedSignal_(signal.data(), signalSize_, patchSize_ - 1, resultChunksize_ - (resultSize_ % result_stride_))
    , state_(State::Uninitialized)
{
    // end of initializer list, now check that len(signal)>=len(patch)
    checkALessEqualB(patchSize_, signalSize_, "OverlapSaveConvolver: len(signal) can't be smaller than len(patch)!");
    // and load the wisdom if required. If unsuccessful, no exception thrown, just print a warning.

    // chunk the signal into strides of same size as padded patch
    // and make complex counterparts too, as well as the corresponding xcorr signals
    for (std::size_t i = 0; i <= paddedSignal_.size() - resultChunksize_; i += result_stride_) {
        inputChunks_.push_back(std::make_unique<DoubleSignal>(&paddedSignal_[i], resultChunksize_));
        inputChunksComplex_.push_back(std::make_unique<ComplexSignal>(resultChunksizeComplex_));
        resultChunks_.push_back(std::make_unique<DoubleSignal>(resultChunksize_));
        resultChunksComplex_.push_back(std::make_unique<ComplexSignal>(resultChunksizeComplex_));
    }
    // make one forward plan per signal chunk, and one for the patch
    // Also backward plans for the xcorr chunks
    forwardPlans_.emplace_back(new FftForwardPlan(paddedPatch_, paddedPatchComplex_));
    for (std::size_t i = 0; i < inputChunks_.size(); i++) {
        forwardPlans_.emplace_back(new FftForwardPlan(*inputChunks_.at(i), *inputChunksComplex_.at(i)));
        backwardPlans_.emplace_back(new FftBackwardPlan(*resultChunksComplex_.at(i), *resultChunks_.at(i)));
    }
}

auto OverlapSaveConvolver::executeConv() -> void
{
    execute(false);
    state_ = State::Conv;
}
auto OverlapSaveConvolver::executeXcorr() -> void
{
    execute(true);
    state_ = State::Xcorr;
}

auto OverlapSaveConvolver::printChunks(std::string const& name) -> void
{
    checkLastExecutedNotNull("printChunks");
    for (std::size_t i = 0; i < resultChunks_.size(); i++) {
        resultChunks_.at(i)->print(name + "_chunk_" + std::to_string(i));
    }
}

// This method implements step 6 of the overlap-save algorithm. In convolution, the first (P-1)
// samples of each chunk are discarded, in xcorr the last (P-1) ones. Therefore, depending on the
// current _state_, the corresponding method is used. USAGE:
// Every time it is called, this function returns a new DoubleSignal instance of size
// len(signal)+len(patch)-1. If the last operation performed was executeConv(), this function
// will return the  convolution of signal and patch. If the last operation performed was
// executeXcorr(), the result will contain the cross-correlation. If none of them was performed
// at the moment of calling this function, an exception will be thrown.
// The indexing will start with the most negative relation, and increase accordingly. Which means:
//   given S:=len(signal), P:=len(patch), T:=S+P-1
// for 0 <= i < T, result[i] will hold dot_product(patch, signal[i-(P-1) : i])
//   where patch will be "reversed" if the convolution was performed. For example:
// Signal :=        [1 2 3 4 5 6 7]    Patch = [1 1 1]
// Result[0] =  [1 1 1]                        => 1*1         = 1  // FIRST ENTRY
// Result[1] =    [1 1 1]                      => 1*1+1*2     = 3
// Result[2] =      [1 1 1]                    => 1*1+1*2+1*3 = 8  // FIRST NON-NEG ENTRY AT P-1
//   ...
// Result[8] =                  [1 1 1]        => 1*7         = 7  // LAST ENTRY
// Note that the returned signal object takes care of its own memory, so no management is needed.
auto OverlapSaveConvolver::extractResult() -> DoubleSignal
{
    // make sure that an operation was called before
    checkLastExecutedNotNull("extractResult");
    // set the offset for the corresponding operation (0 for xcorr).
    std::size_t discardOffset = 0;
    if (state_ == State::Conv) {
        discardOffset = resultChunksize_ - result_stride_;
    }
    // instantiate new signal to be filled with the desired info
    DoubleSignal result(resultSize_);
    double* resultArr = result.data(); // not const because of std::memcpy
    // fill!
    static std::size_t kNumChunks = resultChunks_.size();
    for (std::size_t i = 0; i < kNumChunks; i++) {
        double* xcArr = resultChunks_.at(i)->data();
        std::size_t const kBegin = i * result_stride_;
        // if the last chunk goes above resultSize_, reduce copy size. else copy_size=result_stride_
        std::size_t copySize = result_stride_;
        copySize -= (kBegin + result_stride_ > resultSize_) ? kBegin + result_stride_ - resultSize_ : 0;
        std::memcpy(resultArr + kBegin, xcArr + discardOffset, sizeof(double) * copySize);
    }

    return result;
}

// This private method throws an exception if _state_ is Uninitialized, because that
// means that some "getter" has ben called before any computation has been performed.
auto OverlapSaveConvolver::checkLastExecutedNotNull(std::string const& methodName) -> void
{
    if (state_ == State::Uninitialized) {
        throw std::runtime_error(std::string("[ERROR] OverlapSaveConvolver.") + methodName
            + "() can't be called before executeXcorr() or executeConv()!"
            + " No meaningful data has been computed yet.");
    }
}

// This private method implements steps 3,4,5 of the algorithm. If the given flag is false,
// it will perform a convolution (4a), and a cross-correlation (4b) otherwise.
// Note the parallelization with OpenMP, which increases performance in supporting CPUs.
auto OverlapSaveConvolver::execute(const bool crossCorrelate) -> void
{
    auto operation = (crossCorrelate) ? spectralCorrelation : spectralConvolution;
    // do ffts

    for (auto& forwardPlan : forwardPlans_) {
        forwardPlan->execute();
    }

    // multiply spectra
    for (std::size_t i = 0; i < resultChunks_.size(); i++) {
        operation(*inputChunksComplex_.at(i), this->paddedPatchComplex_, *resultChunksComplex_.at(i));
    }

    // do iffts
    for (std::size_t i = 0; i < resultChunks_.size(); i++) {
        backwardPlans_.at(i)->execute();
        *resultChunks_.at(i) /= resultChunksize_;
    }
}

/// MAIN ROUTINE

// auto main(int argc, char** argv) -> int
// {
//     if (argc < 2)
//     {
//         std::cout << "usage: " << std::string(argv[0]) << " num_rows num_cols num_levels" << '\n';
//         return 1;
//     }

//     auto const count = static_cast<std::size_t>(std::stoi(argv[1]));

//     // // do this just once to configure your system for an optimal FFT
//     // const string kWisdomPatient = "wisdom_real_dft_pow2_patient";
//     // MakeAndExportFftwWisdom(kWisdomPatient, 0, 29, FFTW_PATIENT);

//     std::size_t const kSizeS = 44100;  // 44100*10;
//     auto* sArr               = new double[kSizeS];
//     std::size_t const kSizeP = 44100;  // 44100*1;
//     auto* pArr               = new double[kSizeP];

//     double sum = 0.0f;

//     // for (auto i {0U}; i < count; ++i)
//     // {
//     // auto xi = std::array<double, 3> {1.0f, 2.0f, 3.0f};
//     // auto yi = std::array<double, 3> {0.0f, 1.0f, 0.5f};

//     auto xi = std::array {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
//     auto yi = std::array {0.0f, 1.0f, 0.5f, 1.0f, 1.0f};

//     // create a test signal
//     for (std::size_t i = 0; i < kSizeS; ++i) { sArr[i] = i + 1; }
//     // DoubleSignal s(sArr, kSizeS);
//     DoubleSignal s(xi.data(), xi.size());
//     // s.print("signal");

//     // create a test patch
//     for (std::size_t i = 0; i < kSizeP; ++i) { pArr[i] = i + 1; }
//     // DoubleSignal p(pArr, kSizeP);
//     DoubleSignal p(yi.data(), yi.size());
//     // p.print("patch");

//     // Instantiate convolver with both signals (p can't be bigger than s)
//     OverlapSaveConvolver x(s, p);

//     // CONV
//     x.executeConv();
//     // x.printChunks("conv");
//     auto const conv = x.extractResult();
//     x.extractResult().print("conv");

//     // XCORR
//     x.executeXcorr();
//     // x.printChunks("xcorr");
//     auto const xcorr = x.extractResult();

//     x.extractResult().print("xcorr");
//     // }

//     // clean memory and exit
//     delete[] sArr;
//     delete[] pArr;
//     return 0;
// }