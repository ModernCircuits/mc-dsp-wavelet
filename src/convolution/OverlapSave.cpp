#include "OverlapSave.hpp"

#include <cassert>

namespace {
constexpr auto REAL = 0U; // NOLINT
constexpr auto IMAG = 1U; // NOLINT

auto pow2Ceil(std::size_t x)
{
    return static_cast<std::size_t>(std::pow(2, std::ceil(std::log2(x))));
}
}

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

DoubleSignal::~DoubleSignal() { fftw_free(data_); }

ComplexSignal::ComplexSignal(std::size_t size)
    : Signal(fftw_alloc_complex(size), size)
{
}

ComplexSignal::~ComplexSignal() { fftw_free(data_); }

auto spectralConvolution(ComplexSignal const& a, ComplexSignal const& b, ComplexSignal& result) -> void
{
    assert((a.size() == result.size()) && (b.size() == result.size()));
    for (std::size_t i = 0; i < a.size(); ++i) {
        // a+ib * c+id = ac+iad+ibc-bd = ac-bd + i(ad+bc)
        result[i][REAL] = a[i][REAL] * b[i][REAL] - a[i][IMAG] * b[i][IMAG];
        result[i][IMAG] = a[i][IMAG] * b[i][REAL] + a[i][REAL] * b[i][IMAG];
    }
}

auto spectralCorrelation(ComplexSignal const& a, ComplexSignal const& b, ComplexSignal& result) -> void
{
    assert((a.size() == result.size()) && (b.size() == result.size()));
    for (std::size_t i = 0; i < a.size(); ++i) {
        // a * conj(b) = a+ib * c-id = ac-iad+ibc+bd = ac+bd + i(bc-ad)
        result[i][REAL] = a[i][REAL] * b[i][REAL] + a[i][IMAG] * b[i][IMAG];
        result[i][IMAG] = a[i][IMAG] * b[i][REAL] - a[i][REAL] * b[i][IMAG];
    }
}

FftForwardPlan::FftForwardPlan(DoubleSignal& fs, ComplexSignal& cs)
    : FftPlan(fftw_plan_dft_r2c_1d(fs.size(), fs.data(), cs.data(), FFTW_ESTIMATE))
{
    assert(cs.size() == (fs.size() / 2U + 1U));
}

FftBackwardPlan::FftBackwardPlan(ComplexSignal& cs, DoubleSignal& fs)
    : FftPlan(fftw_plan_dft_c2r_1d(fs.size(), cs.data(), fs.data(), FFTW_ESTIMATE))
{
    assert(cs.size() == (fs.size() / 2U + 1U));
}

OverlapSave::OverlapSave(DoubleSignal& signal, DoubleSignal& patch, std::string const& /*wisdomPath*/)
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
    assert(patchSize_ <= signalSize_);

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

auto OverlapSave::convolute() -> void
{
    execute(false);
    state_ = State::Conv;
}
auto OverlapSave::crossCorrelate() -> void
{
    execute(true);
    state_ = State::Xcorr;
}

// This method implements step 6 of the overlap-save algorithm. In convolution, the first (P-1)
// samples of each chunk are discarded, in xcorr the last (P-1) ones. Therefore, depending on the
// current _state_, the corresponding method is used. USAGE:
// Every time it is called, this function returns a new DoubleSignal instance of size
// len(signal)+len(patch)-1. If the last operation performed was convolute(), this function
// will return the  convolution of signal and patch. If the last operation performed was
// crossCorrelate(), the result will contain the cross-correlation. If none of them was performed
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
auto OverlapSave::extractResult() -> DoubleSignal
{
    assert(state_ != State::Uninitialized);

    auto result = DoubleSignal(resultSize_);
    auto* resultArr = result.data();

    auto const numChunks = resultChunks_.size();
    auto const offset = state_ == State::Conv ? resultChunksize_ - result_stride_ : std::size_t { 0 };
    for (auto i = std::size_t { 0 }; i < numChunks; i++) {
        auto* xcArr = resultChunks_.at(i)->data();
        auto const chunkBegin = i * result_stride_;

        // if the last chunk goes above resultSize_
        //  reduce copy size. else copy_size=result_stride_
        auto copySize = result_stride_;
        copySize -= (chunkBegin + result_stride_ > resultSize_) ? chunkBegin + result_stride_ - resultSize_ : 0;

        std::memcpy(resultArr + chunkBegin, xcArr + offset, sizeof(double) * copySize);
    }

    return result;
}

// This private method implements steps 3,4,5 of the algorithm. If the given flag is false,
// it will perform a convolution (4a), and a cross-correlation (4b) otherwise.
// Note the parallelization with OpenMP, which increases performance in supporting CPUs.
auto OverlapSave::execute(const bool crossCorrelate) -> void
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

        auto& chunk = *resultChunks_.at(i);
        auto divideBy = [this](auto arg) { return arg / resultChunksize_; };
        std::transform(std::begin(chunk), std::end(chunk), std::begin(chunk), divideBy);
    }
}
