#include "OverlapSaveConvolver.hpp"

#include <mc/core/algorithm.hpp>
#include <mc/core/bit.hpp>
#include <mc/core/cassert.hpp>
#include <mc/core/climits.hpp>
#include <mc/core/cstring.hpp>

namespace mc::dsp {
namespace {
auto pow2Ceil(std::size_t x)
{
    return static_cast<std::size_t>(std::pow(2, std::ceil(std::log2(x))));
}
}  // namespace

FloatSignal::FloatSignal(std::size_t size) : data_{fftwf_alloc_real(size)}, size_{size}
{
    std::fill(data_, data_ + size_, float{});
}

FloatSignal::FloatSignal(float* data, size_t size) : FloatSignal(size)
{
    std::copy(data, data + size, data_);
}

FloatSignal::FloatSignal(float* data, size_t size, size_t padBef, size_t padAft)
    : FloatSignal(size + padBef + padAft)
{
    std::copy(data, data + size, data_ + padBef);
}

FloatSignal::~FloatSignal() { fftwf_free(data_); }

ComplexSignal::ComplexSignal(std::size_t size)
    : data_{reinterpret_cast<Complex<float>*>(fftwf_alloc_complex(size))}
    // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    , size_{size}
{
    std::fill(data_, data_ + size_, Complex<float>{});
}

ComplexSignal::~ComplexSignal() { fftwf_free(data_); }

auto spectralConvolution(
    ComplexSignal const& a,
    ComplexSignal const& b,
    ComplexSignal& result
) -> void
{
    MC_ASSERT((a.size() == result.size()) && (b.size() == result.size()));
    for (std::size_t i = 0; i < a.size(); ++i) {
        // a+ib * c+id = ac+iad+ibc-bd = ac-bd + i(ad+bc)
        result[i] = a[i] * b[i];
    }
}

auto spectralCorrelation(
    ComplexSignal const& a,
    ComplexSignal const& b,
    ComplexSignal& result
) -> void
{
    MC_ASSERT((a.size() == result.size()) && (b.size() == result.size()));
    for (std::size_t i{0}; i < a.size(); ++i) {
        // a * conj(b) = a+ib * c-id = ac-iad+ibc+bd = ac+bd + i(bc-ad)
        result[i] = a[i] * std::conj(b[i]);
    }
}

FftForwardPlan::FftForwardPlan(FloatSignal& fs, ComplexSignal& cs)
    : FftPlan(fftwf_plan_dft_r2c_1d(
        static_cast<int>(fs.size()),
        fs.data(),
        reinterpret_cast<fftwf_complex*>(cs.data()),
        FFTW_ESTIMATE
    ))
{
    MC_ASSERT(cs.size() == (fs.size() / 2U + 1U));
}

FftBackwardPlan::FftBackwardPlan(ComplexSignal& cs, FloatSignal& fs)
    : FftPlan(fftwf_plan_dft_c2r_1d(
        static_cast<int>(fs.size()),
        reinterpret_cast<fftwf_complex*>(cs.data()),
        fs.data(),
        FFTW_ESTIMATE
    ))
{
    MC_ASSERT(cs.size() == (fs.size() / 2U + 1U));
}

OverlapSaveConvolver::OverlapSaveConvolver(
    FloatSignal& signal,
    FloatSignal& patch,
    String const& /*wisdomPath*/
)
    : signalSize_{signal.size()}
    , patchSize_{patch.size()}
    , resultSize_{signalSize_ + patchSize_ - 1}
    , paddedPatch_{patch.data(), patchSize_, 0, 2 * pow2Ceil(patchSize_) - patchSize_}
    , resultChunksize_{paddedPatch_.size()}
    , resultChunksizeComplex_{resultChunksize_ / 2 + 1}
    , result_stride_{resultChunksize_ - patchSize_ + 1}
    , paddedPatchComplex_{resultChunksizeComplex_}
    , paddedSignal_{
          signal.data(),
          signalSize_,
          patchSize_ - 1,
          resultChunksize_ - (resultSize_ % result_stride_)}
{
    MC_ASSERT(patchSize_ <= signalSize_);

    // chunk the signal into strides of same size as padded patch
    // and make complex counterparts too, as well as the corresponding xcorr signals
    for (std::size_t i = 0; i <= paddedSignal_.size() - resultChunksize_;
         i += result_stride_) {
        inputChunks_.push_back(
            std::make_unique<FloatSignal>(&paddedSignal_[i], resultChunksize_)
        );
        inputChunksComplex_.push_back(
            std::make_unique<ComplexSignal>(resultChunksizeComplex_)
        );
        resultChunks_.push_back(std::make_unique<FloatSignal>(resultChunksize_));
        resultChunksComplex_.push_back(
            std::make_unique<ComplexSignal>(resultChunksizeComplex_)
        );
    }
    // make one forward plan per signal chunk, and one for the patch
    // Also backward plans for the xcorr chunks
    forwardPlans_.emplace_back(new FftForwardPlan(paddedPatch_, paddedPatchComplex_));
    for (std::size_t i = 0; i < inputChunks_.size(); i++) {
        forwardPlans_.emplace_back(
            new FftForwardPlan(*inputChunks_.at(i), *inputChunksComplex_.at(i))
        );
        backwardPlans_.emplace_back(
            new FftBackwardPlan(*resultChunksComplex_.at(i), *resultChunks_.at(i))
        );
    }
}

auto OverlapSaveConvolver::convolute() -> void
{
    execute(false);
    state_ = State::Conv;
}

auto OverlapSaveConvolver::crossCorrelate() -> void
{
    execute(true);
    state_ = State::Xcorr;
}

// This method implements step 6 of the overlap-save algorithm. In convolution, the first
// (P-1) samples of each chunk are discarded, in xcorr the last (P-1) ones. Therefore,
// depending on the current _state_, the corresponding method is used. USAGE: Every time it
// is called, this function returns a new FloatSignal instance of size
// len(signal)+len(patch)-1. If the last operation performed was convolute(), this function
// will return the  convolution of signal and patch. If the last operation performed was
// crossCorrelate(), the result will contain the cross-correlation. If none of them was
// performed at the moment of calling this function, an exception will be thrown. The
// indexing will start with the most negative relation, and increase accordingly. Which
// means:
//   given S:=len(signal), P:=len(patch), T:=S+P-1
// for 0 <= i < T, result[i] will hold dot_product(patch, signal[i-(P-1) : i])
//   where patch will be "reversed" if the convolution was performed. For example:
// Signal :=        [1 2 3 4 5 6 7]    Patch = [1 1 1]
// Result[0] =  [1 1 1]                        => 1*1         = 1  // FIRST ENTRY
// Result[1] =    [1 1 1]                      => 1*1+1*2     = 3
// Result[2] =      [1 1 1]                    => 1*1+1*2+1*3 = 8  // FIRST NON-NEG ENTRY AT
// P-1
//   ...
// Result[8] =                  [1 1 1]        => 1*7         = 7  // LAST ENTRY
// Note that the returned signal object takes care of its own memory, so no management is
// needed.
auto OverlapSaveConvolver::extractResult() -> FloatSignal
{
    MC_ASSERT(state_ != State::Uninitialized);

    auto result     = FloatSignal(resultSize_);
    auto* resultArr = result.data();

    auto const numChunks = resultChunks_.size();
    auto const offset
        = state_ == State::Conv ? resultChunksize_ - result_stride_ : std::size_t{0};
    for (auto i = std::size_t{0}; i < numChunks; i++) {
        auto* xcArr           = resultChunks_.at(i)->data();
        auto const chunkBegin = i * result_stride_;

        // if the last chunk goes above resultSize_
        // reduce copy size. else copy_size=result_stride_
        auto copySize = result_stride_;
        copySize -= (chunkBegin + result_stride_ > resultSize_)
                      ? chunkBegin + result_stride_ - resultSize_
                      : 0;

        auto const* first = xcArr + offset;
        std::copy(first, first + copySize, resultArr + chunkBegin);
    }

    return result;
}

// This private method implements steps 3,4,5 of the algorithm. If the given flag is false,
// it will perform a convolution (4a), and a cross-correlation (4b) otherwise.
// Note the parallelization with OpenMP, which increases performance in supporting CPUs.
auto OverlapSaveConvolver::execute(bool const crossCorrelate) -> void
{
    for (auto& forwardPlan : forwardPlans_) { forwardPlan->execute(); }

    auto operation = (crossCorrelate) ? spectralCorrelation : spectralConvolution;
    for (std::size_t i = 0; i < resultChunks_.size(); i++) {
        operation(
            *inputChunksComplex_.at(i),
            this->paddedPatchComplex_,
            *resultChunksComplex_.at(i)
        );
    }

    for (std::size_t i = 0; i < resultChunks_.size(); i++) {
        backwardPlans_.at(i)->execute();
        auto& chunk   = *resultChunks_.at(i);
        auto divideBy = [this](auto arg) { return arg / resultChunksize_; };
        std::transform(std::begin(chunk), std::end(chunk), std::begin(chunk), divideBy);
    }
}
}  // namespace mc::dsp
