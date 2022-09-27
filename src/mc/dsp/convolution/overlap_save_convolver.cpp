// SPDX-License-Identifier: GPL-3.0-or-later

#include "overlap_save_convolver.hpp"

#include <mc/dsp/algorithm/spectral_convolution.hpp>
#include <mc/dsp/algorithm/spectral_correlation.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/bit.hpp>
#include <mc/core/cassert.hpp>
#include <mc/core/climits.hpp>
#include <mc/core/cstring.hpp>

namespace mc {
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
    : data_{reinterpret_cast<Complex<float>*>(fftwf_alloc_complex(size))}  // NOLINT
    , size_{size}
{
    std::fill(data_, data_ + size_, Complex<float>{});
}

ComplexSignal::~ComplexSignal() { fftwf_free(data_); }

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

OverlapSaveConvolver::OverlapSaveConvolver(FloatSignal& signal, FloatSignal& patch)
    : _signalSize{signal.size()}
    , _patchSize{patch.size()}
    , _resultSize{_signalSize + _patchSize - 1}
    , _paddedPatch{patch.data(), _patchSize, 0, 2 * pow2Ceil(_patchSize) - _patchSize}
    , _resultChunksize{_paddedPatch.size()}
    , _resultChunksizeComplex{_resultChunksize / 2 + 1}
    , _result_stride{_resultChunksize - _patchSize + 1}
    , _paddedPatchComplex{_resultChunksizeComplex}
    , _paddedSignal{
          signal.data(),
          _signalSize,
          _patchSize - 1,
          _resultChunksize - (_resultSize % _result_stride)}
{
    MC_ASSERT(_patchSize <= _signalSize);

    // chunk the signal into strides of same size as padded patch
    // and make complex counterparts too, as well as the corresponding xcorr signals
    for (std::size_t i = 0; i <= _paddedSignal.size() - _resultChunksize;
         i += _result_stride) {
        _inputChunks.push_back(makeUnique<FloatSignal>(&_paddedSignal[i], _resultChunksize)
        );
        _inputChunksComplex.push_back(makeUnique<ComplexSignal>(_resultChunksizeComplex));
        _resultChunks.push_back(makeUnique<FloatSignal>(_resultChunksize));
        _resultChunksComplex.push_back(makeUnique<ComplexSignal>(_resultChunksizeComplex));
    }
    // make one forward plan per signal chunk, and one for the patch
    // Also backward plans for the xcorr chunks
    _forwardPlans.emplace_back(new FftForwardPlan(_paddedPatch, _paddedPatchComplex));
    for (std::size_t i = 0; i < _inputChunks.size(); i++) {
        _forwardPlans.emplace_back(
            new FftForwardPlan(*_inputChunks.at(i), *_inputChunksComplex.at(i))
        );
        _backwardPlans.emplace_back(
            new FftBackwardPlan(*_resultChunksComplex.at(i), *_resultChunks.at(i))
        );
    }
}

auto OverlapSaveConvolver::convolute() -> void
{
    execute(false);
    _state = State::Conv;
}

auto OverlapSaveConvolver::crossCorrelate() -> void
{
    execute(true);
    _state = State::Xcorr;
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
    MC_ASSERT(_state != State::Uninitialized);

    auto result     = FloatSignal(_resultSize);
    auto* resultArr = result.data();

    auto const numChunks = _resultChunks.size();
    auto const offset
        = _state == State::Conv ? _resultChunksize - _result_stride : std::size_t{0};
    for (auto i = std::size_t{0}; i < numChunks; i++) {
        auto* xcArr           = _resultChunks.at(i)->data();
        auto const chunkBegin = i * _result_stride;

        // if the last chunk goes above resultSize_
        // reduce copy size. else copy_size=result_stride_
        auto copySize = _result_stride;
        copySize -= (chunkBegin + _result_stride > _resultSize)
                      ? chunkBegin + _result_stride - _resultSize
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
    for (auto& forwardPlan : _forwardPlans) { forwardPlan->execute(); }

    auto operation = (crossCorrelate) ? spectralCorrelation : spectralConvolution;
    for (std::size_t i = 0; i < _resultChunks.size(); i++) {
        operation(
            *_inputChunksComplex.at(i),
            this->_paddedPatchComplex,
            *_resultChunksComplex.at(i)
        );
    }

    for (std::size_t i = 0; i < _resultChunks.size(); i++) {
        _backwardPlans.at(i)->execute();
        auto& chunk = *_resultChunks.at(i);
        ranges::transform(chunk, begin(chunk), [this](auto arg) {
            return arg / _resultChunksize;
        });
    }
}
}  // namespace mc
