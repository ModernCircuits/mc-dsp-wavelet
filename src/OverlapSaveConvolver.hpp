#pragma once

#include <fftw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

/// TYPECHECK/ANTIBUGGING

// Given a container or its beginning and end iterables, converts the container to a string
// of the form {a, b, c} (like a basic version of Python's __str__). Usage example:
// vector<string> c1({"foo", "bar"});
// vector<size_t> c2({1});
// list<double> c3({1,2,3,4,5});
// vector<bool> c4({false, true, false});
// list<int> c5;
// std::cout << IterableToString({1.23, 4.56, -789.0}) << '\n';
// std::cout << IterableToString(c1) << '\n';
// std::cout << IterableToString({"hello", "hello"}) << '\n';
// std::cout << IterableToString(c2) << '\n';
// std::cout << IterableToString(c3) << '\n';
// std::cout << IterableToString(c4) << '\n';
// std::cout << IterableToString(c5.begin(), c5.end()) << '\n';
template <typename T>
auto iterableToString(T it, T end) -> std::string
{
    std::stringstream ss;
    ss << "{";
    bool first = true;
    for (; it != end; ++it) {
        if (first) {
            ss << *it;
            first = false;
        } else {
            ss << ", " << *it;
        }
    }
    ss << "}";
    return ss.str();
}
template <typename C> // Overload IterableToString to directly accept any Collection like vector<int>
auto iterableToString(const C& c) -> std::string
{
    return IterableToString(c.begin(), c.end());
}
template <typename T> // Overload IterableToString to directly accept initializer_lists
auto iterableToString(const std::initializer_list<T> c) -> std::string
{
    return iterableToString(c.begin(), c.end());
}

// Given a container or its beginning and end iterables, checks wether all values contained in the
// iterable are equal and raises an exception if not. Usage example:
// vector<size_t> v1({});
// vector<double> v2({123.4, 123.4, 123.4});
// vector<bool> v3({false, false, false});
// vector<size_t> v4({1});
// vector<string> v5({"hello", "hello", "bye"});
// CheckAllEqual({3,3,3,3,3,3,3,3});
// CheckAllEqual(v1);
// CheckAllEqual(v2.begin(), v2.end());
// CheckAllEqual(v3);
// CheckAllEqual(v4);
// CheckAllEqual(v5.begin(), prev(v5.end()));
// CheckAllEqual(v5);
template <typename I>
auto checkAllEqual(I beg, I end, std::string const& message = "CheckAllEqual") -> void
{
    I it = beg;
    bool allEq = true;
    auto last = (it == end) ? end : std::prev(end);
    for (; it != last; ++it) {
        allEq &= (*(it) == *(std::next(it)));
        if (!allEq) {
            throw std::runtime_error(std::string("[ERROR] ") + message + " " + iterableToString(beg, end));
        }
    }
}
template <typename C>
auto checkAllEqual(const C& c, const std::string message = "CheckAllEqual") -> void
{
    CheckAllEqual(c.begin(), c.end(), message);
}
template <typename T>
auto checkAllEqual(const std::initializer_list<T> c, const std::string message = "CheckAllEqual") -> void
{
    checkAllEqual(c.begin(), c.end(), message);
}

// Raises an exception if complex_size!=(real_size/2+1), being "/" an integer division.
auto checkRealComplexRatio(std::size_t realSize, std::size_t complexSize,
    std::string const& funcName = "CheckRealComplexRatio") -> void;

// Abstract function that performs a comparation between any 2 elements, and if the comparation
// returns a truthy value raises an exception with the given message.
template <typename T, class Functor>
auto checkTwoElements(const T a, const T b, const Functor& binaryPredicate, std::string const& message) -> void
{
    if (binaryPredicate(a, b)) {
        throw std::runtime_error(std::string("[ERROR] ") + message + " " + iterableToString({ a, b }));
    }
}

// Raises an exception with the given message if a>b.
auto checkALessEqualB(std::size_t a, std::size_t b, std::string const& message = "a was greater than b!") -> void;

auto pow2Ceil(std::size_t x) -> size_t;

/// This is an abstract base class that provides some basic, type-independent functionality for
/// any container that should behave as a signal. It is not intended to be instantiated directly.
template <typename T>
struct Signal {
    /// Given a size and a reference to an array, it fills the array with <SIZE> zeros.
    /// Therefore, **IT DELETES THE CONTENTS OF THE ARRAY**. It is intended to be passed a newly
    /// allocated array by the classes that inherit from Signal, because it isn't an expensive
    /// operation and avoids memory errors due to non-initialized values.
    explicit Signal(T* data, size_t size)
        : data_(data)
        , size_(size)
    {
        memset(data_, 0, sizeof(T) * size);
    }

    [[nodiscard]] auto size() -> size_t& { return size_; }
    [[nodiscard]] auto size() const -> std::size_t const& { return size_; }
    [[nodiscard]] auto data() -> T* { return data_; }
    [[nodiscard]] auto data() const -> const T* { return data_; }

    [[nodiscard]] auto operator[](std::size_t idx) -> T& { return data_[idx]; }
    [[nodiscard]] auto operator[](std::size_t idx) const -> T& { return data_[idx]; }

    auto print(std::string const& name = "signal") -> void
    {
        std::cout << '\n';
        for (std::size_t i = 0; i < size_; ++i) {
            std::cout << name << "[" << i << "]\t=\t" << data_[i] << '\n';
        }
    }

protected:
    T* data_;
    std::size_t size_;
};

/// This class is a Signal that works on aligned double arrays allocated by FFTW.
/// It also overloads some further operators to do basic arithmetic
struct DoubleSignal : Signal<double> {
    /// the basic constructor allocates an aligned, double array, which is zeroed by the superclass
    explicit DoubleSignal(std::size_t size);
    DoubleSignal(double* data, size_t size);
    DoubleSignal(double* data, size_t size, size_t padBef, size_t padAft);
    ~DoubleSignal();

    auto operator+=(double x) -> void;
    auto operator*=(double x) -> void;
    auto operator/=(double x) -> void;
};

/// This class is a Signal that works on aligned complex (double[2]) arrays allocated by FFTW.
/// It also overloads some further operators to do basic arithmetic
struct ComplexSignal : Signal<fftw_complex> {
    /// the basic constructor allocates an aligned, double[2] array, which is zeroed by the superclass
    explicit ComplexSignal(std::size_t size);
    ~ComplexSignal();

    auto operator*=(double x) -> void;
    auto operator+=(double x) -> void;
    auto operator+=(const fftw_complex x) -> void;

    auto print(std::string const& name = "signal") -> void;
};

/// This free function takes three complex signals a,b,c of the same size and computes the complex
/// element-wise multiplication:   a+ib * c+id = ac+iad+ibc-bd = ac-bd + i(ad+bc)   The computation
/// loop isn't sent to OMP because this function itself is already expected to be called by multiple
/// threads, and it would actually slow down the process.
/// It throuws an exception if
auto spectralConvolution(ComplexSignal const& a, ComplexSignal const& b, ComplexSignal& result) -> void;

/// This function behaves identically to SpectralConvolution, but computes c=a*conj(b) instead
/// of c=a*b:         a * conj(b) = a+ib * c-id = ac-iad+ibc+bd = ac+bd + i(bc-ad)
auto spectralCorrelation(ComplexSignal const& a, ComplexSignal const& b, ComplexSignal& result) -> void;

/// This class is a simple wrapper for the memory management of the fftw plans, plus a
/// parameterless execute() method which is also a wrapper for FFTW's execute.
/// It is not expected to be used directly: rather, to be extended by specific plans, for instance,
/// if working with real, 1D signals, only 1D complex<->real plans are needed.
struct FftPlan {
    explicit FftPlan(fftw_plan p)
        : plan_(p)
    {
    }
    ~FftPlan() { fftw_destroy_plan(plan_); }
    auto execute() { fftw_execute(plan_); }

private:
    fftw_plan plan_;
};

// This forward plan (1D, R->C) is adequate to process 1D doubles (real).
struct FftForwardPlan : FftPlan {
    // This constructor creates a real->complex plan that performs the FFT(real) and saves it into the
    // complex. As explained in the FFTW docs (http://www.fftw.org/#documentation), the size of
    // the complex has to be size(real)/2+1, so the constructor will throw a runtime error if
    // this condition doesn't hold. Since the signals and the superclass already have proper
    // destructors, no special memory management has to be done.
    explicit FftForwardPlan(DoubleSignal& fs, ComplexSignal& cs);
};

// This backward plan (1D, C->R) is adequate to process spectra of 1D doubles (real).
struct FftBackwardPlan : FftPlan {
    // This constructor creates a complex->real plan that performs the IFFT(complex) and saves it
    // complex. As explained in the FFTW docs (http://www.fftw.org/#documentation), the size of
    // the complex has to be size(real)/2+1, so the constructor will throw a runtime error if
    // this condition doesn't hold. Since the signals and the superclass already have proper
    // destructors, no special memory management has to be done.
    explicit FftBackwardPlan(ComplexSignal& cs, DoubleSignal& fs);
};

/// This class performs an efficient version of the spectral convolution/cross-correlation between
/// two 1D double arrays, <SIGNAL> and <PATCH>, called overlap-save:
/// http://www.comm.utoronto.ca/~dkundur/course_info/real-time-DSP/notes/8_Kundur_Overlap_Save_Add.pdf
/// This algorithm requires that the length of <PATCH> is less or equal the length of <SIGNAL>,
/// so an exception is thrown otherwise. The algorithm works as follows:
/// given signal of length S and patch of length P, and being the conv (or xcorr) length U=S+P-1
///   1. pad the patch to X = 2*Pow2Ceil(P). FFTs with powers of 2 are the fastest.
///   2. cut the signal into chunks of size X, with an overlapping section of L=X-(P-1).
///      for that, pad the signal with (P-1) before, and with (X-U%L) after, to make it fit exactly.
///   3. Compute the forward FFT of the padded patch and of every chunk of the signal
///   4. Multiply the FFT of the padded patch with every signal chunk.
///      4a. If the operation is a convolution, perform a complex a*b multiplication
///      4b. If the operation is a cross-correlation, perform a complex a*conj(b) multiplication
///   5. Compute the inverse FFT of every result of step 4
///   6. Concatenate the resulting chunks, ignoring (P-1) samples per chunk
/// Note that steps 3,4,5 may be parallelized with some significant gain in performance.
/// In this class: X = result_chunksize, L = result_stride
struct OverlapSaveConvolver {
    /// The only constructor for the class, receives two signals and performs steps 1 and 2 of the
    /// algorithm on them. The signals are passed by reference but the class works with padded copies
    /// of them, so no care has to be taken regarding memory management.
    /// The wisdomPath may be empty, or a path to a valid wisdom file.
    /// Note that len(signal) can never be smaller than len(patch), or an exception is thrown.
    OverlapSaveConvolver(DoubleSignal& signal, DoubleSignal& patch, std::string const& wisdomPath = "");

    auto executeConv() -> void;
    auto executeXcorr() -> void;
    auto printChunks(std::string const& name = "convolver") -> void;

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
    auto extractResult() -> DoubleSignal;

private:
    // This private method throws an exception if _state_ is Uninitialized, because that
    // means that some "getter" has ben called before any computation has been performed.
    auto checkLastExecutedNotNull(std::string const& methodName) -> void;

    // This private method implements steps 3,4,5 of the algorithm. If the given flag is false,
    // it will perform a convolution (4a), and a cross-correlation (4b) otherwise.
    // Note the parallelization with OpenMP, which increases performance in supporting CPUs.
    auto execute(bool crossCorrelate) -> void;

    // grab input lengths
    std::size_t signalSize_;
    std::size_t patchSize_;
    std::size_t resultSize_;

    // make padded copies of the inputs and get chunk measurements
    DoubleSignal paddedPatch_;
    std::size_t resultChunksize_;
    std::size_t resultChunksizeComplex_;
    std::size_t result_stride_;
    ComplexSignal paddedPatchComplex_;

    // padded copy of the signal
    DoubleSignal paddedSignal_;

    // the deconstructed signal
    std::vector<std::unique_ptr<DoubleSignal>> inputChunks_;
    std::vector<std::unique_ptr<ComplexSignal>> inputChunksComplex_;

    // the corresponding chunks holding convs/xcorrs
    std::vector<std::unique_ptr<DoubleSignal>> resultChunks_;
    std::vector<std::unique_ptr<ComplexSignal>> resultChunksComplex_;

    // the corresponding plans (plus the plan of the patch)
    std::vector<std::unique_ptr<FftForwardPlan>> forwardPlans_;
    std::vector<std::unique_ptr<FftBackwardPlan>> backwardPlans_;

    // Basic state management to prevent getters from being called prematurely.
    // Also to adapt the extractResult getter, since Conv and Xcorr padding behaves differently
    enum class State {
        Uninitialized,
        Conv,
        Xcorr
    };

    State state_; // Uninitialized after instantiation, Conv/Xcorr after respective op.
};
