Experiments for faster implementations of PAR 2.0

See the original spec here:
https://parchive.sourceforge.net/docs/specifications/parity-volume-spec/article-spec.html

Computing recovery slices is equivalent ot multiplying vectors of input data
with a transposed Vandermonde matrix V_ij = c_j ^ i where c_j is the constant
associated with the j-th input slice, and i is the exponent associated with a
recovery slice. (Only 32768 of 65536 input constants are used, for reasons
described in the Par2 spec.)

The current code is based on fast multiplication of a 65536 x 65536 transposed
Vandermonde matrix by a vector of 65536 exponents from 0 through 65535. This
multiplication is performed in O(n log n) time using the additive Fast Fourier
Transform (FFT) that is based on the novel polynomial basis introduced by David
Cantor.

As a result, the current code is only competitive when the number of input
slices is close to the maximum of 32768, and the number of recovery slices is
also large (ideally close to 32768 as well).
