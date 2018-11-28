# Changelog

## [0.1.1] - 2018-11-28
Documentation fix only: functionality is identical
to the 0.1.0 release.

### Fixed
- `galois_2p8` documentation:
  - changed line: "is reduced by a factor of `2^(log(a) - 1)`" to
    "is reduced by a factor of `2^(x - 5)`" for space savings arising
    from the usage of `PrimitivePolynomialField`.
- `galois_2p8::field` documentation:
  - changed line: "ensures that all `2^(n-1)` values from `0` to `2^(n-1) - 1`
    are represented" to "ensures that all `2^n` values from `0` to `2^n - 1`".
    This fixes an off-by-one error in the definition of a degree `n` polynomial.
  - Added monospace formatting to `"simd"` feature.

## [0.1.0] - 2018-11-27
Initial release.

### Added
- Enumeration `IrreducablePolynomial` containing all
  valid irreducable polynomials for `GF(2^8)`
- Arithmetic implementations for both general `GF(2^8)`
  fields (`GeneralField`) and fields over primitive
  polynomials (`PrimitivePolynomialField`)
- (Optional) SIMD-accelerated vector operations in
  both field implementations if compiled with the
  `"simd"` feature