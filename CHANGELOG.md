# Changelog

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