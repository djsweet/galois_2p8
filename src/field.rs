// field.rs: part of the galois_2p8 Crate.
// Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
// directory of this distribution.

//! Implements arithmetic operations over all `GF(2^8)` extensions.
//!
//! Galois (finite) fields are defined in one variable modulo some
//! prime number, or over algebraic extensions, where the members
//! are polynomials with coefficients in the one-variable field modulo
//! some irreducable polynomial.
//!
//! An irreducable polynomial is analogous to a prime number: it cannot
//! be factored as the product of two or more polynomials. Performing
//! polynomial arithmetic modulo an irreducable polynomial of degree `n`
//! ensures that all `2^(n-1)` values from `0` to `2^(n-1) - 1` are represented
//! within the extended field.
//!
//! Algebraic extensions to Galois fields can be expressed as operations
//! modulo several potential irreducable polynomials, except for the special
//! case of `GF(2^2)`, which can only be represented in terms of one
//! irreducable polynomial. This crate implements field arithmetic modulo
//! all possible irreducable polynomials capable of generating `GF(2^8)`.
use std::ptr;

#[cfg(all(feature="simd", target_arch="x86_64"))]
use std::arch::x86_64;

/// Represents an irreducable polynomial of `GF(2^8)`.
///
/// Each polynomial is named according to the nonzero positions of
/// the coefficients. Each digit after the `Poly` prefix corresponds to
/// the exponent of the variable for which a nonzero coefficient is
/// present. Recall that in `GF(2^x)`, the only possible coefficients are
/// either `0` or `1`.
///
/// For example, [`Poly84310`] represents `x^8 + x^4 + x^3 + x + 1`.
///
/// [`Poly84310`]: #variant.Poly84310
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum IrreducablePolynomial {
    Poly84310   = 0x1b,
    Poly84320   = 0x1d, // Primitive polynomial
    Poly85310   = 0x2b, // Primitive polynomial
    Poly85320   = 0x2d, // Primitive polynomial
    Poly85430   = 0x39,
    Poly8543210 = 0x3f,
    Poly86320   = 0x4d, // Primitive polynomial
    Poly8643210 = 0x5f, // Primitive polynomial
    Poly86510   = 0x63, // Primitive polynomial
    Poly86520   = 0x65, // Primitive polynomial
    Poly86530   = 0x69, // Primitive polynomial
    Poly86540   = 0x71, // Primitive polynomial
    Poly8654210 = 0x77,
    Poly8654310 = 0x7b,
    Poly87210   = 0x87, // Primitive polynomial
    Poly87310   = 0x8b,
    Poly87320   = 0x8d, // Primitive polynomial
    Poly8743210 = 0x9f,
    Poly87510   = 0xa3,
    Poly87530   = 0xa9, // Primitive polynomial
    Poly87540   = 0xb1,
    Poly8754320 = 0xbd,
    Poly87610   = 0xc3, // Primitive polynomial
    Poly8763210 = 0xcf, // Primitive polynomial
    Poly8764210 = 0xd7,
    Poly8764320 = 0xdd,
    Poly8765210 = 0xe7, // Primitive polynomial
    Poly8765410 = 0xf3,
    Poly8765420 = 0xf5, // Primitive polynomial
    Poly8765430 = 0xf9,
}

/// Contains all possible irreducable polynomials for `GF(2^8)`.
///
/// This array is exposed primarily for testing purposes, but may be used
/// to enumerate all possible values in the order of their declaration.
pub const POLYNOMIALS: [IrreducablePolynomial; 30] = [
    IrreducablePolynomial::Poly84310,
    IrreducablePolynomial::Poly84320,
    IrreducablePolynomial::Poly85310,
    IrreducablePolynomial::Poly85320,
    IrreducablePolynomial::Poly85430,
    IrreducablePolynomial::Poly8543210,
    IrreducablePolynomial::Poly86320,
    IrreducablePolynomial::Poly8643210,
    IrreducablePolynomial::Poly86510,
    IrreducablePolynomial::Poly86520,
    IrreducablePolynomial::Poly86530,
    IrreducablePolynomial::Poly86540,
    IrreducablePolynomial::Poly8654210,
    IrreducablePolynomial::Poly8654310,
    IrreducablePolynomial::Poly87210,
    IrreducablePolynomial::Poly87310,
    IrreducablePolynomial::Poly87320,
    IrreducablePolynomial::Poly8743210,
    IrreducablePolynomial::Poly87510,
    IrreducablePolynomial::Poly87530,
    IrreducablePolynomial::Poly87540,
    IrreducablePolynomial::Poly8754320,
    IrreducablePolynomial::Poly87610,
    IrreducablePolynomial::Poly8763210,
    IrreducablePolynomial::Poly8764210,
    IrreducablePolynomial::Poly8764320,
    IrreducablePolynomial::Poly8765210,
    IrreducablePolynomial::Poly8765410,
    IrreducablePolynomial::Poly8765420,
    IrreducablePolynomial::Poly8765430,
];

/// Contains the primitive polynomials of `GF(2^8)`.
///
/// In `GF(2^x)` wherein arithmetic operations are performed modulo a
/// polynomial, the polynomial is said to be primitive if for some alpha,
/// each nonzero member value of the field can be uniquely represented by
/// `alpha ^ p` for `p < 2^x`, where alpha is a root of the polynomial. That
/// is, for some root of the polynomial, every member of the field can
/// be represented as the exponentiation of said root.
///
/// In `GF(2^x)`, the only nontrivial prime root of any given
/// [`IrreducablePolynomial`] is two. We say "is primitive" as a shorthand
/// for meaning that the [`IrreducablePolynomial`] is primitive assuming
/// a root of two.
///
/// The use of primitive polynomials confers an immediate performance
/// benefit for single values: we can represent multiplication and
/// division as addition and subtraction within logarithms and exponents.
/// Additionally, some usages of Galois fields, e.g. Reed-Solomon syndrome
/// coding, require primitive polynomials to function properly.
///
/// This table is used by the [`is_primitive`] method. Specifically,
///
/// ```rust
/// # use galois_2p8::field::{IrreducablePolynomial, PRIMITIVES};
/// # let poly = IrreducablePolynomial::Poly84320;
/// # assert_eq!(
/// PRIMITIVES.binary_search(&poly).is_ok() == poly.is_primitive()
/// # , true);
/// ```
///
/// [`IrreducablePolynomial`]: enum.IrreducablePolynomial.html
/// [`is_primitive`]: enum.IrreducablePolynomial.html#method.is_primitive
pub const PRIMITIVES: [IrreducablePolynomial; 16] = [
    IrreducablePolynomial::Poly84320,
    IrreducablePolynomial::Poly85310,
    IrreducablePolynomial::Poly85320,
    IrreducablePolynomial::Poly86320,
    IrreducablePolynomial::Poly8643210,
    IrreducablePolynomial::Poly86510,
    IrreducablePolynomial::Poly86520,
    IrreducablePolynomial::Poly86530,
    IrreducablePolynomial::Poly86540,
    IrreducablePolynomial::Poly87210,
    IrreducablePolynomial::Poly87320,
    IrreducablePolynomial::Poly87530,
    IrreducablePolynomial::Poly87610,
    IrreducablePolynomial::Poly8763210,
    IrreducablePolynomial::Poly8765210,
    IrreducablePolynomial::Poly8765420,
];

impl IrreducablePolynomial {
    /// Converts the [`IrreducablePolynomial`] to its binary representation.
    ///
    /// In `GF(2^x)`, coefficients in the extension polynomial may only take
    /// values of `0` or `1`. As a result, there is a natural mapping of values
    /// in `GF(2^x)` to bits. `0b1011` maps to `x^3 + x^1 + 1`, for example.
    /// The greatest degree that can be encoded by an eight-bit byte is
    /// 7, as a consequence of the least significant bit being treated as
    /// `x^0 == 1`. In order to represent a degree 8 polynomial, there must be
    /// at least 9 bits available, so it can be encoded as a 16-bit value.
    ///
    /// [`IrreducablePolynomial`]: enum.IrreducablePolynomial.html
    pub fn to_u16(&self) -> u16 {
        0x100 + ((*self) as u16)
    }

    /// Determines whether the [`IrreducablePolynomial`] is a primitive polynomial.
    ///
    /// In `GF(2^x)` wherein arithmetic operations are performed modulo a
    /// polynomial, the polynomial is said to be primitive if for some alpha,
    /// each nonzero member value of the field can be uniquely represented by
    /// `alpha ^ p` for `p < 2^x`, where alpha is a root of the polynomial. That
    /// is, for some root of the polynomial, every member of the field can
    /// be represented as the exponentiation of said root.
    ///
    /// In `GF(2^x)`, the only nontrivial prime root of any given
    /// [`IrreducablePolynomial`] is two. We say "is primitive" as a shorthand
    /// for meaning that the [`IrreducablePolynomial`] is primitive assuming
    /// a root of two.
    ///
    /// The use of primitive polynomials confers an immediate performance
    /// benefit for single values: we can represent multiplication and
    /// division as addition and subtraction within logarithms and exponents.
    /// Additionally, some usages of Galois fields, e.g. Reed-Solomon syndrome
    /// coding, require primitive polynomials to function properly.
    ///
    /// This method consults the [`PRIMITIVES`] table to determine if the
    /// [`IrreducablePolynomial`] is actually primitive. That is to say,
    ///
    /// ```rust
    /// # use galois_2p8::field::{IrreducablePolynomial, PRIMITIVES};
    /// # let poly = IrreducablePolynomial::Poly84320;
    /// # assert_eq!(
    /// poly.is_primitive() == PRIMITIVES.binary_search(&poly).is_ok()
    /// # , true);
    /// ```
    ///
    /// [`IrreducablePolynomial`]: enum.IrreducablePolynomial.html
    /// [`PRIMITIVES`]: constant.PRIMITIVES.html
    pub fn is_primitive(&self) -> bool {
        match PRIMITIVES.binary_search(self) {
            Ok(_) => true,
            Err(_) => false
        }
    }
}

/// Establishes `GF(2^8)` arithmetic for scalar and vector operands.
///
/// In all instances of `GF(2^8)`, over every possible [`IrreducablePolynomial`],
/// addition and subtraction is defined as XOR, as in `GF(2)`. Addition and
/// subtraction are accordingly provided as default implementations of
/// this trait.
///
/// Multiplication and division are more complicated, and the optimal strategy
/// for implementing them in a scalar context depends on whether the
/// [`IrreducablePolynomial`] over which the field is implemented is a primitive
/// polynomial.
///
/// Recall that if a `p: IrreducablePolynomial` is primitive, then all members of
/// the field in which operations are performed modulo `p` can be represented
/// as `2^n` for `n in [0..255]`, with the exception of `0`.
///
/// In these cases, we can represent multiplication and division as
/// addition and subtraction within logarithmic representations of the operands.
/// This requires fewer instructions to implement at the scalar level.
/// Note that this cannot be done for an [`IrreducablePolynomial`] that is not
/// also primitive. As a consequence, we provide two concrete implementations
/// of this trait: [`GeneralField`] and [`PrimitivePolynomialField`], where the
/// slightly faster logarithm arithmetic is only used in the latter.
///
/// This trait also exposes operations over vectors containing `GF(2^8)`
/// members.
///
/// Common operations over `GF(2^8)` operands can exploit long-word vector
/// operations as implemented by the target hardware. A trivial example
/// is the addition and subtraction of vectors: this is a simple bitwise
/// XOR across a very long word. This already functions as expected
/// in Rust 1.25 as a consequence of LLVM optimizations. A less trivial
/// example involves multiplication and division: vector processors require
/// a specialized long-word lookup function to implement these operations.
///
/// The `x86_64` architecture mandates SSE4.2 and earlier, as is found in the
/// earlier `x86` architecture; in SSE3, an intrinsic `_mm_shuffle_epi8` was
/// added that allows the entries of a vector register `a` to function as
/// indices of the vector register `b` in the lower four bits, effectively
/// implementing an accelerated 16-entry table lookup. This intrinsic (and the
/// AVX2 32-byte extension `_mm256_shuffle_epi8`) is not yet used, but will be
/// used when Rust gains stable SIMD intrinsic APIs.
///
/// [`IrreducablePolynomial`]: enum.IrreducablePolynomial.html
/// [`GeneralField`]: struct.GeneralField.html
/// [`PrimitivePolynomialField`]: struct.PrimitivePolynomialField.html
pub trait Field {
    /// Returns the polynomial modulo which all operations are performed.
    fn polynomial(&self) -> IrreducablePolynomial;

    /// Returns the result of `src * scale` in this field.
    fn mult(&self, src: u8, scale: u8) -> u8;

    /// Returns the result of `src / scale` in this field.
    ///
    /// Implementations of this method are expected to panic if the `scale`
    /// argument is zero. The contents of the resulting error message are
    /// not defined.
    fn div(&self, src: u8, scale: u8) -> u8;

    /// Returns the result of `2^x` in this field.
    fn two_pow(&self, x: u8) -> u8;

    /// Returns the result of `scale * 2^x` in this field.
    fn mult_two_pow(&self, scale: u8, x: u8) -> u8;

    /// Adds `scale * src[0..len]` into `dst[0..len]` in place.
    unsafe fn add_ptr_scaled_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        scale: u8,
        len: usize
    );

    /// Multiplies `dst[0..len]` by `scale` in place.
    unsafe fn mult_ptr_len(
        &self,
        dst: *mut u8,
        scale: u8,
        len: usize
    );

    /// Divides `dst[0..len]` by `scale` in place.
    ///
    /// Implementations of this method are expected to panic if the `scale`
    /// argument is zero. The contents of the resulting error message are
    /// not defined.
    unsafe fn div_ptr_len(
        &self,
        dst: *mut u8,
        scale: u8,
        len: usize
    );

    // Basic arithmetic

    /// Adds `left` and `right`, returning their sum.
    fn add(&self, left: u8, right: u8) -> u8 {
        left ^ right
    }

    /// Subtracts `right` from `left`, returning the difference.
    fn sub(&self, left: u8, right: u8) -> u8 {
        left ^ right
    }

    // Multiword operations. These will eventually take advantage of SIMD
    // intrinsics, when those become stable.

    /// Adds `src[0..len]` into `dst[0..len]`.
    unsafe fn add_ptr_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        len: usize
    ) {
        for i in 0..len {
            let dst_ptr = dst.offset(i as isize);
            let src_ptr = src.offset(i as isize);
            *dst_ptr ^= *src_ptr;
        }
    }

    /// Adds `src` into `dst` in place, over the smallest common length.
    ///
    /// The length used in operation is set to the minimum of `src.len()` and
    /// `dst.len()`.
    fn add_multiword(
        &self,
        dst: &mut [u8],
        src: &[u8]
    ) {
        let dlen = dst.len();
        let slen = src.len();
        self.add_multiword_len(dst, src, slen.min(dlen));
    }

    /// Adds `src[0..len]` into `dst[0..len]`.
    ///
    /// This method will panic if `src.len()` or `dst.len()` is less than
    /// the supplied `len` parameter.
    fn add_multiword_len(
        &self,
        dst: &mut [u8],
        src: &[u8],
        len: usize
    ) {
        assert!(len <= dst.len());
        assert!(len <= src.len());
        unsafe {
            self.add_ptr_len(dst.as_mut_ptr(), src.as_ptr(), len);
        }
    }

    /// Adds `src * scale` into `dst` in place, over the smallest common length.
    ///
    /// The length used in the operation is set to the minimum of `src.len()` and
    /// `dst.len()`.
    fn add_scaled_multiword(
        &self,
        dst: &mut [u8],
        src: &[u8],
        scale: u8
    ) {
        let dlen = dst.len();
        let slen = src.len();
        self.add_scaled_multiword_len(dst, src, scale, slen.min(dlen));
    }

    /// Adds `src[0..len] * scale` into `dst[0..len]`.
    ///
    /// This method will panic if `src.len()` or `dst.len()` is less than
    /// the supplied `len` parameter.
    fn add_scaled_multiword_len(
        &self,
        dst: &mut [u8],
        src: &[u8],
        scale: u8,
        len: usize
    ) {
        assert!(len <= dst.len());
        assert!(len <= src.len());
        unsafe {
            self.add_ptr_scaled_len(dst.as_mut_ptr(), src.as_ptr(), scale, len);
        }
    }

    /// Subtracts `src[0..len]` from `dst[0..len]` in place.
    unsafe fn sub_ptr_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        len: usize
    ) {
        self.add_ptr_len(dst, src, len);
    }

    /// Subracts `scale * src[0..len]` from `dst[0..len]` in place.
    unsafe fn sub_ptr_scaled_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        scale: u8,
        len: usize
    ) {
        self.add_ptr_scaled_len(dst, src, scale, len);
    }

    /// Subtracts `src` from `dst` in place, over the smallest common length.
    ///
    /// The length used in the operation is set to the minimum of `src.len()` and
    /// `dst.len()`.
    fn sub_multiword(
        &self,
        dst: &mut [u8],
        src: &[u8]
    ) {
        self.add_multiword(dst, src);
    }

    /// Subtracts `src[0..len]` from `dst[0..len]` in place.
    ///
    /// This method will panic if `src.len()` or `dst.len()` is less than
    /// the supplied `len` parameter.
    fn sub_multiword_len(
        &self,
        dst: &mut [u8],
        src: &[u8],
        len: usize
    ) {
        self.add_multiword_len(dst, src, len);
    }

    /// Subtracts `scale * src` from `dst` in place, over the smallest common length.
    ///
    /// The length used in the operation is set to the minimum of `src.len()` and
    /// `dst.len()`.
    fn sub_scaled_multiword(
        &self,
        dst: &mut [u8],
        src: &[u8],
        scale: u8
    ) {
        self.add_scaled_multiword(dst, src, scale);
    }

    /// Subtracts `scale * src[0..len]` from `dst[0..len]` in place.
    ///
    /// This method will panic if `src.len()` or `dst.len()` is less than
    /// the supplied `len` parameter.
    fn sub_scaled_multiword_len(
        &self,
        dst: &mut [u8],
        src: &[u8],
        scale: u8,
        len: usize
    ) {
        self.add_scaled_multiword_len(dst, src, scale, len);
    }

    /// Multiplies `dst` by `scale` in place.
    fn mult_multiword(
        &self,
        dst: &mut [u8],
        scale: u8
    ) {
        unsafe {
            self.mult_ptr_len(dst.as_mut_ptr(), scale, dst.len());
        }
    }

    /// Divides `dst` by `scale` in place.
    ///
    /// This method will panic if `scale` is zero. The contents of the
    /// resulting error message are not defined.
    fn div_multiword(
        &self,
        dst: &mut [u8],
        scale: u8
    ) {
        unsafe {
            self.div_ptr_len(dst.as_mut_ptr(), scale, dst.len());
        }
    }
}

// SIMD support. Currently limited to x86_64 and SSE.
// The simd_scale_vec and simd_scale_vec_into functions are expected
// to be used on all SIMD platforms.

// simd_scale_vec_writeback is x86_64 specific, so we don't use it in
// any client code, only in simd_scale_vec and simd_scale_vec_into.
#[cfg(all(feature="simd", target_arch="x86_64"))]
unsafe fn simd_scale_vec_writeback(
    scale_table: *const u8,
    mut dst: *mut u8,
    mut src: *const u8,
    mut len: usize,
    scale: u8,
    // write_func_32: impl Fn(*mut u8, x86_64::__m256i),
    write_func_16: impl Fn(*mut u8, x86_64::__m128i)
) -> usize {
    let scale_offset = scale_table.offset(scale as isize * 32);
    let scale_reg_lower = x86_64::_mm_loadu_si128(
        scale_offset as *const x86_64::__m128i
    );
    let scale_reg_upper = x86_64::_mm_loadu_si128(
        scale_offset.offset(16) as *const x86_64::__m128i
    );
    let mask = x86_64::_mm_set1_epi8(0x0f);
    // Since Rust 1.27, even though AVX was marked as "stable", the default
    // x86_64 ABI generates broken code that only works with the XMM registers
    // due to a code generation bug at the LLVM level. If AVX is explicitly
    // enabled at the ABI level, this functions correctly. However, because
    // the defaults are currently broken, this AVX usage will have to wait.
    //
    // if is_x86_feature_detected!("avx2") && len >= 32 {
    //     use std::slice;
    //     let srl_256 = x86_64::_mm256_broadcastsi128_si256(scale_reg_lower);
    //     let sru_256 = x86_64::_mm256_broadcastsi128_si256(scale_reg_upper);
    //     let mask_256 = x86_64::_mm256_broadcastsi128_si256(mask);
    //     while len >= 32 {
    //         let window = x86_64::_mm256_loadu_si256(
    //             src as *const x86_64::__m256i
    //         );
    //         let low_portion = x86_64::_mm256_and_si256(
    //             window,
    //             mask_256
    //         );
    //         let low_value = x86_64::_mm256_shuffle_epi8(
    //             srl_256,
    //             low_portion
    //         );
    //         let high_portion = x86_64::_mm256_and_si256(
    //             x86_64::_mm256_srli_epi16(window, 4),
    //             mask_256
    //         );
    //         let high_value = x86_64::_mm256_shuffle_epi8(
    //             sru_256,
    //             high_portion
    //         );
    //         write_func_32(dst, x86_64::_mm256_xor_si256(high_value, low_value));
    //         src = src.offset(32);
    //         dst = dst.offset(32);
    //         len -= 32;
    //     }
    // }
    while len >= 16 {
        let window = x86_64::_mm_loadu_si128(
            src as *const x86_64::__m128i
        );
        let low_portion = x86_64::_mm_and_si128(
            window,
            mask
        );
        let low_value = x86_64::_mm_shuffle_epi8(
            scale_reg_lower,
            low_portion
        );
        let high_portion = x86_64::_mm_and_si128(
            x86_64::_mm_srli_epi16(window, 4),
            mask
        );
        let high_value = x86_64::_mm_shuffle_epi8(
            scale_reg_upper,
            high_portion
        );
        write_func_16(dst, x86_64::_mm_xor_si128(high_value, low_value));
        src = src.offset(16);
        dst = dst.offset(16);
        len -= 16;
    }
    return len;
}

#[cfg(all(feature="simd", target_arch="x86_64"))]
unsafe fn simd_scale_vec(
    scale_table: *const u8,
    dst: *mut u8,
    src: *const u8,
    len: usize,
    scale: u8
) -> usize {
    use std::arch::x86_64;
    simd_scale_vec_writeback(
        scale_table,
        dst,
        src,
        len,
        scale,
        // |dst_ptr, writeback| {
        //     x86_64::_mm256_storeu_si256(
        //         dst_ptr as *mut x86_64::__m256i,
        //         writeback
        //     )
        // },
        |dst_ptr, writeback| {
            x86_64::_mm_storeu_si128(
                dst_ptr as *mut x86_64::__m128i,
                writeback
            )
        }
    )
}

#[cfg(all(feature="simd", target_arch="x86_64"))]
unsafe fn simd_scale_vec_into(
    scale_table: *const u8,
    dst: *mut u8,
    src: *const u8,
    len: usize,
    scale: u8
) -> usize {
    use std::arch::x86_64;
    simd_scale_vec_writeback(
        scale_table,
        dst,
        src,
        len,
        scale,
        // |dst, writeback| {
        //     let dst_ptr = dst as *mut x86_64::__m256i;
        //     let dst_value = x86_64::_mm256_loadu_si256(dst_ptr);
        //     let added = x86_64::_mm256_xor_si256(writeback, dst_value);
        //     x86_64::_mm256_storeu_si256(dst_ptr, added);
        // },
        |dst, writeback| {
            let dst_ptr = dst as *mut x86_64::__m128i;
            let dst_value = x86_64::_mm_loadu_si128(dst_ptr);
            let added = x86_64::_mm_xor_si128(writeback, dst_value);
            x86_64::_mm_storeu_si128(dst_ptr, added);
        }
    )
}

// SIMD operations function best over a traditional multiply/divide table.
// In the case of GeneralField, this is easy enough to support, because
// those implementations use traditional multiply/divide tables anyway.
// In the case of PrimitivePolynomialField, we'll need extra multiply/
// divide tables in addition to the exponent/logarithm tables. However,
// the multiply/divide tables can be constructed the same way for both
// implementations.

fn construct_mult_div_tables(
    poly: IrreducablePolynomial
) -> (Vec<u8>, Vec<u8>) {
    let table_dim = 1 << 13;
    let mut mult_table = Vec::with_capacity(table_dim);
    let mut div_table = Vec::with_capacity(table_dim);
    unsafe {
        // We're fine with pushing to mult table, but div_table
        // needs random access.
        div_table.set_len(table_dim);
    }
    // Handle the mult table.
    for alpha in 0..=255 {
        for low in 0..16 {
            let product = gf2_mult_mod(alpha, low, poly);
            mult_table.push(product);
        }
        for lowshift in 0..16 {
            let high = lowshift << 4;
            let product = gf2_mult_mod(alpha, high, poly);
            mult_table.push(product);
        }
    }
    // These first 32 entries are n/0, which is undefined, but worth
    // clearing out anyway in the event of some OOB access bug.
    for i in 0..32 {
        div_table[i] = 0;
    }
    // Handle the rest of the div table.
    for a in 1..=255 {
        let mt_offset = a * 32;
        // Handle 0/n cases separately in div_table for the high-quad access:
        // we're never going to see 0 fill out the high-quads from in a
        // multiplication so we have to explicitly clear it out here.
        div_table[mt_offset + 16] = 0;
        for b_low_offset in 0..16 {
            let prod_b_low = mult_table[mt_offset + b_low_offset];
            for b_high_offset in 16..32 {
                let prod_b_high = mult_table[mt_offset + b_high_offset];
                let prod = prod_b_low ^ prod_b_high;
                // In here, we have a * b = prod. So we can store
                // both prod / a = b and prod / b = a, which is
                // what we do with the div_b_offset and div_a_offset
                // into the div table.
                let b_high = ((b_high_offset - 16) << 4) as u8;
                let b = b_high + (b_low_offset as u8);
                if prod < 16 {
                    let div_b_offset = (b as usize) * 32 + prod as usize;
                    div_table[div_b_offset] = a as u8;
                    let div_a_offset = mt_offset + prod as usize;
                    div_table[div_a_offset] = b;
                } else if prod & 0x0f == 0 {
                    let prod_offset = (prod >> 4) as usize + 16;
                    let div_b_offset = (b as usize) * 32 + prod_offset;
                    div_table[div_b_offset] = a as u8;
                    let div_a_offset = mt_offset + prod_offset;
                    div_table[div_a_offset] = b;
                }
            }
        }
    }
    (mult_table, div_table)
}

/// Implements field arithmetic compatible with all [`IrreducablePolynomial`]s.
///
/// Recall that there are two strategies for optimizing field arithmetic in
/// `GF(2^8)`: accessing direct multiplication and division tables, and
/// manipulating logarithms and exponentials. The latter method requires fewer
/// operations, but is only possible if the given [`IrreducablePolynomial`] is
/// primitive.
///
/// This struct uses direct multiplication and division tables for its
/// operations, and is compatible with all [`IrreducablePolynomial`]s,
/// including primitive ones. Operations are expected to be less performant
/// than those implemented by [`PrimitivePolynomialField`], with the exception
/// of vector operations, which may be accelerated in both using multiplication
/// and division tables. See the [`Field`] documentation for more details
/// regarding vector operations.
///
/// [`IrreducablePolynomial`]: enum.IrreducablePolynomial.html
/// [`PrimitivePolynomialField`]: struct.PrimitivePolynomialField.html
/// [`Field`]: trait.Field.html
pub struct GeneralField {
    modulo: IrreducablePolynomial,
    mult_table: Vec<u8>,
    div_table: Vec<u8>,
    exp_table: Vec<u8>,
    pmult_table: *const u8,
    pdiv_table: *const u8,
    pexp_table: *const u8
}

// This is ok because instance data is not modified after construction.
unsafe impl Sync for GeneralField {}

impl GeneralField {
    /// Constructs a new `GeneralField` with all tables initialized.
    pub fn new(poly: IrreducablePolynomial) -> Self {
        let mut exp_table = Vec::with_capacity(1 << 8);
        let mut two_x = 1;
        let (mult_table, div_table) = construct_mult_div_tables(poly);
        for alpha in 0..=255 {
            exp_table.push(two_x);
            two_x = gf2_mult_mod(alpha, 2, poly);
        }
        let mut ret = Self {
            modulo: poly,
            mult_table: mult_table,
            div_table: div_table,
            exp_table: exp_table,
            pmult_table: ptr::null(),
            pdiv_table: ptr::null(),
            pexp_table: ptr::null()
        };
        ret.pmult_table = ret.mult_table.as_ptr();
        ret.pdiv_table = ret.div_table.as_ptr();
        ret.pexp_table = ret.exp_table.as_ptr();
        ret
    }
}

impl Field for GeneralField {
    fn polynomial(&self) -> IrreducablePolynomial {
        self.modulo
    }

    // For primitive polynomials, you can take the logarithms
    // and perform addition and subtraction with them, but this
    // is not true of non-primitive polynomials.

    fn mult(&self, src: u8, scale: u8) -> u8 {
        if src == 0 || scale == 0 {
            return 0;
        }
        unsafe {
            let src_low = src & 0x0f;
            let src_high = (src & 0xf0) >> 4;
            let basis = self.pmult_table.offset(scale as isize * 32);
            let prod_low = *basis.offset(src_low as isize);
            let prod_high = *basis.offset(16 + src_high as isize);
            self.add(prod_high, prod_low)
        }
    }

    fn div(&self, src: u8, scale: u8) -> u8 {
        if scale == 0 {
            panic!("Can't divide {}/0", src);
        }
        if src == 0 {
            return 0;
        }
        unsafe {
            let src_low = src & 0x0f;
            let src_high = (src & 0xf0) >> 4;
            let basis = self.pdiv_table.offset(scale as isize * 32);
            let quot_low = *basis.offset(src_low as isize);
            let quot_high = *basis.offset(16 + src_high as isize);
            self.add(quot_low, quot_high)
        }
    }

    // Special support for powers of 2

    // It is important to note that powers of two are most useful when the
    // polynomial is primitive: in these cases, powers of two from 0 through 255
    // generate the field, and 2^256 is congruent to 1, so the powers are
    // cyclical.

    // Again, that's only for primitive polynomials.

    fn two_pow(&self, x: u8) -> u8 {
        unsafe {
            *self.pexp_table.offset(x as isize)
        }
    }

    fn mult_two_pow(&self, scale: u8, x: u8) -> u8 {
        if scale == 0 {
            return 0;
        }
        if x == 0 {
            return scale;
        }
        self.mult(scale, self.two_pow(x))
    }

    unsafe fn add_ptr_scaled_len(
        &self,
        mut dst: *mut u8,
        mut src: *const u8,
        scale: u8,
        mut len: usize
    ) {
        if scale == 0 {
            return;
        }
        if scale == 1 {
            self.add_ptr_len(dst, src, len);
            return;
        }
        #[cfg(feature="simd")]
        {
            let left = simd_scale_vec_into(
                self.pmult_table,
                dst,
                src,
                len,
                scale
            );
            dst = dst.offset((len - left) as isize);
            src = src.offset((len - left) as isize);
            len = left;
        }
        while len > 0 {
            *dst ^= self.mult(scale, *src);
            dst = dst.offset(1);
            src = src.offset(1);
            len -= 1;
        }
    }

    unsafe fn mult_ptr_len(
        &self,
        mut dst: *mut u8,
        scale: u8,
        mut len: usize
    ) {
        if scale == 0 {
            for i in 0..len {
                *dst.offset(i as isize) = 0;
            }
        } else if scale != 1 {
            #[cfg(feature="simd")]
            {
                let left = simd_scale_vec(
                    self.pmult_table,
                    dst,
                    dst,
                    len,
                    scale
                );
                dst = dst.offset((len - left) as isize);
                len = left;
            }
            while len > 0 {
                *dst = self.mult(*dst, scale);
                dst = dst.offset(1);
                len -= 1;
            }
        }
    }

    unsafe fn div_ptr_len(
        &self,
        mut dst: *mut u8,
        scale: u8,
        mut len: usize
    ) {
        if scale == 0 {
            panic!("Can't divide vector by 0");
        } else if scale != 1 {
            #[cfg(feature="simd")]
            {
                let left = simd_scale_vec(
                    self.pdiv_table,
                    dst,
                    dst,
                    len,
                    scale
                );
                dst = dst.offset((len - left) as isize);
                len = left;
            }
            while len > 0 {
                *dst = self.div(*dst, scale);
                dst = dst.offset(1);
                len -= 1;
            }
        }
    }
}

/// Implements field arithmetic compatible with primitive [`IrreducablePolynomial`]s.
///
/// Recall that there are two strategies for optimizing field arithmetic in
/// `GF(2^8)`: accessing direct multiplication and division tables, and
/// manipulating logarithms and exponentials. The latter method requires fewer
/// operations, but is only possible if the given [`IrreducablePolynomial`] is
/// primitive.
///
/// This struct uses exponentiation and logarithm tables, and is only
/// compatible with primitive [`IrreducablePolynomial`]s. For an implementation
/// compatible with all [`IrreducablePolynomial`]s, see [`GeneralField`].
///
/// Note that this implementation may also use multiplication and division
/// tables for its vectorized operations. See the [`Field`] documentation
/// for more details.
///
/// [`IrreducablePolynomial`]: enum.IrreducablePolynomial.html
/// [`GeneralField`]: struct.GeneralField.html
/// [`Field`]: trait.Field.html
pub struct PrimitivePolynomialField {
    modulo: IrreducablePolynomial,
    exp_table: Vec<u8>,
    log_table: Vec<u8>,
    pexp_table: *const u8,
    plog_table: *const u8
}

// This is ok because instance data is not modified after construction.
unsafe impl Sync for PrimitivePolynomialField {}

impl PrimitivePolynomialField {
    /// Constructs a new `PrimitivePolynomialField` with all tables initialized.
    ///
    /// If the given `poly` argument is not primitive, this function returns
    /// `None`; otherwise it returns `Some(f: PrimitivePolynomialField)`.
    /// In situations where the use of `Option<PrimitivePolynomialField>` is
    /// less ideal than incurring a panic, consider [`new_might_panic`].
    ///
    /// [`new_might_panic`]: #method.new_might_panic
    pub fn new(poly: IrreducablePolynomial) -> Option<Self> {
        if !poly.is_primitive() {
            return None;
        }
        let mut exp_table = Vec::with_capacity(510);
        let mut log_table = Vec::with_capacity(255);
        unsafe {
            exp_table.set_len(510);
            log_table.set_len(256);
        }
        let mut member = 1;
        for x in 0..255 {
            exp_table[x] = member;
            exp_table[x + 255] = member;
            log_table[member as usize] = x as u8;
            member = gf2_mult_mod(member, 2, poly);
        }
        // This isn't used
        log_table[0] = 0;
        let mut ret = Self{
            modulo: poly,
            exp_table: exp_table,
            log_table: log_table,
            pexp_table: ptr::null(),
            plog_table: ptr::null()
        };
        ret.pexp_table = ret.exp_table.as_ptr();
        ret.plog_table = ret.log_table.as_ptr();
        Some(ret)
    }

    /// Constructs a new `PrimitivePolynomialField` with all tables initialized.
    ///
    /// If the given `poly` argument is not primitive, this function panics.
    /// The contents of the resulting error message are not defined.
    /// In situations where incurring a panic is less ideal than the use of
    /// `Option<PrimitivePolynomialField>`, consider [`new`].
    ///
    /// [`new`]: #method.new
    pub fn new_might_panic(poly: IrreducablePolynomial) -> Self {
        match Self::new(poly) {
            Some(f) => f,
            None => panic!("Polynomial {:?} is not primitive")
        }
    }
}

impl Field for PrimitivePolynomialField {
    fn polynomial(&self) -> IrreducablePolynomial {
        self.modulo
    }

    fn mult(&self, src: u8, scale: u8) -> u8 {
        if src == 0 || scale == 0 {
            return 0;
        }
        unsafe {
            let log_src = *self.plog_table.offset(src as isize) as isize;
            let log_scale = *self.plog_table.offset(scale as isize) as isize;
            return *self.pexp_table.offset(log_src + log_scale);
        }
    }

    fn div(&self, src: u8, scale: u8) -> u8 {
        if scale == 0 {
            panic!("Can't divide {}/0", src);
        }
        if src == 0 {
            return 0;
        }
        unsafe {
            let log_src = *self.plog_table.offset(src as isize) as isize;
            let log_scale = *self.plog_table.offset(scale as isize) as isize;
            let exp_offset = 255 + log_src - log_scale;
            return *self.pexp_table.offset(exp_offset);
        }
    }

    fn two_pow(&self, x: u8) -> u8 {
        unsafe {
            *self.pexp_table.offset(x as isize)
        }
    }

    fn mult_two_pow(&self, scale: u8, x: u8) -> u8 {
        if scale == 0 {
            return 0;
        }
        unsafe {
            let scale_log = *self.plog_table.offset(scale as isize) as isize;
            return *self.pexp_table.offset(scale_log + x as isize);
        }
    }

    unsafe fn add_ptr_scaled_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        scale: u8,
        len: usize
    ) {
        // TODO: If we support x86_64, which means we support SSE, it's faster
        // to use the classic mult/div tables over 16/32/64 (SSE, AVX2, AVX512)
        // elements than continue to use the exp/log tables.
        if scale == 0 {
            return;
        }
        if scale == 1 {
            self.add_ptr_len(dst, src, len);
            return;
        }
        for i in 0..len {
            let dst_ptr = dst.offset(i as isize);
            let src_ptr = src.offset(i as isize);
            *dst_ptr ^=  self.mult(scale, *src_ptr);
        }
    }

    unsafe fn mult_ptr_len(
        &self,
        dst: *mut u8,
        scale: u8,
        len: usize
    ) {
        // TODO: If we support x86_64, which means we support SSE, it's faster
        // to use the classic mult/div tables over 16/32/64 (SSE, AVX2, AVX512)
        // elements than continue to use the exp/log tables. 
        if scale == 0 {
            for i in 0..len {
                let dst_ptr = dst.offset(i as isize);
                *dst_ptr = 0;
            }
        } else if scale != 1 {
            for i in 0..len {
                let dst_ptr = dst.offset(i as isize);
                *dst_ptr = self.mult(*dst_ptr, scale);
            }
        }
    }

    unsafe fn div_ptr_len(
        &self,
        dst: *mut u8,
        scale: u8,
        len: usize
    ) {
        // TODO: If we support x86_64, which means we support SSE, it's faster
        // to use the classic mult/div tables over 16/32/64 (SSE, AVX2, AVX512)
        // elements than continue to use the exp/log tables.
        if scale == 0 {
            panic!("Cannot divide vector by 0");
        } else if scale != 1 {
            for i in 0..len {
                let dst_ptr = dst.offset(i as isize);
                *dst_ptr = self.div(*dst_ptr, scale);
            }
        }
    }
}

// Slow variants of GF(2^8) arithmetic, used for generating tables used
// by the fast variants.

// Exposed as pub(crate) for testing purposes.
pub(crate) fn gf2_mult_mod(
    left: u8,
    right: u8,
    modulo: IrreducablePolynomial
) -> u8 {
    gf2_mod(gf2_mult(left, right), modulo)
}

fn gf2_mult(left: u8, right: u8) -> u16 {
    let mut ret: u16 = 0;
    let mut left_working = left as u16;
    for _ in 0..8 {
        ret <<= 1;
        if left_working & 0x80 != 0 {
            ret ^= right as u16;
        }
        left_working <<= 1;
    }
    ret
}

fn gf2_mod(mut quotient: u16, poly: IrreducablePolynomial) -> u8 {
    let sub = poly.to_u16();
    let mut degree = log2_floor(quotient);
    while degree > 7 {
        let shift_by = degree - 8;
        let sub_shifted = sub << shift_by;
        quotient ^= sub_shifted;
        degree = log2_floor(quotient);
    }
    quotient as u8
}

// Note that we treat 0 as 1 for the purposes of this function
fn log2_floor(mut target: u16) -> u8 {
    let mut ret = 0;
    for _ in 0..16 {
        target >>= 1;
        if target != 0 {
            ret += 1;
        } else {
            break;
        }
    }
    ret
}
