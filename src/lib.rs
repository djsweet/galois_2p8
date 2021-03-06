// lib.rs: part of the galois_2p8 Crate.
// Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
// directory of this distribution.

//! Provides operations over all `GF(2^8)` extensions.
//!
//! # Fields
//!
//! Fields are sets of values over which addition, subtraction, multiplication
//! and division are defined, such that the results of these operations are
//! still members of the defined set and the following properties hold:
//!
//! 1. Associativity: `a + (b + c) == (a + b) + c`.
//! 2. Commutativity: `a + b` == `b + a` and `a * b == b * a`.
//! 3. Identity: There exists a `0` and `1` element such that
//!    `a + 0 == a` and `a * 1 == a`.
//! 4. Additive inverse: For every `a` there is some `-a` such that
//!    `a + -a == 0`.
//! 5. Multiplicative inverse: For every `a` there is some `a^-1` such that
//!    `a * a^-1 == 1`.
//! 6. Distributivity: `a * (b + c) == a * b + b * c`.
//!
//! The set of real numbers constitutes a field, but said set is not the only
//! possible field.
//!
//! # Galois Fields
//!
//! Galois fields are finite fields: a limited number of possible members are
//! defined. The number of unique members of a Galois field is called the order
//! of the field. Galois fields only exist for integers that are also powers
//! of prime numbers. The prime base of an order is called the characteristic
//! of the field. For any Galois field with a characteristic of `c`, for all
//! members of the field `m`, `let sum = 0; for i in 0..c { sum += m };
//! sum == 0`.
//!
//! Galois fields with non-prime orders are defined in terms of polynomials with
//! the coefficients being members of the field defined by the characteristic as
//! the order. This maps directly to the representation of numbers in terms of
//! their base. The resulting field is considered an _extension_ of the field
//! generated by the characteristic.
//!
//! For example, a Galois field with an order of 256 necessarily has a
//! characteristic of 2, because 2 is the least prime base of 256. As a result,
//! `GF(256) == GF(2^8)` is defined in terms of an order 8 polynomial, taking
//! the form
//!
//! `k_8*x^8 + k_7*x^7 + k_6*x^6 + k_5*x^5 + k_4*x^4 + k_3*x^3 + k_2*x^2 + k_1*x + k_0`
//!
//! where `k_n in [0, 1] for n in [0..8]`.
//!
//! Similar to non-extension fields being defined over prime numbers, extensions
//! are defined over prime polynomials, that is, multiplication and division are
//! performed modulo some polynomial that cannot be represented as the result
//! of the multiplication of two polynomial factors. These prime polynomials are
//! known as _irreducable polynomials_, because they cannot be reduced to
//! factors.
//!
//! # Primitive Polynomials
//!
//! Consider the exponential operation `a ^ b`. In Galois fields, this is still
//! defined in terms of repeated multiplication of `a` times itself over `b`
//! iterations. Consider also a Galois field defined over an irreducable
//! polynomial `p`. If `a` is prime, but also the member of a given field, and
//! for every member `b` in the field, `a ^ b` corresponds to exactly one other
//! member of the field, then we say that `p` is a _primitive polynomial_,
//! because the entirety of the field members can be represented in terms of
//! exponents and logarithms about the base `a`.
//!
//! # `GF(2^x)` on Hardware
//!
//! The definition of finite fields with a characteristic of 2 results in addition
//! and subtraction mapping directly to bitwise XOR. Consider that:
//!
//! - `0 +/- 0 == 0; 0 XOR 0 == 0`
//! - `1 +/- 0 == 1; 1 XOR 0 == 1`
//! - `0 +/- 1 == 1; 0 XOR 1 == 1`
//! - `1 +/- 1 == 0; 1 XOR 1 == 0`
//!
//! Multiplication and division are slightly more difficult. Multiplication is
//! still defined as repeated addition of coefficients within the field, but
//! is performed modulo the primitive polynomial, which requires taking the
//! remainder of divison. Division within the field is purely defined as the
//! inverse of multiplication, with no canonical polynomial-time algorithm.
//! (This ostensible incongruity of computational complexity between
//! multiplication and division forms the basis of several areas of study
//! in cryptography).
//!
//! If the size of the field is small enough, multiplication and division can
//! be performed by way of table lookups. For non-primitive polynomials, the
//! size of the lookup tables can be reduced by decomposing operations
//! according to the field properties. For example, `a * b == (a_1 + a_2) * b
//! == a_1 * b + a_2 * b`, where a_1 represents the upper bits of `a` and `a_2`
//! represents the lower bits of `a`. As a result, the required storage space
//! is reduced by a factor of `2 ^ (x - 5)`. For primitive polynomials,
//! the operations can be performed in terms of exponentials and logarithms,
//! which only requires two to five times as many entries as members of the
//! field, depending on the implementation. `galois_2p8` currently uses
//! tables requiring three times the size of the field for multiplication
//! and division for primitive polynomials.
//!
//! # The `galois_2p8` Crate
//!
//! This crate implements `GF(2^8)` arithmetic for all possible irreducable
//! polynomials. The possible irreducable polynomials are represented as members
//! of the [`IrreducablePolynomial`] enumeration. Each [`IrreducablePolynomial`]
//! is passed as an argument to the constructor of a struct implementing the
//! [`Field`] trait. The [`GeneralField`] struct implements the [`Field`] trait
//! over all [`IrreducablePolynomial`]s, even non-primitive ones. The
//! [`PrimitivePolynomialField`] struct and implementation are specialized for
//! primitive polynomials.
//!
//! Future releases of this crate will support optimized linear algebra
//! primitives implemented across potentially heterogeneous compute
//! environments.
//!
//! [`IrreducablePolynomial`]: field/enum.IrreducablePolynomial.html
//! [`Field`]: field/trait.Field.html
//! [`GeneralField`]: field/struct.GeneralField.html
//! [`PrimitivePolynomialField`]: field/struct.PrimitivePolynomialField.html

// "The general internet" says that this _has_ to happen before
// more module definitions. It'd be real nice if the spec warned us. :/
#[cfg(test)]
#[macro_use]
extern crate proptest;
#[cfg(test)]
#[macro_use]
extern crate lazy_static;

pub use field::{
    IrreducablePolynomial,
    Field,
    GeneralField,
    PrimitivePolynomialField
};

pub mod field;

#[cfg(test)]
#[macro_use]
mod field_tests;
#[cfg(test)]
mod field_multiword_tests;
