# galois_2p8
Basic Arithmetic over Galois (finite) fields with `2^8 == 256` members.

This library currently implements addition, subtraction, multiplication, and
division over members of a `GF(2^8) == GF(256)` field. All possible `GF(2^8)`
fields are supported by these implementations.

## Galois Fields: A Primer

### Fields
Fields are mathematical objects over which the following operations are defined:

- Addition: `a + b`
- Subtraction `a - b`
- Multiplication `a * b`
- Division `a / b`

Where Subtraction is the inverse of Addition:

    a + b - b == a
    a + b - a == b
    
And Division is the inverse of Multiplication:

    a * b / b == a
    a * b / a == b
    
With a few more properties as well:

- Associativity: `a + (b + c) == (a + b) + c` and `a * (b * c) == (a * b) * c`
- Commutativity: `a + b == b + a` and `a * b == b * a`
- Identity: `a + 0 == a` and `a * 1 == a`
- Distributivity: `a * (b + c) == a * b + a * c`

Real numbers constitute a field, but there are many others. In particular,
_finite_ fields exist.

### Finite Fields
Finite fields only contain a certain number of distinct members: operations
within these finite fields become cyclical. Finite fields are also known as
Galois fields, hence the crate name `galois_2p8`.

The canonical definition of a Galois field requires that the cardinality of the
membership set be a prime number. That is, one cannot have a finite field
containing a composite count of distinct members under the canonical definition.

Galois fields with a membership set of a composite cardinality are possible
under arithmetic extension. That is, if the cardinality of the membership set
may be composite if said cardinality is a power of a prime number. The prime
number basis is known as the _characteristic_ of the resulting field, and
the exponent is known as the _order_. In the case of this crate, the fields
implemented have a characteristic of two and an order of 8. The shorthand
notation for Galois fields of a given characteristic and order is
`GF(characteristic ^ order)`, where `characteristic ^ order` may be concretely
written as the result of the exponentiation. This crate implements operations
over `GF(2^8) == GF(256)`.

Arithmetic extensions of Galois fields treat all members as a polynomial of a
degree equal to the extension's order minus one, ensuring that there may be
exactly `characteristic ^ order` members of the field when accounting for 0.
Within this crate, field members are treated as a degree 8 polynomial, where
the coefficients of said polynomial are members of the Galois field with two
possible members: 0 and 1. The coefficients naturally map into bits, and the
polynomial as a whole maps into unsigned 8-bit integers.

Arithmetic within these algebraic extensions follows the general definition of
arithmetic over polynomials. Addition and subtraction of two polynomials results
in a polynomial whose coefficients are the addition or subtraction of the
corresponding coefficients in the source polynomial. Division is defined purely
in terms of the multiplicative inverse, and multiplication is performed using
the distributive property of fields, with one important caveat: in order to
maintain the field properties, multiplications are performed modulo an
_irreducable polynomial_.

### Irreducable Polynomials and Multiplication
A given polynomial is said to be irreducable in general if it cannot be
represented as the result of the multiplication of two other polynomials.
This is analogous to prime numbers in terms of integers. The irreducable
polynomials used in algebraic extensions have a degree equal to the order
of the extension.

Multiplication within algebraic extensions of Galois fields requires that the
usual distributive multiplication of polynomials be performed modulo an
irreducable polynomial. That is to say, the result of the overall operation
must be equal to the remainder of the non-modulo operation divided by the
target irreducable polynomial.

In general, algebraic extensions allow for the existence of multiple valid
irreducable polynomials, modulo which multiplication can be performed. Each
distinct irreducable polynomial defines its own distinct Galois field.

This crate implements all degree 8 irreducable polynomials with a characteristic
of 2 in general, with special support for _primitive polynomials_.

### Primitive Polynomials
An irreducable polynomial is said to be _primitive_ if the resulting field it
generates is structured such that there exists some basis `alpha` for which the
set of all powers, exponentiated up to the cardinality of the membership set, is
equal to the set of all members of the finite field. That is, for any arbitrary
`n` in the Galois field, `alpha ^ n` is a unique member of the same Galois field.
In the case of `GF(2 ^ n)` for arbitrary `n`, `alpha` is 2.

Multiplication and division may be represented in terms of logarithms and
exponentiation if the irreducable polynomial is a primitive polynomial.
Implementing multiplication and division in terms of logarithms and
exponentiation requires either less storage or fewer operations than using
rendered multiplication and division tables. In this crate, the multiplication
and division tables are designed to minimize storage requirements, so the
implementation of multiplication and division for single elements is faster for
fields over primitive polynomials.

Note that only a subset of the irreducable polynomials are primitive.

## The `galois_2p8` Crate
This crate contains implementations of general `GF(2^8)` field arithmetic across
all irreducable polynomials. The irreducable polynomials have been precomputed
with an auxiliary utility script (see `utils/polys.py`), and the primitive
polynomials have been identified by this precomputation as well. Special
implementations of field arithmetic are also available for these primitive
polynomials, taking advantage of the ability to accurately represent
multiplication and division in terms of exponentiation and logarithms.

### Vector Operations
The following vector operation kernels are provided in addition to
single-element arithmetic:

- Adding and subtracting vectors
- Multiplying and dividing vectors by one scalar
- Adding to or subtracting from vectors, where the non-destination vector
  is scaled by a provided scalar.

The vector kernels concerning addition and subtraction are expected to result
in optimized code if compiled with the right optimization level, but other
vector kernels default to simply being loops over scalar operations unless
SIMD support is enabled with the `"simd"` feature.

### SIMD Support
Common operations over vectors of `GF(2^8)` members lend themselves to
acceleration with SIMD intrinsics, as per James Plank's [Screaming Fast Galois
Field Arithmetic Using Intel SIMD Instructions]
(http://web.eecs.utk.edu/~plank/plank/papers/FAST-2013-GF.html). This crate
currently utilizes SSE 3 on `x86_64` only, and only if compiled with the
`"simd"` feature. AVX intrinsics are currently not used due to
[code generation bugs in `rustc`]
(https://github.com/rust-lang/rust/issues/50154),
but will be enabled when these bugs are fixed.

### Future Direction
As of the current version of this crate (`v0.1.0`), only arithmetic and vector
kernels are provided. The long-term goal of this crate is to provide optimized
matrix storage and basic matrix operations, optimized matrix inversion, and
possibly GPGPU implementations of matrix operations.