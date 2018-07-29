// sbox.rs: Generate the Rijndael S-box.
// A usage example of, and part of, the galois_2p8 Crate.
// Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
// directory of this distribution.

extern crate galois_2p8;

use galois_2p8::Field;

fn rotate(n: u8, shift: usize) -> u8 {
    if shift > 8 {
        return n;
    }
    (n << shift) | (n >> (8 - shift))
}

fn main() {
    // The Rijndael field is generated with the polynomial
    // x^8 + x^4 + x^3 + x^1 + 1, which is Poly84310. This polynomial is
    // _not_ primitive, so we have to use the GeneralField implementation
    // instead.
    let rijndael_field = galois_2p8::GeneralField::new(
        galois_2p8::IrreducablePolynomial::Poly84310
    );
    // The S-box has an entry for every possible 8-bit number.
    let mut sbox = Vec::with_capacity(256);
    for n in 0..=255 {
        // Each entry in the S-box is an affine transformation of the
        // multiplicative inverse of the corresponding entry offset,
        // with the exception of the 0th entry, which is defined as 0.
        let mut m_inv = if n == 0 { 0 } else { rijndael_field.div(1, n) };
        let mut result = 0;
        for _ in 0..5 {
            result ^= m_inv;
            // We have to rotate, not just shift, due to the transformation
            // definition.
            m_inv = rotate(m_inv, 1);
        }
        // The final step of creating an S-box entry is adding the column vector
        // [0, 1, 1, 0, 0, 0, 1, 1] to the polynomial entry within GF(2 ^ 8).
        // This is equivalent to simply adding the corresponding number within
        // the field.
        // 0xc6 == 0b01100011
        result = rijndael_field.add(result, 0x63);
        sbox.push(result);
    }
    // Print entries in a 16x16 grid.
    let mut column = 0;
    for ent in sbox.iter() {
        print!("{:02x}", ent);
        if column == 15 {
            println!("");
            column = 0;
        } else {
            print!(", ");
            column += 1;
        }
    }
}
