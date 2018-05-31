// field.rs: part of the galois2p8 crate.
// Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
// directory of this project.
use std::ptr;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
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
    pub fn to_u16(&self) -> u16 {
        0x100 + ((*self) as u16)
    }

    pub fn is_primitive(&self) -> bool {
        match PRIMITIVES.binary_search(self) {
            Ok(_) => true,
            Err(_) => false
        }
    }
}

pub trait Field {
    fn polynomial(&self) -> IrreducablePolynomial;
    fn mult(&self, src: u8, scale: u8) -> u8;
    fn div(&self, src: u8, scale: u8) -> u8;
    fn two_pow(&self, x: u8) -> u8;
    fn mult_two_pow(&self, scale: u8, x: u8) -> u8;

    unsafe fn add_ptr_scaled_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        scale: u8,
        len: usize
    );

    unsafe fn mult_ptr_len(
        &self,
        dst: *mut u8,
        scale: u8,
        len: usize
    );

    unsafe fn div_ptr_len(
        &self,
        dst: *mut u8,
        scale: u8,
        len: usize
    );

    // Basic arithmetic

    fn add(&self, left: u8, right: u8) -> u8 {
        left ^ right
    }

    fn sub(&self, left: u8, right: u8) -> u8 {
        left ^ right
    }

    // Multiword operations. These will eventually take advantage of SIMD
    // intrinsics, when those become stable.

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

    fn add_multiword(
        &self,
        dst: &mut [u8],
        src: &[u8]
    ) {
        let dlen = dst.len();
        let slen = src.len();
        self.add_multiword_len(dst, src, slen.min(dlen));
    }

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

    unsafe fn sub_ptr_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        len: usize
    ) {
        self.add_ptr_len(dst, src, len);
    }

    unsafe fn sub_ptr_scaled_len(
        &self,
        dst: *mut u8,
        src: *const u8,
        scale: u8,
        len: usize
    ) {
        self.add_ptr_scaled_len(dst, src, scale, len);
    }

    fn sub_multiword(
        &self,
        dst: &mut [u8],
        src: &[u8]
    ) {
        self.add_multiword(dst, src);
    }

    fn sub_multiword_len(
        &self,
        dst: &mut [u8],
        src: &[u8],
        len: usize
    ) {
        self.add_multiword_len(dst, src, len);
    }

    fn sub_scaled_multiword(
        &self,
        dst: &mut [u8],
        src: &[u8],
        scale: u8
    ) {
        self.add_scaled_multiword(dst, src, scale);
    }

    fn sub_scaled_multiword_len(
        &self,
        dst: &mut [u8],
        src: &[u8],
        scale: u8,
        len: usize
    ) {
        self.add_scaled_multiword_len(dst, src, scale, len);
    }

    fn mult_multiword(
        &self,
        dst: &mut [u8],
        scale: u8
    ) {
        unsafe {
            self.mult_ptr_len(dst.as_mut_ptr(), scale, dst.len());
        }
    }

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

// TODO: Optionally support SIMD. This can reuse mult_table and div_table
// on x86_64 with SSE 4.
pub struct GeneralField {
    modulo: IrreducablePolynomial,
    mult_table: Vec<u8>,
    div_table: Vec<u8>,
    exp_table: Vec<u8>,
    pmult_table: *const u8,
    pdiv_table: *const u8,
    pexp_table: *const u8
}

impl GeneralField {
    pub fn new(poly: IrreducablePolynomial) -> Self {
        let table_dim = 1 << 13;
        let mut mult_table = Vec::with_capacity(table_dim);
        let mut div_table = Vec::with_capacity(table_dim);
        let mut exp_table = Vec::with_capacity(1 << 8);
        let mut two_x = 1;
        unsafe {
            // We're fine with pushing to mult_table, but
            // div_table needs random access
            div_table.set_len(table_dim);
        }
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
            exp_table.push(two_x);
            two_x = gf2_mult_mod(alpha, 2, poly);
        }
        for a in 0..=255 {
            let mt_offset = a * 32;
            for b_low_offset in 0..16 {
                let prod_b_low = mult_table[mt_offset + b_low_offset];
                for b_high_offset in 16..32 {
                    let b_high = ((b_high_offset - 16) << 4) as u8;
                    let b = b_high + (b_low_offset as u8);
                    let prod_b_high = mult_table[mt_offset + b_high_offset];
                    let prod = prod_b_low ^ prod_b_high;
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
        dst: *mut u8,
        src: *const u8,
        scale: u8,
        len: usize
    ) {
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
            *dst_ptr ^= self.mult(scale, *src_ptr);
        }
    }

    unsafe fn mult_ptr_len(
        &self,
        dst: *mut u8,
        scale: u8,
        len: usize
    ) {
        if scale == 0 {
            for i in 0..len {
                *dst.offset(i as isize) = 0;
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
        if scale == 0 {
            panic!("Can't divide vector by 0");
        } else if scale != 1 {
            for i in 0..len {
                let dst_ptr = dst.offset(i as isize);
                *dst_ptr = self.div(*dst_ptr, scale);
            }
        }
    }
}

pub struct PrimitivePolynomialField {
    modulo: IrreducablePolynomial,
    exp_table: Vec<u8>,
    log_table: Vec<u8>,
    pexp_table: *const u8,
    plog_table: *const u8
}

impl PrimitivePolynomialField {
    // Returns None if the polynomial is not primitive
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
        // to use the classic mult/div tables over 16/32/64 (SSE, AVX, AVX2)
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
        // to use the classic mult/div tables over 16/32/64 (SSE, AVX, AVX2)
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
        // to use the classic mult/div tables over 16/32/64 (SSE, AVX, AVX2)
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
