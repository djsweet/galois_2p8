// field_multiword_tests.rs: part of the test suite of the galois_2p8 Crate.
// Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
// directory of this distribution.

use field::Field;

use field_tests::{
    all_fields_do
};

use proptest::{collection, num, strategy::Strategy};

// add_ptr_len
// add_ptr_scaled_len
// add_multiword
// add_multiword_len
// add_scaled_multiword
// add_scaled_multiword_len
// sub_ptr_len
// sub_ptr_scaled_len
// sub_multiword
// sub_multiword_len
// sub_scaled_multiword
// sub_scaled_multiword_len

// mult_ptr_len
fn test_mult_ptr_len_params(
    field: &impl Field,
    src: &[u8],
    scale: u8
) {
    let mut test_dst = Vec::with_capacity(src.len());
    let mut ctl_dst = Vec::with_capacity(src.len());
    test_dst.extend(src);
    unsafe {
        ctl_dst.set_len(src.len());
        let test_dst_ptr = test_dst.as_mut_ptr();
        field.mult_ptr_len(test_dst_ptr, scale, src.len());
    }
    for i in 0..src.len() {
        ctl_dst[i] = field.mult(src[i], scale);
    }
    assert_eq!(test_dst, ctl_dst);
}

proptest!{
    #[test]
    fn test_mult_ptr_len(
        ref src in collection::vec(num::u8::ANY, 1..1024usize),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_mult_ptr_len_params(f, src, scale)
        );
    }
}

// mult_multiword
fn test_mult_multiword_params(
    field: &impl Field,
    src: &[u8],
    scale: u8
) {
    let mut test_dst = Vec::with_capacity(src.len());
    let mut ctl_dst = Vec::with_capacity(src.len());
    test_dst.extend(src);
    field.mult_multiword(&mut test_dst, scale);
    for i in 0..src.len() {
        ctl_dst.push(field.mult(src[i], scale));
    }
    assert_eq!(test_dst, ctl_dst);
}

proptest! {
    #[test]
    fn test_mult_multiword(
        ref src in collection::vec(num::u8::ANY, 1..1024usize),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_mult_multiword_params(f, src, scale)
        );
    }
}

// div_ptr_len
fn test_div_ptr_len_params(
    field: &impl Field,
    src: &[u8],
    scale: u8
) {
    let mut test_dst = Vec::with_capacity(src.len());
    let mut ctl_dst = Vec::with_capacity(src.len());
    test_dst.extend(src);
    unsafe {
        ctl_dst.set_len(src.len());
        let test_dst_ptr = test_dst.as_mut_ptr();
        field.div_ptr_len(test_dst_ptr, scale, src.len());
    }
    for i in 0..src.len() {
        ctl_dst[i] = field.div(src[i], scale);
    }
    assert_eq!(test_dst, ctl_dst);
}

proptest! {
    #[test]
    fn test_div_ptr_len(
        ref src in collection::vec(num::u8::ANY, 1..1024usize),
        scale in num::u8::ANY.prop_filter(
            "Denominator can't be zero",
            |&n| n > 0
        )
    ) {
        all_fields_do_both!(
            |f| test_div_ptr_len_params(f, src, scale)
        );
    }
}

// div_multiword
fn test_div_multiword_params(
    field: &impl Field,
    src: &[u8],
    scale: u8
) {
    let mut test_dst = Vec::with_capacity(src.len());
    let mut ctl_dst = Vec::with_capacity(src.len());
    test_dst.extend(src);
    field.div_multiword(&mut test_dst, scale);
    for i in 0..src.len() {
        ctl_dst.push(field.div(src[i], scale));
    }
    assert_eq!(test_dst, ctl_dst);
}

proptest! {
    #[test]
    fn test_div_multiword(
        ref src in collection::vec(num::u8::ANY, 1..1024usize),
        scale in num::u8::ANY.prop_filter(
            "Denominator can't be zero",
            |&n| n > 0
        )
    ) {
        all_fields_do_both!(
            |f| test_div_multiword_params(f, src, scale)
        );
    }
}
