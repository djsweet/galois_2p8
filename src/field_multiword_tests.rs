// field_multiword_tests.rs: part of the test suite of the galois_2p8 Crate.
// Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
// directory of this distribution.

use field::Field;

use field_tests::{
    all_fields_do
};

use proptest::{collection, num, strategy::Strategy};

fn test_ptr_ops(
    dst: &[u8],
    src: &[u8],
    vector_op: impl Fn(*mut u8, *const u8, usize),
    scalar_op: impl Fn(u8, u8) -> u8
) {
    assert_eq!(dst.len(), src.len());
    let mut test_dst = Vec::with_capacity(src.len());
    let mut ctl_dst = Vec::with_capacity(src.len());
    test_dst.extend(dst);
    let test_dst_ptr = test_dst.as_mut_ptr();
    let src_ptr = src.as_ptr();
    vector_op(test_dst_ptr, src_ptr, src.len());
    for i in 0..src.len() {
        ctl_dst.push(scalar_op(dst[i], src[i]));
    }
    assert_eq!(test_dst, ctl_dst);
}

fn test_slice_ops(
    dst: &[u8],
    src: &[u8],
    vector_op: impl Fn(&mut [u8], &[u8]),
    scalar_op: impl Fn(u8, u8) -> u8
) {
    let lower_len = src.len().min(dst.len());
    let mut test_dst = Vec::with_capacity(lower_len);
    let mut ctl_dst = Vec::with_capacity(lower_len);
    test_dst.extend(dst);
    vector_op(&mut test_dst, src);
    for i in 0..lower_len {
        ctl_dst.push(scalar_op(dst[i], src[i]));
    }
    for i in lower_len..dst.len() {
        ctl_dst.push(dst[i]);
    }
    assert_eq!(test_dst, ctl_dst);
}

fn test_slice_ops_upto(
    dst: &[u8],
    src: &[u8],
    upto: usize,
    vector_op: impl Fn(&mut [u8], &[u8], usize),
    scalar_op: impl Fn(u8, u8) -> u8
) {
    assert_eq!(dst.len(), src.len());
    assert!(upto <= dst.len());
    let mut test_dst = Vec::with_capacity(dst.len());
    let mut ctl_dst = Vec::with_capacity(dst.len());
    test_dst.extend(dst);
    ctl_dst.extend(dst);
    vector_op(&mut test_dst, src, upto);
    for i in 0..upto {
        ctl_dst[i] = scalar_op(dst[i], src[i]);
    }
    assert_eq!(test_dst, ctl_dst);
}

prop_compose! {
    fn equal_length_vecs()(
        (dst, src) in (1..1024usize).prop_flat_map(
            |len| (
                collection::vec(num::u8::ANY, len),
                collection::vec(num::u8::ANY, len)
            )
        )
    ) -> (Vec<u8>, Vec<u8>) {
        (dst, src)
    }
}

prop_compose! {
    fn equal_length_vecs_offset()(
        (dst, src, off) in (2..1024usize).prop_flat_map(
            |len| (
                collection::vec(num::u8::ANY, len),
                collection::vec(num::u8::ANY, len),
                1..len
            )
        )
    ) -> (Vec<u8>, Vec<u8>, usize) {
        (dst, src, off)
    }
}

prop_compose! {
    fn unequal_length_vecs()(
        dst in collection::vec(num::u8::ANY, 1..1024usize),
        src in collection::vec(num::u8::ANY, 1..1024usize)
    ) -> (Vec<u8>, Vec<u8>) {
        (dst, src)
    }
}

// add_ptr_len
proptest! {
    #[test]
    fn test_add_ptr_len(
        (ref dst, ref src) in equal_length_vecs()
    ) {
        all_fields_do_both!(
            |f| test_ptr_ops(
                dst,
                src,
                |d, s, l| unsafe { f.add_ptr_len(d, s, l) },
                |d, s| f.add(d, s)
            )
        );
    }
}

// add_ptr_scaled_len
proptest! {
    #[test]
    fn test_add_ptr_scaled_len(
        (ref dst, ref src) in equal_length_vecs(),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_ptr_ops(
                dst,
                src,
                |d, s, l| unsafe { f.add_ptr_scaled_len(d, s, scale, l) },
                |d, s| f.add(d, f.mult(s, scale))
            )
        );
    }
}

// add_multiword
proptest! {
    #[test]
    fn test_add_multiword(
        (ref dst, ref src) in unequal_length_vecs()
    ) {
        all_fields_do_both!(
            |f| test_slice_ops(
                dst,
                src,
                |d, s| f.add_multiword(d, s),
                |d, s| f.add(d, s)
            )
        );
    }
}

// add_multiword_len
proptest! {
    #[test]
    fn test_add_multiword_len(
        (ref dst, ref src, off) in equal_length_vecs_offset()
    ) {
        all_fields_do_both!(
            |f| test_slice_ops_upto(
                dst,
                src,
                off,
                |d, s, o| f.add_multiword_len(d, s, o),
                |d, s| f.add(d, s)
            )
        );
    }
}

// add_scaled_multiword
proptest! {
    #[test]
    fn test_add_scaled_multiword(
        (ref dst, ref src) in unequal_length_vecs(),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_slice_ops(
                dst,
                src,
                |d, s| f.add_scaled_multiword(d, s, scale),
                |d, s| f.add(d, f.mult(s, scale))
            )
        );
    }
}

// add_scaled_multiword_len
proptest! {
    #[test]
    fn test_add_scaled_multiword_len(
        (ref dst, ref src, off) in equal_length_vecs_offset(),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_slice_ops_upto(
                dst,
                src,
                off,
                |d, s, o| f.add_scaled_multiword_len(d, s, scale, o),
                |d, s| f.add(d, f.mult(s, scale))
            )
        );
    }
}

// sub_ptr_len
proptest! {
    #[test]
    fn test_sub_ptr_len(
        (ref dst, ref src) in equal_length_vecs()
    ) {
        all_fields_do_both!(
            |f| test_ptr_ops(
                dst,
                src,
                |d, s, l| unsafe { f.sub_ptr_len(d, s, l) },
                |d, s| f.sub(d, s)
            )
        );
    }
}

// sub_ptr_scaled_len
proptest! {
    #[test]
    fn test_sub_ptr_scaled_len(
        (ref dst, ref src) in equal_length_vecs(),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_ptr_ops(
                dst,
                src,
                |d, s, l| unsafe { f.sub_ptr_scaled_len(d, s, scale, l) },
                |d, s| f.sub(d, f.mult(s, scale))
            )
        );
    }
}

// sub_multiword
proptest! {
    #[test]
    fn test_sub_multiword(
        (ref dst, ref src) in unequal_length_vecs()
    ) {
        all_fields_do_both!(
            |f| test_slice_ops(
                dst,
                src,
                |d, s| f.sub_multiword(d, s),
                |d, s| f.sub(d, s)
            )
        );
    }
}

// sub_multiword_len
proptest! {
    #[test]
    fn test_sub_multiword_len(
        (ref src, ref dst, off) in equal_length_vecs_offset()
    ) {
        all_fields_do_both!(
            |f| test_slice_ops_upto(
                dst,
                src,
                off,
                |d, s, o| f.sub_multiword_len(d, s, o),
                |d, s| f.sub(d, s)
            )
        );
    }
}

// sub_scaled_multiword
proptest! {
    #[test]
    fn test_sub_scaled_multiword(
        (ref dst, ref src) in unequal_length_vecs(),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_slice_ops(
                dst,
                src,
                |d, s| f.sub_scaled_multiword(d, s, scale),
                |d, s| f.sub(d, f.mult(s, scale))
            )
        );
    }
}

// sub_scaled_multiword_len
proptest! {
    #[test]
    fn test_sub_scaled_multiword_len(
        (ref dst, ref src, off) in equal_length_vecs_offset(),
        scale in num::u8::ANY
    ) {
        all_fields_do_both!(
            |f| test_slice_ops_upto(
                dst,
                src,
                off,
                |d, s, o| f.sub_scaled_multiword_len(d, s, scale, o),
                |d, s| f.sub(d, f.mult(s, scale))
            )
        );
    }
}

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
