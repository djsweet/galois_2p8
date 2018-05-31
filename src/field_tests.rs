// field_tests.rs: part of the test suite of the galois2p8 Crate.
// Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
// directory of this project.

use field;
use field::Field;

fn general_fields_do(func: impl Fn(field::GeneralField)) {
    for &poly in field::POLYNOMIALS.iter() {
        let test_field = field::GeneralField::new(poly);
        func(test_field);
    }
}

fn primitive_fields_do(func: impl Fn(field::PrimitivePolynomialField)) {
    for &poly in field::PRIMITIVES.iter() {
        match field::PrimitivePolynomialField::new(poly) {
            Some(test_field) => func(test_field),
            None => panic!("{:?} wasn't primitive", poly)
        }
    }
}

fn all_fields_do(
    general: impl Fn(field::GeneralField),
    primitive: impl Fn(field::PrimitivePolynomialField)
) {
    general_fields_do(general);
    primitive_fields_do(primitive);
}

fn test_field_inverses_typed(field: impl field::Field) {
    // Additive inverses
    // We're cheating here by way of knowing that
    // this is all just xor under the hood, and as a result,
    // every element is its own additive inverse.
    for i in 0..=255 {
        assert_eq!(field.sub(i, i), 0, "{:?}", field.polynomial());
    }
    // We aren't so lucky with multiplicative inverses.
    // Fields, by definition, only guarantee that they exist
    // for every item, not that it will be trivial to find them.
    let mut found_inverse = Vec::<bool>::new();
    found_inverse.resize(256, false);
    // We're skipping 0 for two reasons:
    // 1. 0 / x = 0 for all x != 0
    // 2. x / 0 is undefined.
    for i in 1..=255 {
        for j in 1..=255 {
            if field.mult(i, j) == 1 {
                found_inverse[i as usize] = true;
                let inv_i = field.div(1, i);
                let inv_j = field.div(1, j);
                assert_eq!(inv_i, j);
                assert_eq!(inv_j, i);
                break;
            }
        }
    }
    for i in 1..=255 {
        assert!(
            found_inverse[i],
            "No multiplicative inverse found for {} in {:?}!",
            i,
            field.polynomial()
        );
    }
}

#[test]
fn test_field_inverses_general() {
    general_fields_do(test_field_inverses_typed);
}

#[test]
fn test_field_inverses_primitive() {
    primitive_fields_do(test_field_inverses_typed);
}

fn test_field_associativity_typed(field: impl field::Field) {
    // Looping over all possible elements in the field is actually
    // pretty fast, usually less than a second in C even without
    // compiler optimizations, and no more than 30 seconds on
    // 2013 laptop hardware in Python.
    for a in 0..=255 {
        for b in 0..=255 {
            for c in 0..=255 {
                let a_plus_b = field.add(a, b);
                let b_plus_c = field.add(b, c);
                let ab_plus_c = field.add(a_plus_b, c);
                let a_plus_bc = field.add(a, b_plus_c);
                assert_eq!(
                    ab_plus_c,
                    a_plus_bc,
                    "({} + {} + {}) {:?}",
                    a,
                    b,
                    c,
                    field.polynomial()
                );
                let a_by_b = field.mult(a, b);
                let b_by_c = field.mult(b, c);
                let ab_by_c = field.mult(a_by_b, c);
                let a_by_bc = field.mult(a, b_by_c);
                assert_eq!(
                    ab_by_c,
                    a_by_bc,
                    "({} * {} * {}) {:?}",
                    a,
                    b,
                    c,
                    field.polynomial()
                );
            }
        }
    }
}

#[test]
fn test_field_associativity() {
    all_fields_do(
        test_field_associativity_typed,
        test_field_associativity_typed
    )
}

fn test_field_commutivity_typed(field: impl field::Field) {
    for a in 0..=255 {
        for b in 0..=255 {
            let a_plus_b = field.add(a, b);
            let b_plus_a = field.add(b, a);
            assert_eq!(
                a_plus_b,
                b_plus_a,
                "({} + {}) {:?}",
                a,
                b,
                field.polynomial()
            );
            let a_by_b = field.mult(a, b);
            let b_by_a = field.mult(b, a);
            assert_eq!(
                a_by_b,
                b_by_a,
                "{:?}",
                field.polynomial()
            );
        }
    }
}

#[test]
fn test_field_commutivity() {
    all_fields_do(
        test_field_commutivity_typed,
        test_field_commutivity_typed
    )
}

fn test_field_distributivity_typed(field: impl field::Field) {
    for a in 0..=255 {
        for b in 0..=255 {
            for c in 0..=255 {
                let a_by_b = field.mult(a, b);
                let a_by_c = field.mult(a, c);
                let b_plus_c = field.add(b, c);
                let a_by_b_plus_c = field.mult(a, b_plus_c);
                let ab_plus_ac = field.add(a_by_b, a_by_c);
                assert_eq!(
                    a_by_b_plus_c,
                    ab_plus_ac,
                    "({} * {} + {}) {:?}",
                    a,
                    b,
                    c,
                    field.polynomial()
                );
            }
        }
    }
}

#[test]
fn test_field_distributivity() {
    all_fields_do(
        test_field_distributivity_typed,
        test_field_distributivity_typed
    )
}

fn test_primitive_general_equiv_poly(poly: field::IrreducablePolynomial) {
    assert_eq!(poly.is_primitive(), true);
    let general = field::GeneralField::new(poly);
    let primitive = field::PrimitivePolynomialField::new_might_panic(poly);
    for a in 0..=255 {
        for b in 0..=255 {
            let gen_add = general.add(a, b);
            let prim_add = primitive.add(a, b);
            assert_eq!(gen_add, prim_add, "{} + {} in {:?}", a, b, poly);
            let gen_sub = general.sub(a, b);
            let prim_sub = primitive.sub(a, b);
            assert_eq!(gen_sub, prim_sub, "{} - {} in {:?}", a, b, poly);
            let gen_mult = general.mult(a, b);
            let prim_mult = primitive.mult(a, b);
            assert_eq!(gen_mult, prim_mult, "{} * {} in {:?}", a, b, poly);
        }
        for b in 1..=255 {
            let gen_div = general.div(a, b);
            let prim_div = primitive.div(a, b);
            assert_eq!(gen_div, prim_div, "{} / {} in {:?}", a, b, poly);
        }
    }
}

#[test]
fn test_primitive_general_equiv() {
    for &poly in field::PRIMITIVES.iter() {
        test_primitive_general_equiv_poly(poly);
    }
}

fn test_typed_slow_equiv_field(field: impl field::Field) {
    for a in 0..=255 {
        for b in 0..=255 {
            let field_mult = field.mult(a, b);
            let slow_mult = field::gf2_mult_mod(a, b, field.polynomial());
            assert_eq!(
                field_mult,
                slow_mult,
                "({} * {}) {:?}",
                a,
                b,
                field.polynomial()
            );
        }
    }
}

#[test]
fn test_general_slow_equiv() {
    general_fields_do(test_typed_slow_equiv_field);
}

#[test]
fn test_primitive_slow_equiv() {
    primitive_fields_do(test_typed_slow_equiv_field);
}
