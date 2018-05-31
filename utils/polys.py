#!/usr/bin/env python
# polys.py: a utility script for the development of the galois2p8 Crate.
# Copyright 2018 Daniel Sweet. See the COPYRIGHT file at the top-level
# directory of this project.

# Determine the irreducable polynomials of GF(2^8). We are going to attempt
# to support all of them.
from __future__ import print_function

# Can't use the old log(x)/log(2) trick here on platforms with Sun Fast Math
# (which includes Solaris, most of the BSDs, macOS, & such):
# log((1 << 48) - 1) was returning 48 even though it's supposed to be 47.
# I don't know if glibc has the same issue.
def log2_floor(x):
    accum = 0
    while x > 0:
        x >>= 1
        if x > 0:
            accum += 1
    return accum

def gf2_mult(left, right):
    result = 0
    for _ in range(0, 8):
        result <<= 1
        if left & 0x80:
            result ^= right
        left <<= 1
    return result

def gf2_mod(num, poly):
    degree = log2_floor(num)
    while degree > 7:
        deg_diff = degree - 8
        num ^= (poly << deg_diff)
        degree = log2_floor(num)
    return num

COMPOSITES = [False] * 256

def find_primes():
    global COMPOSITES
    for left in range(0, 256):
        for right in range(0, 256):
            res = gf2_mult(left, right)
            if res < 0x100:
                continue
            if res >= 0x200:
                break # right
            res -= 0x100
            # If either of the left or right is 1, that's still not a
            # reducable polynomial. You're just performing the identity
            # operation there.
            COMPOSITES[res] = left != 1 and right != 1
    # We've now filled the COMPOSITES table. We don't need to process
    # it further, because its entries are all between 0x100 and 0x200
    # exclusive; any polynomial in that range multiplied by any other 
    # polynomial except for 0x1 will put the result out of that range.
    ret = []
    for i in range(0, len(COMPOSITES)):
        if not COMPOSITES[i]:
            ret.append(i)
    return ret

def find_primitives(polys):
    prims = set()
    # An irreducable polynomial f(x) is primitive iff the smallest
    # polynomial p(x) congruent to 0 mod f(x) is x^(p^m - 1) - 1, or in
    # our case 2 ^ 255 - 1. Thankfully, Python does bignums implicitly.
    for poly in polys:
        fullpoly = 0x100 + poly
        was_primitive = True
        for x in range(1, 255):
            mod_test = (1 << x) - 1
            if gf2_mod(mod_test, fullpoly) == 0:
                was_primitive = False
                break
        if was_primitive:
            prims.add(poly)
    return prims

def print_poly(num):
    powers = ["8"]
    cur_pow = 7
    onum = num
    while cur_pow >= 0:
        if num & 0x80:
            powers.append("{}".format(cur_pow))
        cur_pow -= 1
        num <<= 1
    return "Poly{}".format("".join(powers))

if __name__ == '__main__':
    primes = find_primes()
    primitives = find_primitives(primes)
    print("pub enum IrreducablePolynomial {")
    for prime in primes:
        if prime in primitives:
            trailer = ", // Primitive polynomial"
        else:
            trailer = ","
        print(
            "    {:11}".format(print_poly(prime)) 
            + " = " 
            + hex(prime) 
            + trailer
        )
    print("}")
    print()
    print(
        "pub const POLYNOMIALS: [IrreducablePolynomial; %s] = [" % len(primes)
    )
    for prime in primes:
        print("    IrreducablePolynomial::" + print_poly(prime) + ",")
    print("];")
    print()
    print(
        "pub const PRIMITIVES: [IrreducablePolynomial; %s] = [" % \
                len(primitives)
    )
    for prime in primes:
        if prime not in primitives:
            continue
        print("    IrreducablePolynomial::" + print_poly(prime) + ",")
    print("];")
