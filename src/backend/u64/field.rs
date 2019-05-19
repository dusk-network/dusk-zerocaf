//! Field arithmetic modulo `2^252 + 27742317777372353535851937790883648493` 
//! using 64-bit limbs with 128-bit products

use core::fmt::Debug;
use std::default::Default;

use core::ops::{Index, IndexMut};
use core::ops::Add;
use core::ops::Sub;
use core::ops::Mul;

use num::Integer;

use crate::backend::u64::constants;
use crate::scalar::Ristretto255Scalar;

/// A `FieldElement` represents an element into the field 
/// `2^252 + 27742317777372353535851937790883648493`
/// 
/// In the 64-bit backend implementation, the `FieldElement is 
/// represented in radix `2^52`

#[derive(Copy, Clone, PartialOrd)]
pub struct FieldElement(pub [u64;5] );

impl Debug for FieldElement {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "FieldElement({:?})", &self.0[..])
    }
}

impl Index<usize> for FieldElement {
    type Output = u64;
    fn index(&self, _index: usize) -> &u64 {
        &(self.0[_index])
    }
}

impl IndexMut<usize> for FieldElement {
    fn index_mut(&mut self, _index: usize) -> &mut u64 {
        &mut (self.0[_index])
    }
}

impl<'a, 'b> Add<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    /// Compute `a + b` (mod l)
    fn add(self, b: &'b FieldElement) -> FieldElement {
        let mut sum = FieldElement::zero();
        let mask = (1u64 << 52) - 1;

        // a + b
        let mut carry: u64 = 0;
        for i in 0..5 {
            carry = self.0[i] + b[i] + (carry >> 52);
            sum[i] = carry & mask;
        }
        // subtract l if the sum is >= l
        sum.sub(&constants::FIELD_L)
    }
}
 
impl<'a, 'b> Sub<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    /// Compute `a - b (mod l)` 
    fn sub(self, b: &'b FieldElement) -> FieldElement {
        let mut sub = 0u64; 
        let mut difference: FieldElement = FieldElement::zero(); 
        let mask = (1u64 << 52) - 1;
        // Save wrapping_sub result. Store as a reminder on the next limb.
        for i in 0..5 { 
            sub = self.0[i].wrapping_sub(b[i] + (sub >> 63));  
            difference[i] = sub & mask;
        }
        // Conditionaly add l, if difference is negative.
        // Be aware that here `sub` tells us the most significant bit of the last limb
        // so then we know if it is greater than `l` or not.
        let underflow_mask = ((sub >> 63) ^ 1).wrapping_sub(1);
        let mut carry = 0u64;
        for i in 0..5 {
            carry = (carry >> 52) + difference[i] + (constants::FIELD_L[i] & underflow_mask);
            difference[i] = carry & mask;
        } 
        difference 
    }   
}

impl<'a, 'b> Mul<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn mul(self, _rhs: &'b FieldElement) -> FieldElement {
        let prod = FieldElement::montgomery_reduce(&FieldElement::mul_internal(self, _rhs)); 
        FieldElement::montgomery_reduce(&FieldElement::mul_internal(&prod, &constants::RR_FIELD))
    }
}

impl Default for FieldElement {
    ///Returns the default value for a FieldElement = Zero.
    fn default() -> FieldElement {
        FieldElement::zero()
    }
}

impl<'a> From<&'a Ristretto255Scalar> for FieldElement {
    /// Given a Ristretto255Scalar on canonical bytes representation
    /// get it's FieldElement equivalent value as 5 limbs and 
    /// radix-52.
    fn from(origin: &'a Ristretto255Scalar) -> FieldElement {
        let origin_bytes = origin.to_bytes();
        FieldElement::from_bytes(&origin_bytes)
    }
}

impl Into<Ristretto255Scalar> for FieldElement {
    /// Given a FieldElement reference get it's
    /// Ristretto255Scalar Equivalent on it's 
    /// canonical bytes representation.
    fn into(self) -> Ristretto255Scalar {
        Ristretto255Scalar::from_canonical_bytes(self.to_bytes()).unwrap()
    }
}

/// u64 * u64 = u128 inline func multiply helper
#[inline]
fn m(x: u64, y: u64) -> u128 {
    (x as u128) * (y as u128)
}


impl FieldElement {

    /// Construct zero.
    pub fn zero() -> FieldElement {
        FieldElement([ 0, 0, 0, 0, 0 ])
    }

    /// Construct one.
    pub fn one() -> FieldElement {
        FieldElement([ 1, 0, 0, 0, 0 ])
    }

    /// Construct -1 (mod l).
    pub fn minus_one() -> FieldElement {
        FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416])
    }

    /// Evaluate if a `FieldElement` is even or not.
    pub fn is_even(self) -> bool {
        self.0[0].is_even() 
    }

    /// Load a `FieldElement` from the low 253b   bits of a 256-bit
    /// input. So Little Endian representation in bytes of a FieldElement.
    // @TODO: Macro for Inline load8 function as has variadic arguments.
    #[warn(dead_code)]
    fn from_bytes(bytes: &[u8;32]) -> Self {

        let load8 = |input: &[u8]| -> u64 {
               (input[0] as u64)
            | ((input[1] as u64) << 8)
            | ((input[2] as u64) << 16)
            | ((input[3] as u64) << 24)
            | ((input[4] as u64) << 32)
            | ((input[5] as u64) << 40)
            | ((input[6] as u64) << 48)
            | ((input[7] as u64) << 56)
        };

        let low_52_bit_mask = (1u64 << 52) - 1;
        
        FieldElement(
        // load bits [  0, 64), no shift
        [  load8(&bytes[ 0..])        & low_52_bit_mask
        // load bits [ 48,112), shift to [ 52,112)
        , (load8(&bytes[ 6..]) >>  4) & low_52_bit_mask
        // load bits [ 96,160), shift to [104,160)
        , (load8(&bytes[12..]) >>  8) & low_52_bit_mask
        // load bits [152,216), shift to [156,216)
        , (load8(&bytes[19..]) >>  4) & low_52_bit_mask
        // load bits [192,256), shift to [208,256)
        , (load8(&bytes[24..]) >> 16) & low_52_bit_mask
        ])
    }

    /// Serialize this `FieldElement` to a 32-byte array.  The
    /// encoding is canonical.
    pub fn to_bytes(self) -> [u8; 32] {

        let mut res = [0u8; 32];

        res[0]  =  (self.0[0] >> 0)                        as u8;
        res[1]  =  (self.0[0] >> 8)                        as u8;
        res[2]  =  (self.0[0] >> 16)                       as u8;
        res[3]  =  (self.0[0] >> 24)                       as u8;
        res[4]  =  (self.0[0] >> 32)                       as u8;
        res[5]  =  (self.0[0] >> 40)                       as u8;
        // Satisfy radix 52 with the next limb value shifted according the needs
        res[6]  =  ((self.0[0] >> 48) | (self.0[1] << 4))  as u8;
        res[7]  =  (self.0[1] >> 4)                        as u8;
        res[8]  =  (self.0[1] >> 12)                       as u8;
        res[9]  =  (self.0[ 1] >> 20)                      as u8;
        res[10] =  (self.0[ 1] >> 28)                      as u8;
        res[11] =  (self.0[ 1] >> 36)                      as u8;
        res[12] =  (self.0[ 1] >> 44)                      as u8;
        res[13] =  (self.0[ 2] >>  0)                      as u8;
        res[14] =  (self.0[ 2] >>  8)                      as u8;
        res[15] =  (self.0[ 2] >> 16)                      as u8;
        res[16] =  (self.0[ 2] >> 24)                      as u8;
        res[17] =  (self.0[ 2] >> 32)                      as u8;
        res[18] =  (self.0[ 2] >> 40)                      as u8;
        res[19] = ((self.0[ 2] >> 48) | (self.0[ 3] << 4)) as u8;
        res[20] =  (self.0[ 3] >>  4)                      as u8;
        res[21] =  (self.0[ 3] >> 12)                      as u8;
        res[22] =  (self.0[ 3] >> 20)                      as u8;
        res[23] =  (self.0[ 3] >> 28)                      as u8;
        res[24] =  (self.0[ 3] >> 36)                      as u8;
        res[25] =  (self.0[ 3] >> 44)                      as u8;
        res[26] =  (self.0[ 4] >>  0)                      as u8;
        res[27] =  (self.0[ 4] >>  8)                      as u8;
        res[28] =  (self.0[ 4] >> 16)                      as u8;
        res[29] =  (self.0[ 4] >> 24)                      as u8;
        res[30] =  (self.0[ 4] >> 32)                      as u8;
        res[31] =  (self.0[ 4] >> 40)                      as u8;

        // High bit should be zero.
        debug_assert!((res[31] & 0b1000_0000u8) == 0u8);
        res
    }

    /// Compute `a * b` with the function multiplying helper
    #[inline]
    pub fn mul_internal(a: &FieldElement, b: &FieldElement) -> [u128; 9] {
        let mut res = [0u128; 9];
        // Note that this is just the normal way of performing a product.
        // We need to store the results on u128 as otherwise we'll end
        // up having overflowings.
        res[0] = m(a[0],b[0]);
        res[1] = m(a[0],b[1]) + m(a[1],b[0]);
        res[2] = m(a[0],b[2]) + m(a[1],b[1]) + m(a[2],b[0]);
        res[3] = m(a[0],b[3]) + m(a[1],b[2]) + m(a[2],b[1]) + m(a[3],b[0]);
        res[4] = m(a[0],b[4]) + m(a[1],b[3]) + m(a[2],b[2]) + m(a[3],b[1]) + m(a[4],b[0]);
        res[5] =                m(a[1],b[4]) + m(a[2],b[3]) + m(a[3],b[2]) + m(a[4],b[1]);
        res[6] =                               m(a[2],b[4]) + m(a[3],b[3]) + m(a[4],b[2]);
        res[7] =                                              m(a[3],b[4]) + m(a[4],b[3]);
        res[8] =                                                             m(a[4],b[4]);

        res
    }

    /// Compute `limbs/R` (mod l), where R is the Montgomery modulus 2^260
    #[inline]
    pub (crate) fn montgomery_reduce(limbs: &[u128; 9]) -> FieldElement {

        #[inline]
        fn adjustment_fact(sum: u128) -> (u128, u64) {
            let p = (sum as u64).wrapping_mul(constants::LFACTOR_FIELD) & ((1u64 << 52) - 1);
            ((sum + m(p,constants::FIELD_L[0])) >> 52, p)
        }

        #[inline]
        fn montg_red_res(sum: u128) -> (u128, u64) {
            let w = (sum as u64) & ((1u64 << 52) - 1);
            (sum >> 52, w)
        }

        // FIELD_L[3] = 0 so we can skip these products.
        let l = &constants::FIELD_L;
        
        // the first half computes the Montgomery adjustment factor n, and begins adding n*l to make limbs divisible by R
        let (carry, n0) = adjustment_fact(        limbs[0]);
        let (carry, n1) = adjustment_fact(carry + limbs[1] + m(n0,l[1]));
        let (carry, n2) = adjustment_fact(carry + limbs[2] + m(n0,l[2]) + m(n1,l[1]));
        let (carry, n3) = adjustment_fact(carry + limbs[3]              + m(n1,l[2]) + m(n2,l[1]));
        let (carry, n4) = adjustment_fact(carry + limbs[4] + m(n0,l[4])              + m(n2,l[2]) + m(n3,l[1]));

        // limbs is divisible by R now, so we can divide by R by simply storing the upper half as the result
        let (carry, r0) = montg_red_res(carry + limbs[5]              + m(n1,l[4])              + m(n3,l[2]) + m(n4,l[1]));
        let (carry, r1) = montg_red_res(carry + limbs[6]                           + m(n2,l[4])              + m(n4,l[2]));
        let (carry, r2) = montg_red_res(carry + limbs[7]                                        + m(n3,l[4])             );
        let (carry, r3) = montg_red_res(carry + limbs[8]                                                     + m(n4,l[4]));
        let         r4 = carry as u64;

        // result may be >= r, so attempt to subtract l
        FieldElement([r0,r1,r2,r3,r4]).sub(l)
    }

    /// Compute `a^-1 (mod l)` using the The Montgomery Modular Inverse - Revisited
    /// E. Sava¸s, C¸. K. Ko¸: https://ieeexplore.ieee.org/document/863048
    #[inline]
    pub (crate) fn inverse(a: &FieldElement) -> FieldElement {
        
        #[inline]
        fn phase1(a: &FieldElement) -> (FieldElement, u64) {
            // Declare L = 2^252 + 27742317777372353535851937790883648493
            let p = FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551]);
            let mut u = p.clone();
            let mut v = a.clone();
            let mut r = FieldElement::zero();
            let mut s = FieldElement::one();
            let two = FieldElement([2, 0, 0, 0, 0]);
            let mut k = 0u64;

            while v > FieldElement::zero() {
                match(u.is_even(), v.is_even(), u > v, v > u) {
                    // u is even
                    (true, _, _, _) => {
                        for i in 0..5 {
                            u[i] = u[i] >> 1;
                        };
                        s = &s * &two;
                    },
                    // u isn't even but v is even
                    (false, true, _, _) => {
                        for i in 0..5 {
                            v[i] = v[i] >> 1;
                        };
                        r = &r * &two;
                    },
                    // u and v aren't even and u > v
                    (false, false, true, _) => {
                        u = &u - &v;
                        for i in 0..5 {
                            u[i] = u[i] >> 1;
                        };
                        r = &r + &s;
                        s = &s * &two;
                    },
                    // u and v aren't even and v > u
                    (false, false, false, true) => {
                        v = &v - &u;
                        for i in 0..5 {
                            v[i] = v[i] >> 1;
                        };
                        s = &r + &s;
                        r = &r * &two;
                    }, 
                    (false, false, false, false) => panic!("InverseMod does not exist"),
                }

                k += 1u64;
                if r >= p { r = &r - &p; }
                }
                (&p - &r, k)
            }
            
        #[inline]
        fn phase2(r: &FieldElement, k: &u64) -> FieldElement {
            let mut rr = r.clone();
            let mut p = FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551]);

            // Maybe 253, need to review it since it's the result of the log(base 2) of `FIELD_L`
            for i in 1..(k-252) {
                match rr.is_even() {
                    true => {
                        for i in 0..5 {
                            rr[i] = rr[i] >> 1;
                        };
                    },
                    false => {
                        rr = &rr + &p;
                        for i in 0..5 {
                            rr[i] = rr[i] >> 1;
                        };
                    }
                }
            }
            rr
        }

        // Implementation of McIvor, Corinne & Mcloone, Maire & Mccanny, J.V.. (2004). 
        //Improved Montgomery modular inverse algorithm. 
        //Electronics Letters. 40. 1110 - 1112. 10.1049/el:20045610. 
        
        let (mut r, mut z) = phase1(&a.clone());
        if z == 52u64 {
            return FieldElement::montgomery_reduce(&FieldElement::mul_internal(&r, &FieldElement::one()));
        } else {
            let exp = &(2u64 * 52u64) - &z;

            r = FieldElement::montgomery_reduce(&FieldElement::mul_internal(&r, &FieldElement::one()));
        }


        unimplemented!()
    }
}



pub mod tests {

    use crate::backend::u64::field::FieldElement;
    use crate::scalar::Ristretto255Scalar;

    /// Bytes representation of `-1 (mod l) = 7237005577332262213973186563042994240857116359379907606001950938285454250988`
    pub(crate) static MINUS_ONE_BYTES: [u8; 32] = [236, 211, 245, 92, 26, 99, 18, 88, 214, 156, 247, 162, 222, 249, 222, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16];
    
    /// `A = 182687704666362864775460604089535377456991567872`
    pub static A: FieldElement = FieldElement([0, 0, 0, 2, 0]);
 
    /// `B = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static B: FieldElement = FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    /// `C = 2009874587549`
    pub static C: FieldElement = FieldElement([2009874587549, 0, 0, 0, 0]);

    /// `A + B (mod l) = 904625697166532776746648320380374280088526716493097995792780030332043239911`
    pub static A_PLUS_B: FieldElement = FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551]); 

    /// `A - B (mod l) = 6332379880165729437226538243027995370101315372437730818388241662867394146822`
    pub static A_MINUS_B: FieldElement = FieldElement([2409288332882438, 4182428486726422, 2114509, 4, 15393162788864]);

    /// `B - A (mod l) = 904625697166532776746648320014998870755800986942176787613709275418060104167`
    pub static B_MINUS_A: FieldElement = FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370491, 2199023255551]);  
    
    /// `A * B (mod l) = 918847811638530094170030839746468112210851935758749834752998326598248143582`
    pub static A_TIMES_B: FieldElement = FieldElement([2201910185007838, 1263014888683320, 1977367609994094, 4238575041099341, 2233595300724]);

    /// `A * C (mod l) = 367179375066579585494548942140953299433414959963106839625728`
    pub static A_TIMES_C: FieldElement = FieldElement([0, 0, 0, 4019749175098, 0]);

    #[test]
    fn addition_with_modulo() {
        let res = &FieldElement::minus_one() + &FieldElement::one();
        for i in 0..5 {
            assert!(res[i] == FieldElement::zero()[i]);
        }
    }

    #[test]
    fn addition_without_modulo() {
        let res = &A + &B;
        for i in 0..5 {
            assert!(res[i] == A_PLUS_B[i]);
        }
    }

    #[test]
    fn subtraction_with_mod() {
        let res = &A - &B;
        for i in 0..5 {
            assert!(res[i] == A_MINUS_B[i]);
        }
    }

    #[test]
    fn subtraction_without_mod() {
        let res = &B - &A;
        for i in 0..5 {
            assert!(res[i] == B_MINUS_A[i]);
        }
    }

    #[test]
    fn subtract_equals() {
        let res = &B - &B;
        for i in 0..5 {
            assert!(res[i] == FieldElement::zero()[i]);
        }
    }

    #[test]
    fn mul_with_modulo() {
        let res = &A * &B;
        for i in 0..5 {
            assert!(res[i] == A_TIMES_B[i]);
        }
    }

    #[test]
    fn mul_without_modulo() {
        let res = &A * &C;
        for i in 0..5 {
            assert!(res[i] == A_TIMES_C[i]);
        }
    }
 
    #[test]
    fn from_bytes_conversion() {
        let num = FieldElement::from_bytes(&MINUS_ONE_BYTES);
        for i in 0..5 {
            assert!(num[i] == FieldElement::minus_one()[i]);
        }
    }

    #[test]
    fn to_bytes_conversion() {
        let bytes = FieldElement::minus_one().to_bytes();
        for i in 0..32 {
            assert!(bytes[i] == MINUS_ONE_BYTES[i]);
        }
    }

    #[test]
    fn from_ristretto255scalar() {
        // a = `2238329342913194256032495932344128051776374960164957527413114840482143558222` = res.
        let a: Ristretto255Scalar = Ristretto255Scalar::from_canonical_bytes([0x4e, 0x5a, 0xb4, 0x34, 0x5d, 0x47, 0x08, 0x84,
                                                      0x59, 0x13, 0xb4, 0x64, 0x1b, 0xc2, 0x7d, 0x52,
                                                      0x52, 0xa5, 0x85, 0x10, 0x1b, 0xcc, 0x42, 0x44,
                                                      0xd4, 0x49, 0xf4, 0xa8, 0x79, 0xd9, 0xf2, 0x04]).unwrap();
        let a_conv = FieldElement::from(&a);
        let res = FieldElement([2330265455450702, 481909309544512, 146945097235906, 1298816433963441, 5441077225716]);

        for i in 0..5 {
            assert!(a_conv[i] == res[i]);
        }
    }

    #[test]
    fn into_ristretto255scalar() {
        // a = `2238329342913194256032495932344128051776374960164957527413114840482143558222` = res.
        let a: Ristretto255Scalar = Ristretto255Scalar::from_canonical_bytes([0x4e, 0x5a, 0xb4, 0x34, 0x5d, 0x47, 0x08, 0x84,
                                                      0x59, 0x13, 0xb4, 0x64, 0x1b, 0xc2, 0x7d, 0x52,
                                                      0x52, 0xa5, 0x85, 0x10, 0x1b, 0xcc, 0x42, 0x44,
                                                      0xd4, 0x49, 0xf4, 0xa8, 0x79, 0xd9, 0xf2, 0x04]).unwrap();
        let res: Ristretto255Scalar = FieldElement([2330265455450702, 481909309544512, 146945097235906, 1298816433963441, 5441077225716]).into();

        for i in 0..32 {
            assert!(a[i] == res[i]);
        }
    }


}