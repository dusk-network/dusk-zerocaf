//! Arithmetic mod `2^249 - 15145038707218910765482344729778085401`
//! with five 52-bit unsigned limbs
//! represented in radix `2^52`.


use core::fmt::Debug;
use core::ops::{Index, IndexMut};
use core::ops::{Add, Sub, Mul, Neg};

use std::cmp::{PartialOrd, Ordering, Ord};

use num::Integer;

use crate::backend::u64::constants;
use crate::traits::Identity;
use crate::traits::ops::*;

/// The `Scalar` struct represents an Scalar over the modulo
/// `2^249 - 15145038707218910765482344729778085401` as 5 52-bit limbs
/// represented in radix `2^52`.
#[derive(Copy,Clone)]
pub struct Scalar(pub [u64; 5]);

impl Debug for Scalar {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "Scalar: {:?}", &self.0[..])
    }
}

impl Index<usize> for Scalar {
    type Output = u64;
    fn index(&self, _index: usize) -> &u64 {
        &(self.0[_index])
    }
}

impl IndexMut<usize> for Scalar {
    fn index_mut(&mut self, _index: usize) -> &mut u64 {
        &mut (self.0[_index])
    }
}

impl PartialOrd for Scalar {
    fn partial_cmp(&self, other: &Scalar) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl Ord for Scalar {
    fn cmp(&self, other: &Self) -> Ordering {
        for i in (0..5).rev() {
            if self[i] > other[i] {
                return Ordering::Greater;
            }else if self[i] < other[i] {
                return Ordering::Less;
            }
        }
        Ordering::Equal
    }
}

//-------------- From Implementations -----------------//
impl<'a> From<&'a u8> for Scalar {
    /// Performs the conversion.
    fn from(_inp: &'a u8) -> Scalar {
        let mut res = Scalar::zero();
        res[0] = *_inp as u64;
        res
    }
}

impl<'a> From<&'a u16> for Scalar {
    /// Performs the conversion.
    fn from(_inp: &'a u16) -> Scalar {
        let mut res = Scalar::zero();
        res[0] = *_inp as u64;
        res
    }
}

impl<'a> From<&'a u32> for Scalar {
    /// Performs the conversion.
    fn from(_inp: &'a u32) -> Scalar {
        let mut res = Scalar::zero();
        res[0] = *_inp as u64;
        res
    }
}

impl<'a> From<&'a u64> for Scalar {
    /// Performs the conversion.
    fn from(_inp: &'a u64) -> Scalar {
        let mut res = Scalar::zero();
        let mask = (1u64 << 52) - 1;
        res[0] = _inp & mask;
        res[1] = _inp >> 52;
        res
    }
}

impl<'a> From<&'a u128> for Scalar {
    /// Performs the conversion.
    fn from(_inp: &'a u128) -> Scalar {
        let mut res = Scalar::zero();
        let mask = (1u128 << 52) - 1;

        // Since 128 / 52 < 4 , we only need to care
        // about the first three limbs.
        res[0] = (_inp & mask) as u64;
        res[1] = ((_inp >> 52) & mask) as u64;
        res[2] = (_inp >> 104) as u64;

        res
    }
}

impl<'a> Neg for &'a Scalar {
    type Output = Scalar;
    /// Performs the negate operation over the
    /// sub-group modulo l. 
    fn neg(self) -> Scalar {
        &Scalar::zero() - &self
    }
}

impl Neg for Scalar {
    type Output = Scalar;
    /// Performs the negate operation over the
    /// sub-group modulo l.
    fn neg(self) -> Scalar {
        -&self
    }
}

impl Identity for Scalar {
    /// Returns the `Identity` element for `Scalar`
    /// which equals `1 (mod l)`. 
    fn identity() -> Scalar {
        Scalar::one()
    }
}

impl<'a, 'b> Add<&'b Scalar> for &'a Scalar {
    type Output = Scalar;
    /// Compute `a + b (mod l)`.
    fn add(self, b: &'b Scalar) -> Scalar {
        let mut sum = Scalar::zero();
        let mask = (1u64 << 52) - 1;

        // a + b
        let mut carry: u64 = 0;
        for i in 0..5 {
            carry = self.0[i] + b[i] + (carry >> 52);
            sum[i] = carry & mask;
        }
        // subtract l if the sum is >= l
        sum - constants::L
    }
} 

impl Add<Scalar> for Scalar {
    type Output = Scalar;
    /// Compute `a + b (mod l)`.
    fn add(self, b: Scalar) -> Scalar {
        &self + &b
    }
}

impl<'a, 'b> Sub<&'b Scalar> for &'a Scalar {
    type Output = Scalar;
    /// Compute `a - b (mod l)`.
    fn sub(self, b: &'b Scalar) -> Scalar {
        let mut difference = Scalar::zero();
        let mask = (1u64 << 52) - 1;

        // a - b
        let mut borrow: u64 = 0;
        // Save the wrapping_sub in borrow and add the remainder to the next limb.
        for i in 0..5 {
            // Borrow >> 63 so the Most Significant Bit of the remainder (2^64) can be carried to the next limb.
            borrow = self.0[i].wrapping_sub(b[i] + (borrow >> 63));
            difference[i] = borrow & mask;
        }

        // conditionally add `l` if the difference is negative.
        // Note that here borrow tells us the Most Signif Bit of the last limb so then we know if it's greater than `l`.
        let underflow_mask = ((borrow >> 63) ^ 1).wrapping_sub(1); // If isn't greater, we will not add it as XOR = 0.
        let mut carry: u64 = 0;
        for i in 0..5 {
            carry = (carry >> 52) + difference[i] + (constants::L[i] & underflow_mask);
            difference[i] = carry & mask;
        }

        difference
    }
}

impl Sub<Scalar> for Scalar {
    type Output = Scalar;
    /// Compute `a - b (mod l)`.
    fn sub(self, b: Scalar) -> Scalar {
        &self - &b
    }
}

impl<'a, 'b> Mul<&'a Scalar> for &'b Scalar {
    type Output = Scalar;
    /// This `Mul` implementation returns a double precision result.
    /// The result of the standard mul is stored on a [u128; 9].
    /// 
    /// Then, we apply the Montgomery Reduction function to perform
    /// the modulo and the reduction to the `Scalar` format: [u64; 5].
    fn mul(self, b: &'a Scalar) -> Scalar {
        let ab = Scalar::montgomery_reduce(&Scalar::mul_internal(self, b)); 
        Scalar::montgomery_reduce(&Scalar::mul_internal(&ab, &constants::RR))
    }
}

impl Mul<Scalar> for Scalar {
    type Output = Scalar;
    /// This `Mul` implementation returns a double precision result.
    /// The result of the standard mul is stored on a [u128; 9].
    /// 
    /// Then, we apply the Montgomery Reduction function to perform
    /// the modulo and the reduction to the `Scalar` format: [u64; 5].
    fn mul(self, b: Scalar) -> Scalar {
        &self * &b
    }
}

impl<'a> Square for &'a Scalar {
    type Output = Scalar;
    /// This `Square` implementation returns a double precision result.
    /// The result of the standard mul is stored on a [u128; 9].
    /// 
    /// Then, we apply the Montgomery Reduction function to perform
    /// the modulo and the reduction to the `Scalar` format: [u64; 5]. 
    fn square(self) -> Scalar {
        let aa = Scalar::montgomery_reduce(&Scalar::square_internal(self)); 
        Scalar::montgomery_reduce(&Scalar::mul_internal(&aa, &constants::RR))
    }
}

impl<'a> Half for &'a Scalar {
    type Output = Scalar;
    /// Give the half of the Scalar value (mod l).
    /// 
    /// This op **SHOULD ONLY** be used with even 
    /// `Scalars` otherways, can produce erroneus
    /// results.
    /// 
    /// The implementation for `Scalar` has indeed 
    /// an `assert!` statement to check this.
    #[inline]
    fn half(self) -> Scalar {
        assert!(self.is_even(), "The Scalar has to be even.");
        let mut res = self.clone();
        let mut remainder = 0u64;
        for i in (0..5).rev() {
            res[i] = res[i] + remainder;
            match(res[i] == 1, res[i].is_even()){
                (true, _) => {
                    remainder = 4503599627370496u64;
                }
                (_, false) => {
                    res[i] = res[i] - 1u64;
                    remainder = 4503599627370496u64;
                }
                (_, true) => {
                    remainder = 0;
                }
            }
            res[i] = res[i] >> 1;
        };
        res
    }
}

/// u64 * u64 = u128 inline func multiply helper.
#[inline]
fn m(x: u64, y: u64) -> u128 {
    (x as u128) * (y as u128)
}

/// u64 * u64 = u128 macro multiply helper
macro_rules! m {
    ($x:expr, $y:expr) => {
        $x as u128 * $y as u128 
    }
}


impl Scalar {

    /// Return a Scalar with value = `0`.
    pub fn zero() -> Scalar {
        Scalar([0,0,0,0,0])
    }

    /// Return a Scalar with value = `1`.
    pub fn one() -> Scalar {
        Scalar([1,0,0,0,0])
    }

    /// Return a Scalar with value = `-1 (mod l)`.
    pub fn minus_one() -> Scalar {
        Scalar([2766226127823334, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551])
    }

    /// Evaluate if a `Scalar` is even or not.
    pub fn is_even(self) -> bool {
        self.0[0].is_even()
    }

    /// Unpack a 32 byte / 256 bit Scalar into 5 52-bit limbs.
    pub fn from_bytes(bytes: &[u8; 32]) -> Scalar {
        let mut words = [0u64; 4];
        for i in 0..4 {
            for j in 0..8 {
                words[i] |= (bytes[(i * 8) + j] as u64) << (j * 8);
            } 
        }

        let mask = (1u64 << 52) - 1;
        let top_mask = (1u64 << 48) - 1;
        let mut s = Scalar::zero();

        s[0] =   words[0]                            & mask;
        // Get the 64-52 = 12 bits and add words[1] (shifting 12 to the left) on the front with `|` then apply mask.
        s[1] = ((words[0] >> 52) | (words[1] << 12)) & mask; 
        s[2] = ((words[1] >> 40) | (words[2] << 24)) & mask;
        s[3] = ((words[2] >> 28) | (words[3] << 36)) & mask;
        // Shift 16 to the right to get the 52 bits of the scalar on that limb. Then apply top_mask.
        s[4] =  (words[3] >> 16)                     & top_mask;

        s
    }

    /// Reduce a 64 byte / 512 bit scalar mod l
    pub fn from_bytes_wide(_bytes: &[u8; 64]) -> Scalar {
       // We could provide 512 bit scalar support using Montgomery Reduction. 
       // But first we need to finnish the 256-bit implementation.
       unimplemented!()
    }

    /// Pack the limbs of this `Scalar` into 32 bytes
    pub fn to_bytes(&self) -> [u8; 32] {
        let mut res = [0u8; 32];

        res[0]  =  (self.0[0] >> 0)                        as u8;
        res[1]  =  (self.0[0] >> 8)                        as u8;
        res[2]  =  (self.0[0] >> 16)                       as u8;
        res[3]  =  (self.0[0] >> 24)                       as u8;
        res[4]  =  (self.0[0] >> 32)                       as u8;
        res[5]  =  (self.0[0] >> 40)                       as u8;
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

    /// Given a `k`: u64, compute `2^k` giving the resulting result
    /// as a `Scalar`.
    /// 
    /// See that the input must be between the range => 0..250.
    /// 
    /// NOTE: This function implements an `assert!` statement that
    /// checks the correctness of the exponent provided as param.
    pub fn two_pow_k(exp: &u64) -> Scalar {
        
        // Check that exp has to be less than 260.
        // Note that a Scalar can be as much
        // `2^249 - 15145038707218910765482344729778085401` so we pick
        // 250 knowing that 249 will be lower than the prime of the
        // sub group.
        assert!(exp < &253u64, "Exponent can't be greater than 260");
        
        let mut res = Scalar::zero();
        match exp {
            0...51 => {
               res[0]  = 1u64 << exp;
            },
            52...103 => {
                res[1] = 1u64 << (exp - 52);
            },
            104...155 => {
                res[2] = 1u64 << (exp - 104);
            },
            156...207 => {
                res[3] = 1u64 << (exp - 156);
            },
            _ => {
                res[4] = 1u64 << (exp - 208);
            }
        }
        res
    }

    /// Compute `a * b`.
    /// Note that this is just the normal way of performing a product. 
    /// This operation returns back a double precision result stored
    /// on a `[u128; 9] in order to avoid overflowings.
    #[inline]
    pub(self) fn mul_internal(a: &Scalar, b: &Scalar) -> [u128; 9] {
        let mut res = [0u128; 9];
        
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

    #[allow(dead_code)]
    #[inline]
    /// Compute `a * b`.
    /// Note that this is just the normal way of performing a product. 
    /// This operation returns back a double precision result stored
    /// on a `[u128; 9] in order to avoid overflowings.
    pub(self) fn mul_internal_macros(a: &Scalar, b: &Scalar) -> [u128; 9] {
        let mut res = [0u128; 9];
        
        res[0] = m!(a[0],b[0]);
        res[1] = m!(a[0],b[1]) + m!(a[1],b[0]);
        res[2] = m!(a[0],b[2]) + m!(a[1],b[1]) + m!(a[2],b[0]);
        res[3] = m!(a[0],b[3]) + m!(a[1],b[2]) + m!(a[2],b[1]) + m!(a[3],b[0]);
        res[4] = m!(a[0],b[4]) + m!(a[1],b[3]) + m!(a[2],b[2]) + m!(a[3],b[1]) + m!(a[4],b[0]);
        res[5] =                                 m!(a[1],b[4]) + m!(a[2],b[3]) + m!(a[3],b[2]) + m!(a[4],b[1]);
        res[6] =                                                 m!(a[2],b[4]) + m!(a[3],b[3]) + m!(a[4],b[2]);
        res[7] =                                                                 m!(a[3],b[4]) + m!(a[4],b[3]);
        res[8] =                                                                                 m!(a[4],b[4]);

        res
    }

    /// Compute `a^2`. 
    /// 
    /// This operation returns a double precision result. 
    /// So it gives back a `[u128; 9]` with the result of the squaring.
    #[inline]
    pub(self) fn square_internal(a: &Scalar) -> [u128; 9] {
        let a_sqrt = [
            a[0]*2,
            a[1]*2,
            a[2]*2,
            a[3]*2,
        ];

        [
            m(a[0],a[0]),
            m(a_sqrt[0],a[1]),
            m(a_sqrt[0],a[2]) + m(a[1],a[1]),
            m(a_sqrt[0],a[3]) + m(a_sqrt[1],a[2]),
            m(a_sqrt[0],a[4]) + m(a_sqrt[1],a[3]) + m(a[2],a[2]),
                                m(a_sqrt[1],a[4]) + m(a_sqrt[2],a[3]),
                                                    m(a_sqrt[2],a[4]) + m(a[3],a[3]),
                                                                        m(a_sqrt[3],a[4]),
                                                                        m(a[4],a[4])
        ]
    }

    /// Give the half of the Scalar value (mod l).
    /// 
    /// This op **SHOULD NEVER** be used by the end-user
    /// since it's designed to allow some behaviours
    /// needed on certain points of algorithm implementations.
    #[inline]
    #[doc(hidden)]
    pub(crate) fn inner_half(self) -> Scalar {
        let mut res = self.clone();
        let mut remainder = 0u64;
        for i in (0..5).rev() {
            res[i] = res[i] + remainder;
            match(res[i] == 1, res[i].is_even()){
                (true, _) => {
                    remainder = 4503599627370496u64;
                }
                (_, false) => {
                    res[i] = res[i] - 1u64;
                    remainder = 4503599627370496u64;
                }
                (_, true) => {
                    remainder = 0;
                }
            }
            res[i] = res[i] >> 1;
        };
        res
    }

    /// Compute `limbs/R` (mod l), where R is the Montgomery modulus 2^260
    #[inline]
    pub(self) fn montgomery_reduce(limbs: &[u128; 9]) -> Scalar {

        #[inline]
        fn adjustment_fact(sum: u128) -> (u128, u64) {
            let p = (sum as u64).wrapping_mul(constants::LFACTOR) & ((1u64 << 52) - 1);
            ((sum + m(p,constants::L[0])) >> 52, p)
        }

        #[inline]
        fn montg_red_res(sum: u128) -> (u128, u64) {
            let w = (sum as u64) & ((1u64 << 52) - 1);
            (sum >> 52, w)
        }

        let l = &constants::L;
        
        // the first half computes the Montgomery adjustment factor n, and begins adding n*l to make limbs divisible by R
        let (carry, n0) = adjustment_fact(        limbs[0]);
        let (carry, n1) = adjustment_fact(carry + limbs[1] + m(n0,l[1]));
        let (carry, n2) = adjustment_fact(carry + limbs[2] + m(n0,l[2]) + m(n1,l[1]));
        let (carry, n3) = adjustment_fact(carry + limbs[3] + m(n0,l[3]) + m(n1,l[2]) + m(n2,l[1]));
        let (carry, n4) = adjustment_fact(carry + limbs[4] + m(n0,l[4]) + m(n1,l[3]) + m(n2,l[2]) + m(n3,l[1]));

        // limbs is divisible by R now, so we can divide by R by simply storing the upper half as the result
        let (carry, r0) = montg_red_res(carry + limbs[5]              + m(n1,l[4]) + m(n2,l[3]) + m(n3,l[2]) + m(n4,l[1]));
        let (carry, r1) = montg_red_res(carry + limbs[6]                           + m(n2,l[4]) + m(n3,l[3]) + m(n4,l[2]));
        let (carry, r2) = montg_red_res(carry + limbs[7]                                        + m(n3,l[4]) + m(n4,l[3]));
        let (carry, r3) = montg_red_res(carry + limbs[8]                                                     + m(n4,l[4]));
        let         r4 = carry as u64;

        // result may be >= r, so attempt to subtract l
        &Scalar([r0,r1,r2,r3,r4]) - l
    }

    /// Compute `(a * b) / R` (mod l), where R is the Montgomery modulus 2^260
    #[inline]
    #[allow(dead_code)]
    pub(self) fn montgomery_mul(a: &Scalar, b: &Scalar) -> Scalar {
        Scalar::montgomery_reduce(&Scalar::mul_internal(a, b))
    }

    /// Puts a Scalar into Montgomery form, i.e. computes `a*R (mod l)`
    #[inline]
    #[allow(dead_code)]
    pub(self) fn to_montgomery(&self) -> Scalar {
        Scalar::montgomery_mul(self, &constants::RR)
    }

    /// Takes a Scalar out of Montgomery form, i.e. computes `a/R (mod l)`
    #[inline]
    #[allow(dead_code)]
    pub(self) fn from_montgomery(&self) -> Scalar {
        let mut limbs = [0u128; 9];
        for i in 0..5 {
            limbs[i] = self[i] as u128;
        }
        Scalar::montgomery_reduce(&limbs)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
        
    /// `A = 182687704666362864775460604089535377456991567872`.
    pub static A: Scalar = Scalar([0, 0, 0, 2, 0]);

    /// `B = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static B: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    /// `AB = A - B = -904625697166532776746648320014998870755800986942176787613709275418060104167 (mod l)`.
    /// which is equal to: `365375409332725729550921208179070754913983135744`.
    pub static AB: Scalar = Scalar([0, 0, 0, 4, 0]);

    /// `BA = B - A = 904625697166532776746648320014998870755800986942176787613709275418060104167`.
    pub static BA: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370491, 2199023255551]);

    /// `A * AB (mod l). Result expected of the product mentioned before.
    pub static A_TIMES_AB: [u128; 9] = [0,0,0,0,0,0,0,8,0];

    /// `B * BA` computed in Sage limb by limb. (Since we don't have any other way to verify it.)
    pub static B_TIMES_BA: [u128; 9] = 
        [7652006990252481706224970522225, 
        23445622381543053554951959203660, 
        42875199347605145563220777152894, 
        63086978359456741425512297249892, 
        58465604036906492621308128018971, 
        40583457398062310210466901672404, 
        20302216644276907411437105105337, 
        19807040628557059606945202184, 
        4835703278454118652313601];

    /// A in Montgomery domain; `A_MONT = (A * R) (mod l) = 295345389055300509611653655730781949282003822754281035405592286100742960688`
    pub static A_MONT: Scalar = Scalar([946644518663728, 4368868487057990, 2524289321948647, 594442899788814, 717944870444]);

    /// `X = 1809251394333065553493296640760748560207343510400633813116524750123642650623`
    pub static X: Scalar = Scalar([4503599627370495, 4503599627370495, 4503599627370495, 4503599627370495, 4398046511103]);
    
    /// `Y = 717350576871794411262215878514291949349241575907629849852603275827191647632`.
    pub static Y: Scalar = Scalar([138340288859536, 461913478537005, 1182880083788836, 1688835920473363, 1743782656037]);

    /// `Y^2 (mod l) = 351405481126033478170820083848817267677692781462667634162278835827410585241`. 
    pub static Y_SQ: Scalar = Scalar([2359521284310681, 3823495160731511, 2863901539039406, 2131140264591444, 854219405379]);

    /// `Y/2 = 358675288435897205631107939257145974674620787953814924926301637913595823816`.
    pub static Y_HALF: Scalar = Scalar([2320969958115016, 230956739268502, 2843239855579666, 3096217773921929, 871891328018]);
    
    /// Y in Montgomery domain; `Y_MONT = (Y * R) (mod l) = 682963356548663143913382285837893622221394109239214830065314998385324548003`
    pub static Y_MONT: Scalar = Scalar([2328716356837283, 1997480944140188, 4481133454453893, 3196446152249575, 1660191953914]);
    
    /// `(X * Y)/R (mod l) = 781842614815424000988673591006250240924873016371899513350486876443913409068`
    pub static X_TIMES_Y_MONT: Scalar = Scalar([1458967730377260, 963769115966027, 34859148282403, 2124040828839810, 1900554115968]); 

    /// `X * Y (mod l) = 890263784947025690345271110799906008759402458672628420828189878638015362081`
    pub static X_TIMES_Y: Scalar = Scalar([3414372756436001, 1500062170770321, 4341044393209371, 2791496957276064, 2164111380879]); 


    //------------------ Tests ------------------//

    #[test]
    fn partial_ord_and_eq() {
        assert!(Y.is_even());
        assert!(!X.is_even());

        assert!(A_MONT < Y);
        assert!(Y < X);

        assert!(Y >= Y);
        assert!(X == X);
    }

    #[test]
    fn add_with_modulo() {
        let res = A + B;
        let zero = Scalar::zero();;
        for i in 0..5 {
            assert!(res[i] == zero[i]);
        }
    }

    #[test]
    fn sub_with_modulo() {
        let res = A - B;
        for i in 0..5 {
            assert!(res[i] == AB[i]);
        }
    }

    #[test]
    fn sub_without_modulo() {
        let res = B - A;
        for i in 0..5 {
            assert!(res[i] == BA[i]);
        }
    }

    #[test]
    fn mul_internal() {
        let easy_res = Scalar::mul_internal(&A, &AB);
        for i in 0..5 {
            assert!(easy_res[i] == A_TIMES_AB[i]);
        }

        let res = Scalar::mul_internal(&B, &BA);
        for i in 0..9 {
            assert!(res[i] == B_TIMES_BA[i]);
        }
    }

    #[test]
    fn square_internal() {
        let easy_res = Scalar::square_internal(&AB);
        let res_correct: [u128; 9] = [0,0,0,0,0,0,16,0,0];
        for i in 0..5 {
            assert!(easy_res[i] == res_correct[i]);
        }
    }

    #[test]
    fn to_montgomery_conversion() {
        let a = Scalar::to_montgomery(&A);
        for i in 0..5 {
            assert!(a[i] == A_MONT[i]);
        }
    }

    #[test]
    fn from_montgomery_conversion() {
        let y = Scalar::from_montgomery(&Y_MONT);
        for i in 0..5 {
            assert!(y[i] == Y[i]);
        }
    }

    #[test]
    fn scalar_mul() {
        let res = &X * &Y;
        for i in 0..5 {
            assert!(res[i] == X_TIMES_Y[i]);
        }
    }

    #[test]
    fn mul_by_identity() {
        let res = &Y * &Scalar::identity();

        println!("{:?}", res);
        for i in 0..5 {
            assert!(res[i] == Y[i]);
        }
    }

    #[test]
    fn mul_by_zero() {
        let res = &Y * &Scalar::zero();
        for i in 0..5 {
            assert!(res[i] == Scalar::zero()[i]);
        }
    }

    #[test]
    fn montgomery_mul() {
        let res = Scalar::montgomery_mul(&X, &Y);
        for i in 0..5 {
            assert!(res[i] == X_TIMES_Y_MONT[i]);
        }
    }

    #[test]
    fn square() {
        let res = &Y.square();

        for i in 0..5 {
            assert!(res[i] == Y_SQ[i]);
        }
    }

    #[test]
    fn square_zero_and_identity() {
        let zero = &Scalar::zero().square();
        let one = &Scalar::identity().square();

        for i in 0..5 {
            assert!(zero[i] == Scalar::zero()[i]);
            assert!(one[i] == Scalar::one()[i]);
        }
    }

    #[test]
    fn half() {
        let res = &Y.half();
        for i in 0..5 {
            assert!(res[i] == Y_HALF[i]);
        }

        let a_half = Scalar([0, 0, 0, 1, 0]);
        let a_half_half = Scalar([0, 0, 2251799813685248, 0, 0]);

        for i in 0..5 {
            assert!(a_half[i] == A.half()[i]);
            assert!(a_half_half[i] == A.half().half()[i]);
        }
    }

    #[test]
    fn even_scalar() {
        assert!(Y.is_even());
        assert!(!X.is_even());
        assert!(Scalar::zero().is_even());
    }
}