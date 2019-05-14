//! Field arithmetic modulo `2^252 + 27742317777372353535851937790883648493` 
//! using 64-bit limbs with 128-bit products

use core::fmt::Debug;
use core::ops::Neg;
use core::ops::{Index, IndexMut};
use core::ops::{Add, AddAssign};
use core::ops::{Sub, SubAssign};
use core::ops::{Mul, MulAssign};
use crate::backend::u64::constants;


/// A `FieldElement` represents an element into the field 
/// `2^252 + 27742317777372353535851937790883648493`
/// 
/// In the 64-bit backend implementation, the `FieldElement is 
/// represented in radix `2^52`

#[derive(Copy, Clone)]
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
        // subtract r if the sum is >= l
        sum.sub(&constants::FIELD_L)
    }
}

impl<'a, 'b> Sub<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn sub(self, b: &'b FieldElement) -> FieldElement {
    unimplemented!()
    }    
}



impl<'a> Neg for &'a FieldElement {
    type Output = FieldElement;
    fn neg(self) -> FieldElement {
        let mut output = *self;
        output.negate();
        output
    }
}

impl<'b> MulAssign<&'b FieldElement> for FieldElement {
    fn mul_assign(&mut self, _rhs: &'b FieldElement) {
        let result = (self as &FieldElement) * _rhs;
        self.0 = result.0;
    }
}

impl<'a, 'b> Mul<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn mul(self, _rhs: &'b FieldElement) -> FieldElement {
        unimplemented!()
    }
}


impl FieldElement {
    /// Invert the sign of this field element
    pub fn negate(&mut self) {
        unimplemented!()
    }

    /// Construct zero.
    pub fn zero() -> FieldElement {
        FieldElement([ 0, 0, 0, 0, 0 ])
    }

    /// Construct one.
    pub fn one() -> FieldElement {
        FieldElement([ 1, 0, 0, 0, 0 ])
    }

    /// Construct -1.
    pub fn minus_one() -> FieldElement {
        FieldElement([ 671914833335276, 1077929209154306, 5471207, 0, 281474976710656 ])
    }

    /// Given 64-bit input limbs, reduce to enforce the bound 2^(51 + epsilon).
    #[inline(always)]
    fn reduce(mut limbs: [u64; 5]) -> FieldElement {
        const LOW_51_BIT_MASK: u64 = (1u64 << 51) - 1;

        // Since the input limbs are bounded by 2^64, the biggest
        // carry-out is bounded by 2^13.
        //
        // The biggest carry-in is c4 * 19, resulting in
        //
        // 2^51 + 19*2^13 < 2^51.0000000001
        //
        // Because we don't need to canonicalize, only to reduce the
        // limb sizes, it's OK to do a "weak reduction", where we
        // compute the carry-outs in parallel.

        let c0 = limbs[0] >> 51;
        let c1 = limbs[1] >> 51;
        let c2 = limbs[2] >> 51;
        let c3 = limbs[3] >> 51;
        let c4 = limbs[4] >> 51;

        limbs[0] &= LOW_51_BIT_MASK;    
        limbs[1] &= LOW_51_BIT_MASK;
        limbs[2] &= LOW_51_BIT_MASK;
        limbs[3] &= LOW_51_BIT_MASK;
        limbs[4] &= LOW_51_BIT_MASK;

        // @TODO: Figure out how to build c0 with our -27742317777372353535851937790883648493 mod p
        unimplemented!()
    }

    /// Load a `FieldElement` from the low 253b   bits of a 256-bit
    /// input.
    pub fn from_bytes(bytes: &[u8;32]) -> Self {
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

        let low_51_bit_mask = (1u64 << 51) - 1;
        
        FieldElement(
        // load bits [  0, 64), no shift
        [  load8(&bytes[ 0..])        & low_51_bit_mask
        // load bits [ 48,112), shift to [ 51,112)
        , (load8(&bytes[ 6..]) >>  3) & low_51_bit_mask
        // load bits [ 96,160), shift to [102,160)
        , (load8(&bytes[12..]) >>  6) & low_51_bit_mask
        // load bits [152,216), shift to [153,216)
        , (load8(&bytes[19..]) >>  1) & low_51_bit_mask
        // load bits [192,256), shift to [204,112)
        , (load8(&bytes[24..]) >> 12) & low_51_bit_mask
        ])
    }

    /// Serialize this `FieldElement` to a 32-byte array.  The
    /// encoding is canonical.
    pub fn to_bytes(self) -> [u8; 32] {

        // @TODO: Get the part needed to implement, replace correct values and update it.
        // Wee need to prepare the limbs and reduce them before we encode them as bytes.


        let limbs = self.0;
        // Now arrange the bits of the limbs.
        let mut s = [0u8;32];
        s[ 0] =   limbs[0]        as u8;
        s[ 1] =  (limbs[0] >>  8) as u8;
        s[ 2] =  (limbs[0] >> 16) as u8;
        s[ 3] =  (limbs[0] >> 24) as u8;
        s[ 4] =  (limbs[0] >> 32) as u8;
        s[ 5] =  (limbs[0] >> 40) as u8;
        s[ 6] = ((limbs[0] >> 48) | (limbs[1] << 3)) as u8;
        s[ 7] =  (limbs[1] >>  5) as u8;
        s[ 8] =  (limbs[1] >> 13) as u8;
        s[ 9] =  (limbs[1] >> 21) as u8;
        s[10] =  (limbs[1] >> 29) as u8;
        s[11] =  (limbs[1] >> 37) as u8;
        s[12] = ((limbs[1] >> 45) | (limbs[2] << 6)) as u8;
        s[13] =  (limbs[2] >>  2) as u8;
        s[14] =  (limbs[2] >> 10) as u8;
        s[15] =  (limbs[2] >> 18) as u8;
        s[16] =  (limbs[2] >> 26) as u8;
        s[17] =  (limbs[2] >> 34) as u8;
        s[18] =  (limbs[2] >> 42) as u8;
        s[19] = ((limbs[2] >> 50) | (limbs[3] << 1)) as u8;
        s[20] =  (limbs[3] >>  7) as u8;
        s[21] =  (limbs[3] >> 15) as u8;
        s[22] =  (limbs[3] >> 23) as u8;
        s[23] =  (limbs[3] >> 31) as u8;
        s[24] =  (limbs[3] >> 39) as u8;
        s[25] = ((limbs[3] >> 47) | (limbs[4] << 4)) as u8;
        s[26] =  (limbs[4] >>  4) as u8;
        s[27] =  (limbs[4] >> 12) as u8;
        s[28] =  (limbs[4] >> 20) as u8;
        s[29] =  (limbs[4] >> 28) as u8;
        s[30] =  (limbs[4] >> 36) as u8;
        s[31] =  (limbs[4] >> 44) as u8;

        // High bit should be zero.
        debug_assert!((s[31] & 0b1000_0000u8) == 0u8);

        s
    }
}

