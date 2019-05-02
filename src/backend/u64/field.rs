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

    /// Construct -1 (mod l).
    pub fn minus_one() -> FieldElement {
        FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416])
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



pub mod tests {

    use crate::backend::u64::field::FieldElement;
    use crate::backend::u64::constants;

    /// Bytes representation of `-1 (mod l) = 7237005577332262213973186563042994240857116359379907606001950938285454250988`
    pub(crate) static MINUS_ONE_BYTES: [u8; 32] = [236, 211, 245, 92, 26, 99, 18, 88, 214, 156, 247, 162, 222, 249, 222, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16];

    #[test]
    #[ignore]
        fn it_adds_correctly() {
        let a = FieldElement([929955233495203, 466365720129213, 1662059464998953, 2033849074728123, 1442794654840575]);
        let a_2 = FieldElement([1859910466990425, 932731440258426, 1072319116312658, 1815898335770999, 633789495995903]);
        let res = &a + &a;
        assert_eq!(res.to_bytes(), a_2.to_bytes());
    }

    #[test]
    fn from_bytes_conversion() {
        let num = FieldElement::from_bytes(&MINUS_ONE_BYTES);
        for i in 0..5 {
            assert!(num[i] == FieldElement::minus_one()[i]);
        }
    }

    #[test]
    #[ignore]
    fn to_bytes_conversion() {
        let bytes = FieldElement::minus_one().to_bytes();
        println!("{:?}", bytes);
        for i in 0..32 {
            assert!(bytes[i] == MINUS_ONE_BYTES[i]);
        }
    }

}