//! Arithmetic mod \\(2\^{249} - 15145038707218910765482344729778085401\\)
//! with five \\(52\\)-bit unsigned limbs.
//!
//! \\(51\\)-bit limbs would cover the desired bit range (\\(253\\)
//! bits), but isn't large enough to reduce a \\(512\\)-bit number with
//! Montgomery multiplication, so \\(52\\) bits is used instead.  To see
//! that this is safe for intermediate results, note that the largest
//! limb in a \\(5\times 5\\) product of \\(52\\)-bit limbs will be
//!
//! ```text
//! (0xfffffffffffff^2) * 5 = 0x4ffffffffffff60000000000005 (107 bits).
//! ```

use core::fmt::Debug;
use core::ops::{Index, IndexMut};
use core::ops::Add;
use core::ops::Sub;
use crate::backend::u64::constants;

/// The `Scalar` struct represents an element in
/// \\(\mathbb Z / \ell \mathbb Z\\) as 5 \\(52\\)-bit limbs.
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

impl<'b> Add<&'b Scalar> for Scalar {
    type Output = Scalar;
    /// Compute `a + b` (mod l)
    fn add(self, b: &'b Scalar) -> Scalar {
        let mut sum = Scalar::zero();
        let mask = (1u64 << 52) - 1;

        // a + b
        let mut carry: u64 = 0;
        for i in 0..5 {
            carry = self.0[i] + b[i] + (carry >> 52);
            sum[i] = carry & mask;
        }
        // subtract r if the sum is >= r
        sum.sub(&constants::R)
    }
} 

impl<'b> Sub<&'b Scalar> for Scalar {
    type Output = Scalar;
    /// Compute `a - b` (mod r)
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

        // conditionally add `r` if the difference is negative.
        // Note that here borrow tells us the Most Signif Bit of the last limb so then we know if it's greater than `r`.
        let underflow_mask = ((borrow >> 63) ^ 1).wrapping_sub(1); // If isn't greater, we will not add it as XOR = 0.
        let mut carry: u64 = 0;
        for i in 0..5 {
            carry = (carry >> 52) + difference[i] + (constants::R[i] & underflow_mask);
            difference[i] = carry & mask;
        }

        difference
    }
}

/// u64 * u64 = u128 macro multiply helper

macro_rules! m {
    ($x:expr, $y:expr) => {
        ($x as u128 * $y as u128) 
    }
}

/// u64 * u64 = u128 inline func multiply helper
#[inline]
fn m(x: u64, y: u64) -> u128 {
    (x as u128) * (y as u128)
}

impl Scalar {
    /// Return the zero Scalar
    pub fn zero() -> Scalar {
        Scalar([0,0,0,0,0])
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

        s[ 0] =   words[0]                            & mask;
        // Get the 64-52 = 12 bits and add words[1] (shifting 12 to the left) on the front with `|` then apply mask.
        s[ 1] = ((words[0] >> 52) | (words[1] << 12)) & mask; 
        s[ 2] = ((words[1] >> 40) | (words[2] << 24)) & mask;
        s[ 3] = ((words[2] >> 28) | (words[3] << 36)) & mask;
        // Shift 16 to the right to get the 52 bits of the scalar on that limb. Then apply top_mask.
        s[ 4] =  (words[3] >> 16)                     & top_mask;

        s
    }

    /// Reduce a 64 byte / 512 bit scalar mod l
    pub fn from_bytes_wide(bytes: &[u8; 64]) -> Scalar {
       // We could provide 512 bit scalar support using Montgomery Reduction. 
       //But first we need to finnish the 256-bit implementation.
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

    /// Compute `a * b`
    #[inline]
    pub fn mul_internal(a: &Scalar, b: &Scalar) -> [u128; 9] {
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

    /// Compute `a * b` with macros
    #[inline]
    pub fn mul_internal_macros(a: &Scalar, b: &Scalar) -> [u128; 9] {
        let mut res = [0u128; 9];
        // Note that this is just the normal way of performing a product.
        // We need to store the results on u128 as otherwise we'll end
        // up having overflowings.
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

    /// Compute `a^2`
    #[inline]
    pub fn square_internal(a: &Scalar) -> [u128; 9] {
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

    /// Compute `limbs/R` (mod l), where R is the Montgomery modulus 2^260
    #[inline]
    pub (crate) fn montgomery_reduce(limbs: &[u128; 9]) -> Scalar {

        #[inline]
        fn carry_and_res(sum: u128) -> (u128, u64) {
            let p = (sum as u64).wrapping_mul(constants::LFACTOR) & ((1u64 << 52) - 1);
            ((sum + m(p,constants::R[0])) >> 52, p)
        }

        #[inline]
        fn mont_res(sum: u128) -> (u128, u64) {
            let w = (sum as u64) & ((1u64 << 52) - 1);
            (sum >> 52, w)
        }

        let r = &constants::R;

        // the first half computes the Montgomery adjustment factor n, and begins adding n*l to make limbs divisible by R
        let (carry, n0) = carry_and_res(        limbs[0]);
        let (carry, n1) = carry_and_res(carry + limbs[1] + m(n0,r[1]));
        let (carry, n2) = carry_and_res(carry + limbs[2] + m(n0,r[2]) + m(n1,r[1]));
        let (carry, n3) = carry_and_res(carry + limbs[3] + m(n0,r[3]) + m(n1,r[2]) + m(n2,r[1]));
        let (carry, n4) = carry_and_res(carry + limbs[4] + m(n0,r[4]) + m(n1,r[3]) + m(n2,r[2]) + m(n3,r[1]));

        // limbs is divisible by R now, so we can divide by R by simply storing the upper half as the result
        let (carry, r0) = mont_res(carry + limbs[5]              + m(n1,r[4]) + m(n2,r[3]) + m(n3,r[2]) + m(n4,r[1]));
        let (carry, r1) = mont_res(carry + limbs[6]                           + m(n2,r[4]) + m(n3,r[3]) + m(n4,r[2]));
        let (carry, r2) = mont_res(carry + limbs[7]                                        + m(n3,r[4]) + m(n4,r[3]));
        let (carry, r3) = mont_res(carry + limbs[8]                                                     + m(n4,r[4]));
        let         r4 = carry as u64;

        // result may be >= l, so attempt to subtract l
        Scalar::sub(Scalar([r0,r1,r2,r3,r4]), r)
    }

    /// Compute `a * b` (mod l)
    #[inline]
    pub fn mul(a: &Scalar, b: &Scalar) -> Scalar {
        let ab = Scalar::montgomery_reduce(&Scalar::mul_internal(a, b)); 
        Scalar::montgomery_reduce(&Scalar::mul_internal(&ab, &constants::RR))
    }

    /// Compute `(a * b) / R` (mod l), where R is the Montgomery modulus 2^260
    #[inline]
    pub fn montgomery_mul(a: &Scalar, b: &Scalar) -> Scalar {
        Scalar::montgomery_reduce(&Scalar::mul_internal(a, b))
    }

    /// Puts a Scalar into Montgomery form, i.e. computes `a*R (mod l)`
    #[inline(never)]
    pub fn to_montgomery(&self) -> Scalar {
        Scalar::montgomery_mul(self, &constants::RR)
    }

    /// Takes a Scalar out of Montgomery form, i.e. computes `a/R (mod l)`
    #[inline(never)]
    pub fn from_montgomery(&self) -> Scalar {
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
    /// Computed A and B, some Scalars defined over the modulo `r = 2^{249} - 15145038707218910765482344729778085401`.
    /// A = 182687704666362864775460604089535377456991567872
    pub static A: Scalar = Scalar([0, 0, 0, 2, 0]);

    /// B = 904625697166532776746648320197686575422163851717637391703244652875051672039
    pub static B: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    /// AB = A - B = -904625697166532776746648320014998870755800986942176787613709275418060104167 (mod r)
    /// which is equal to: 365375409332725729550921208179070754913983135744
    pub static AB: Scalar = Scalar([0, 0, 0, 4, 0]);

    /// BA = B - A = 904625697166532776746648320014998870755800986942176787613709275418060104167
    pub static BA: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370491, 2199023255551]);

    /// A * AB. Result expected of the product mentioned before.
    pub static A_TIMES_AB: [u128; 9] = [0,0,0,0,0,0,0,8,0];

    /// B * BA computed in Sage limb by limb. (Since we don't have any other way to verify it.)
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


    /// Montgomery x
    pub static X_MONT: Scalar = Scalar([3967855409315932, 4468715158741414, 1532785699107429, 3431559202417771, 864832837719]);

    pub static X: Scalar = Scalar([4503599627370495, 4503599627370495, 4503599627370495, 4503599627370495, 4398046511103]);
    /// Y
    pub static Y: Scalar = Scalar([3224898173688058, 3370928136179116, 1182880079308587, 1688835920473363, 14937922189349]);

    /// `(x * y)/R (mod l)`
    pub static X_TIMES_Y_MONT: Scalar = Scalar([2506269041400604, 1231103222380683, 3428697373767643, 1095442560240545, 1294911592878]); 

    /// `x * y (mod l)`
    pub static X_TIMES_Y: Scalar = Scalar([3414372756436001, 1500062170770321, 4341044393209371, 2791496957276064, 2164111380879]); 


    #[test]
    fn add_with_modulo() {
        let res = A.add(&B);
        let zero = Scalar::zero();;
        for i in 0..5 {
            assert!(res[i] == zero[i]);
        }
    }

    #[test]
    fn sub_with_modulo() {
        let res = A.sub(&B);
        for i in 0..5 {
            assert!(res[i] == AB[i]);
        }
    }

    #[test]
    fn sub_without_modulo() {
        let res = B.sub(&A);
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

    // Needs to be reviewed
    #[test]
    fn to_montgomery_conversion() {
        let mont_x = Scalar::to_montgomery(&X);
        for i in 0..5 {
            assert!(mont_x[i] == X_MONT[i]);
        }
    }

    // Need to review and compute correct values
    #[test]
    fn from_montgomery_conversion() {
        let x = Scalar::from_montgomery(&X_MONT);
        for i in 0..5 {
            assert!(x[i] == X[i]);
        }
    }

    // Need to review the values and the algorithm
    #[test]
    fn full_montgomery_mul() {
        let res = Scalar::mul(&X_MONT, &Y);
        for i in 0..5 {
            assert!(res[i] == X_TIMES_Y[i]);
        }
    }

    // Need to review the values and the algorithm
    #[test]
    fn montgomery_mul() {
        let res = Scalar::montgomery_mul(&X_MONT, &Scalar::to_montgomery(&Y));
        for i in 0..5 {
            assert!( res[i] == X_TIMES_Y_MONT[i]);
        }
    }
}