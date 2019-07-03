//! Edwards Point operation implementations and definitions.
//! Encoding/decoding processes implementation 
//! and support for all kind of interactions with them.

use crate::field::FieldElement;
use crate::scalar::Scalar;
use crate::montgomery::MontgomeryPoint;
use crate::constants;
use crate::traits::*;


use subtle::Choice;
use subtle::ConstantTimeEq;

use std::default::Default;
use std::fmt::Debug;

use core::ops::{Index, IndexMut};
use std::ops::{Add, Sub, Mul, Neg};


/// The first 255 bits of a `CompressedEdwardsY` represent the
/// (y)-coordinate.  The high bit of the 32nd byte gives the sign of (x).
#[derive(Copy, Clone, Eq, PartialEq)]
pub struct CompressedEdwardsY(pub [u8; 32]);

impl ConstantTimeEq for CompressedEdwardsY {
    fn ct_eq(&self, other: &CompressedEdwardsY) -> Choice {
        self.to_bytes().ct_eq(&other.to_bytes())
    }
}

impl Index<usize> for CompressedEdwardsY {
    type Output = u8;
    fn index(&self, _index: usize) -> &u8 {
        &(self.0[_index])
    }
}

impl IndexMut<usize> for CompressedEdwardsY {
    fn index_mut(&mut self, _index: usize) -> &mut u8 {
        &mut (self.0[_index])
    }
}

impl Debug for CompressedEdwardsY {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "CompressedEdwardsY: {:?}", self.to_bytes())
    }
}

impl Default for CompressedEdwardsY {
    /// Returns the identity for `CompressedEdwardsY` point.
    fn default() -> CompressedEdwardsY {
        CompressedEdwardsY::identity()
    }
}

impl<'a> Neg for &'a CompressedEdwardsY {
    type Output = CompressedEdwardsY;
    /// Negates an `CompressedEdwardsY` by decompressing
    /// it, negating over Twisted Edwards Extended 
    /// Projective Coordinates and compressing it back.
    fn neg(self) -> CompressedEdwardsY {
        (-&self.decompress().unwrap()).compress()
    }
}

impl Neg for CompressedEdwardsY {
    type Output = CompressedEdwardsY;
    /// Negates an `CompressedEdwardsY` by decompressing
    /// it, negating over Twisted Edwards Extended 
    /// Projective Coordinates and compressing it back.
    fn neg(self) -> CompressedEdwardsY {
        -& self
    }
}

impl Identity for CompressedEdwardsY {
    /// Returns the `CompressedEdwardsY` identity point value 
    /// that corresponds to `1 (mod l)`
    /// with the sign bit setted to `0`.
    fn identity() -> CompressedEdwardsY {
        CompressedEdwardsY([1, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0])
    }
}

impl CompressedEdwardsY {
    /// Construct a `CompressedEdwardsY` from a slice of bytes.
    pub fn from_slice(bytes: &[u8]) -> CompressedEdwardsY {
        let mut tmp = [0u8; 32];

        tmp.copy_from_slice(bytes);

        CompressedEdwardsY(tmp)
    }

    /// Return the `CompressedEdwardsY` as an array of bytes (it's cannonical state).
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }

    /// Attempt to decompress to an `EdwardsPoint`.
    ///
    /// Returns `Err` if the input is not the Y-coordinate of a
    /// curve point.
    pub fn decompress(&self) -> Option<EdwardsPoint> {
        unimplemented!()
    }
}



/// An `EdwardsPoint` represents a point on the Doppio Curve expressed
/// over the Twisted Edwards Extended Coordinates.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct EdwardsPoint {
    pub X: FieldElement,
    pub Y: FieldElement,
    pub Z: FieldElement,
    pub T: FieldElement,
}
/*
impl ConstantTimeEq for EdwardsPoint {
    fn ct_eq(&self, other: &EdwardsPoint) -> Choice {
        self.compress().ct_eq(&other.compress())
    }
}

impl PartialEq for EdwardsPoint {
    fn eq(&self, other: &EdwardsPoint) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}
*/
impl Default for EdwardsPoint {
    /// Returns the default EdwardsPoint Extended Coordinates: (0, 1, 1, 0). 
    fn default() -> EdwardsPoint {
        EdwardsPoint::identity()
    }
}

impl Identity for EdwardsPoint {
    /// Returns the Edwards Point identity value = `(0, 1, 1, 0)`.
    fn identity() -> EdwardsPoint {
        EdwardsPoint {
            X: FieldElement::zero(),
            Y: FieldElement::one(),
            Z: FieldElement::one(),
            T: FieldElement::zero()
        }
    }
}

impl<'a> Neg for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Negates an `EdwardsPoint` giving it as a result.
    /// Since the negative of a point is (-X:Y:Z:-T), it
    /// gives as a result: `(-X, Y, Z, -T)`.
    fn neg(self) -> EdwardsPoint {
       EdwardsPoint{
           X: -&self.X,
           Y:   self.Y,
           Z:   self.Z,
           T: -&self.T,
       }
    }
}

impl Neg for EdwardsPoint {
    type Output = EdwardsPoint;
    /// Negates an `EdwardsPoint` giving it as a result
    fn neg(self) -> EdwardsPoint {
        -&self
    }
}

impl<'a, 'b> Add<&'b EdwardsPoint> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Add two EdwardsPoints and give the resulting `EdwardsPoint`.
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// [Source: 2008 Hisil–Wong–Carter–Dawson], 
    /// (http://eprint.iacr.org/2008/522), Section 3.1.
    /// 
    /// Trick 2k' is used to reduce 1D and 1M. `2d' = k = 2*(-a/EDWARDS_D)`.
    #[inline]
    fn add(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        let k = constants::EDWARDS_2_D_PRIME;
        let two = FieldElement::from(&2u8);

        let A = &(&self.Y - &self.X) * &(&other.Y - &other.X);
        let B = &(&self.Y + &self.X) * &(&other.Y + &other.X);
        let C = &k * &(&self.T * &other.T);
        let D = &two * &(&self.Z * &other.Z);
        
        EdwardsPoint {
            X: &(&B - &A) * &(&D - &C),
            Y: &(&D + &C) * &(&B + &A),
            Z: &(&D - &C) * &(&D + &C),
            T: &(&B - &A) * &(&B + &A),
        }
    }
}

impl<'a, 'b> Sub<&'b EdwardsPoint> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Substract two EdwardsPoints and give the resulting `EdwardsPoint`
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson, 
    /// http://eprint.iacr.org/2008/522, Section 3.1.
    /// 
    /// Trick 2k' is used to reduce 1D and 1M. `2d' = k = 2*(-a/EDWARDS_D)`.
    /// 
    /// The only thing we do is to negate the second `EdwardsPoint`
    /// and add it following the same addition algorithm.
    fn sub(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        let other_neg = -other;
        let k = constants::EDWARDS_2_D_PRIME;
        let two = FieldElement::from(&2u8);

        let A = &(&self.Y - &self.X) * &(&other_neg.Y - &other_neg.X);
        let B = &(&self.Y + &self.X) * &(&other_neg.Y + &other_neg.X);
        let C = &k * &(&self.T * &other_neg.T);
        let D = &two * &(&self.Z * &other_neg.Z);
        
        EdwardsPoint {
            X: &(&B - &A) * &(&D - &C),
            Y: &(&D + &C) * &(&B + &A),
            Z: &(&D - &C) * &(&D + &C),
            T: &(&B - &A) * &(&B + &A),
        }
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Scalar multiplication: compute `scalar * self`.
    fn mul(self, scalar: &'b Scalar) -> EdwardsPoint {
        unimplemented!()
    }
}

impl<'a, 'b> Mul<&'b EdwardsPoint> for &'a Scalar {
    type Output = EdwardsPoint;
    /// Scalar multiplication: compute `scalar * self`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    fn mul(self, point: &'b EdwardsPoint) -> EdwardsPoint {
        unimplemented!()
    }
}



impl EdwardsPoint {
    /// Double the given point following:
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson, 
    /// http://eprint.iacr.org/2008/522, Section 3.1.
    /// Cost: 4M+ 4S+ 1D
    pub fn double(&self) -> EdwardsPoint {
        let two = FieldElement::from(&2u8);
        let a = FieldElement::minus_one();

        let A = self.X.square();
        let B = self.Y.square();
        let C = &two * &self.Z.square();
        let D = &a * &A;
        let E = &((&(&self.X + &self.Y).square()) - &A) - &B;
        let G = &D + &B;
        let F = &G - &C;
        let H = &D - &B;

        EdwardsPoint {
            X: &E * &F,
            Y: &G * &H,
            Z: &F * &G,
            T: &E * &H,
        }

    }

    /// Convert this `EdwardsPoint` on the Edwards model to the
    /// corresponding `MontgomeryPoint` on the Montgomery model.
    pub fn to_montgomery(&self) -> MontgomeryPoint {
       unimplemented!()
    }

    /// Compress this point to `CompressedEdwardsY` format.
    pub fn compress(&self) -> CompressedEdwardsY {
        unimplemented!()
    }
    

    /// Multiply by the cofactor: return (8 P).
    pub fn mul_by_cofactor(&self) -> EdwardsPoint {
        unimplemented!()
    }

    /// Compute ([2^k] P)
    pub fn mul_by_pow_2(&self, k: u32) -> EdwardsPoint {
        unimplemented!()
    }
}

/// Module used for tesing `EdwardsPoint` operations and implementations.
/// Also used to check the correctness of the transformation functions.
pub mod tests {
    use super::*;
    use constants::*;

    #[test]
    fn edwards_extended_coords_neg() {

        let inv_a: EdwardsPoint = EdwardsPoint{
           X: FieldElement::minus_one(),
           Y: FieldElement::zero(),
           Z: FieldElement::zero(),
           T: FieldElement::minus_one(),
        };

        let a: EdwardsPoint = EdwardsPoint{
           X: FieldElement::one(),
           Y: FieldElement::zero(),
           Z: FieldElement::zero(),
           T: FieldElement::one(),
        };

        let res = -a;
        assert!(res == inv_a);
    }

    #[test]
    fn edwards_extended_coords_neg_identity() {
        let res = - &EdwardsPoint::identity();

        assert!(res == EdwardsPoint::identity())
    }
}
