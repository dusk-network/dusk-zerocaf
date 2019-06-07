/// Edwards Points support and encoding/decoding proces.

use crate::field::FieldElement;
use crate::scalar::Scalar;
use crate::montgomery::MontgomeryPoint;
use crate::constants;
use crate::traits::*;


use subtle::Choice;
use subtle::ConstantTimeEq;

use std::default::Default;
use std::ops::{Add, Sub, Mul, Neg};
use std::fmt::Debug;

/// The first 255 bits of a `CompressedEdwardsY` represent the
/// \\(y\\)-coordinate.  The high bit of the 32nd byte gives the sign of \\(x\\).
#[derive(Copy, Clone, Eq, PartialEq)]
pub struct CompressedEdwardsY(pub [u8; 32]);

impl ConstantTimeEq for CompressedEdwardsY {
    fn ct_eq(&self, other: &CompressedEdwardsY) -> Choice {
        self.to_bytes().ct_eq(&other.to_bytes())
    }
}

impl Debug for CompressedEdwardsY {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "CompressedEdwardsY: {:?}", self.to_bytes())
    }
}

impl CompressedEdwardsY {

    /// Return the `CompressedEdwardsY` as an array of bytes (it's cannonical state).
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }

    /// Attempt to decompress to an `EdwardsPoint`.
    ///
    /// Returns `Err` if the input is not the \\(y\\)-coordinate of a
    /// curve point.
    pub fn decompress(&self) -> Option<EdwardsPoint> {
        unimplemented!();
    }
}

impl Identity for CompressedEdwardsY {
    fn identity() -> CompressedEdwardsY {
        CompressedEdwardsY([1, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0])
    }
}

impl Default for CompressedEdwardsY {
    fn default() -> CompressedEdwardsY {
        CompressedEdwardsY::identity()
    }
}

impl CompressedEdwardsY {
    /// Construct a `CompressedEdwardsY` from a slice of bytes.
    pub fn from_slice(bytes: &[u8]) -> CompressedEdwardsY {
        let mut tmp = [0u8; 32];

        tmp.copy_from_slice(bytes);

        CompressedEdwardsY(tmp)
    }
}


/// An `EdwardsPoint` represents a point on the Edwards form of Doppio Curve.
#[derive(Copy, Clone)]
pub struct EdwardsPoint {
    pub(crate) X: FieldElement,
    pub(crate) Y: FieldElement,
    pub(crate) Z: FieldElement,
    pub(crate) T: FieldElement,
}

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
    /// Negates an `EdwardsPoint` giving it as a result
    fn neg(self) -> EdwardsPoint {
       unimplemented!()
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
    /// Cost: 9M + 1*a + 7add.
    /// Cost: 9M + 1*a + 6add dependent upon the first point.
    /// This implementation is speciffic for curves with `a = -1` as Doppio is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson, http://eprint.iacr.org/2008/522, Section 3.1.
    #[inline]
    fn add(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        let A: FieldElement = &self.X * &other.X;
        let B: FieldElement = &self.Y * &other.Y;
        let C: FieldElement = &self.Z * &other.T;
        let D: FieldElement = &self.T * &other.Z;
        let E: FieldElement = &D + &C;
        let F: FieldElement = &(&(&(&self.X - &self.Y) * &(&other.X - &other.Y)) + &A) + &B;
        let G: FieldElement = &(&B + &constants::EDWARDS_A) * &A;
        let H: FieldElement = &D - &C;

        EdwardsPoint{
            X: &E * &F,
            Y: &G * &H,
            Z: &F * &G,
            T: &E * &H
        }
    }
}

impl<'a, 'b> Sub<&'b EdwardsPoint> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Substract two EdwardsPoints and give the resulting `EdwardsPoint`
    fn sub(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        unimplemented!()
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
    fn mul(self, point: &'b EdwardsPoint) -> EdwardsPoint {
        unimplemented!()
    }
}



impl EdwardsPoint {
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
    pub(crate) fn mul_by_pow_2(&self, k: u32) -> EdwardsPoint {
        unimplemented!()
    }
}

pub mod tests {
    use super::*;
    use constants::*;


}
