#![allow(non_snake_case)]
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
#[derive(Copy, Clone, Eq, PartialEq)]
pub struct EdwardsPoint {
    pub X: FieldElement,
    pub Y: FieldElement,
    pub Z: FieldElement,
    pub T: FieldElement,
}

impl Debug for EdwardsPoint {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "
        EdwardsPoint {{
            X: {:?},
            Y: {:?},
            Z: {:?},
            T: {:?}
        }};", self.X, self.Y, self.Z, self.T)
    }
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

impl Eq for EdwardsPoint {}
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
    #[inline]
    fn add(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        let A = &self.X * &other.X;
        let B = &self.Y * &other.Y;
        let C = &constants::EDWARDS_D * &(&self.T * &other.T);
        let D = &self.Z * &other.Z;
        let E = &(&(&(&self.X + &self.Y) * &(&other.X + &other.Y)) - &A) - &B;
        let F = &D - &C;
        let G = &D + &C;
        let H = &B - &(&constants::EDWARDS_A * &A);

        EdwardsPoint {
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
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson, 
    /// http://eprint.iacr.org/2008/522, Section 3.1.
    /// 
    /// The only thing we do is to negate the second `EdwardsPoint`
    /// and add it following the same addition algorithm.
    fn sub(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        let other_neg = -other;
        
        let A = &self.X * &other_neg.X;
        let B = &self.Y * &other_neg.Y;
        let C = &constants::EDWARDS_D * &(&self.T * &other_neg.T);
        let D = &self.Z * &other_neg.Z;
        let E = &(&(&(&self.X + &self.Y) * &(&other_neg.X + &other_neg.Y)) - &A) - &B;
        let F = &D - &C;
        let G = &D + &C;
        let H = &B - &(&constants::EDWARDS_A * &A);

        EdwardsPoint {
            X: &E * &F,
            Y: &G * &H,
            Z: &F * &G,
            T: &E * &H
        }
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Scalar multiplication: compute `scalar * self`.
    fn mul(self, scalar: &'b Scalar) -> EdwardsPoint {
        self.double_and_add(scalar)
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
        point.double_and_add(self)
    }
}


impl EdwardsPoint {
    /// Double the given point following:
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson, 
    /// http://eprint.iacr.org/2008/522, Section 3.1.
    /// Cost: 4M+ 4S+ 1D
    pub fn double(&self) -> EdwardsPoint {
        // This algorithm will be using point addition as base
        // until we find out what problem are we experiencing
        // with the doubling formulae on Extended Coords.
        //
        // TODO: Fix this as prio to reduce ops number.

        /*let two: FieldElement = FieldElement::from(&2u8);

        let A = &self.X * &self.X;
        let B = &self.Y * &self.Y;
        let C = &two * &(&self.Z * &self.Z);
        // a = -1
        let D = -&A;
        let E = &(&(&(&self.X + &self.Y) * &(&self.X + &self.Y)) - &A) - &B;
        let G = &D + &B;
        let F = &G - &C;
        let H = &D - &B;

        EdwardsPoint {
            X: &E * &F,
            Y: &G * &H,
            Z: &F * &G,
            T: &E * &H
        }*/

        self + self
    }

    /// Return the `EdwardsPoint` with Extended Coordinates
    /// eq to: {0, 0, 0, 0}.
    pub fn zero() -> EdwardsPoint {
        EdwardsPoint {
            X: FieldElement::zero(),
            Y: FieldElement::zero(),
            Z: FieldElement::zero(),
            T: FieldElement::zero(),
        }
    }

    /// Convert this `EdwardsPoint` on the Edwards model to the
    /// corresponding `MontgomeryPoint` on the Montgomery model.
    pub fn to_montgomery(&self) -> MontgomeryPoint {
       unimplemented!()
    }

    /// Compress this point to `CompressedEdwardsY` format.
    pub fn compress(&self) -> CompressedEdwardsY {
        let recip = self.Z.savas_koc_inverse();
        let x = &self.X * &recip;
        let y = &self.Y * &recip;
        let mut s: [u8; 32];

        s = y.to_bytes();
        // The 31st byte is always even, so we can't have a one by default.
        // That's why we can play with this bit to represent the sign info.
        s[31] ^= x.is_negative().unwrap_u8() << 7;
        CompressedEdwardsY(s)
    }

    /// Implementation of the standard algorithm of `double_and_add`.
    /// 
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
    /// Guide to Elliptic Curve Cryptography. 
    /// Springer Professional Computing. New York: Springer-Verlag. 
    /// 
    /// We implement this and not `windowing algorithms` because we 
    /// prioritize less constraints on R1CS over the computational 
    /// costs of the algorithm. Since, building circuits for ZK proofs
    /// it'll be more important.
    pub fn double_and_add(&self, s: &Scalar) -> EdwardsPoint {
        let mut N = self.clone();
        let mut n = s.clone();
        let mut Q = EdwardsPoint::identity();

        while n != Scalar::zero() {
            if !n.is_even() {
                Q = &Q + &N;
            };

            N = N.double();
            n = n.half();
        }  

        Q
    }
    

    /// Multiply by the cofactor: return (8 P).
    pub fn mul_by_cofactor(&self) -> EdwardsPoint {
        unimplemented!()
    }

    /// Compute ([2^k] P)
    pub fn mul_by_pow_2(&self, _k: u32) -> EdwardsPoint {
        unimplemented!()
    }
}

/// Module used for tesing `EdwardsPoint` operations and implementations.
/// Also used to check the correctness of the transformation functions.
pub mod tests {
    use super::*;

    pub static P1: EdwardsPoint = EdwardsPoint {
        X: FieldElement([23, 0, 0, 0, 0]),
        Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([4351986304670635, 4020128726404030, 674192131526433, 1158854437106827, 6468984742885])
    };

    pub static P2: EdwardsPoint = EdwardsPoint {
        X: FieldElement([68, 0, 0, 0, 0]),
        Y: FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([3505259403500377, 292342788271022, 2608000066641474, 796697979921534, 2995435405555])
    };

    /// `P4 = P1 + P2` over Twisted Edwards Extended Coordinates. 
    pub static P4: EdwardsPoint = EdwardsPoint {
        X: FieldElement([1933054591726350, 4403792816774408, 3029093546253310, 1134491999944368, 8146384875494]),
        Y: FieldElement([3369514960042642, 4355698098047571, 1547650124195635, 1314697306673062, 12051634278308]),
        Z: FieldElement([576719056868307, 1763329757922533, 3184642959683715, 2550235128581121, 11094626825862]),
        T: FieldElement([4345938036071968, 1280347559053553, 3286762790776823, 3577757860131876, 6505793015434])
    };

    /// `P3 = 2P1` over Twisted Edwards Extended Coordinates. 
    pub static P3: EdwardsPoint = EdwardsPoint {
        X: FieldElement([3851124475403222, 3539758816612178, 1146717153316815, 2152796892260637, 5956037993247]),
        Y: FieldElement([980361497621373, 1671502808757874, 2143986549518967, 1109176323830729, 9039277193734]),
        Z: FieldElement([2942534902618579, 3556685217095302, 1974617438797742, 1657071371119364, 16635295697052]),
        T: FieldElement([2487305805734419, 681684275336734, 499518740758148, 156812857292600, 3978688323434])
    };


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

    #[test]
    fn point_addition_extended_coords() {
        let res = &P1 + &P2;

        assert!(res == P4);
    }

    #[test]
    fn point_doubling_double_adding_extended_coords() {
        let res: EdwardsPoint = &P1 + &P1;
        
        assert!(res == P3);
    }
}
