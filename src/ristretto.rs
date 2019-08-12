//! Contains the implementation of the Ristretto protocol.


use crate::constants;
use crate::scalar::Ristretto255Scalar;
use crate::edwards::EdwardsPoint;
use crate::field::FieldElement;
use crate::traits::ops::*;
use crate::traits::Identity;

use core::ops::{Add, Sub, Mul, Neg, Index};

use std::fmt::Debug;

use subtle::{Choice, ConditionallyNegatable, ConditionallySelectable, ConstantTimeEq};


/// Ristretto Point expressed in wire format.
/// Since the Ristretto bytes encoding is canonical,
/// two points are equal if their encodin form is equal. 
#[derive(Debug, Clone, Copy)]
pub struct CompressedRistretto(pub [u8; 32]);

impl Index<usize> for CompressedRistretto {
    type Output = u8;
    fn index(&self, _index: usize) -> &u8 {
        &(self.0[_index])
    }
}

impl ConstantTimeEq for CompressedRistretto {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.as_bytes().ct_eq(&other.as_bytes())
    }
}

impl PartialEq for CompressedRistretto {
    fn eq(&self, other: &CompressedRistretto) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl Eq for CompressedRistretto {}

impl CompressedRistretto {
    /// Get the bytes of the `CompressedRistretto` point.
    pub fn as_bytes(&self) -> [u8; 32] {
        self.0
    }

    #[allow(non_snake_case)]
    /// Attempt to decompress a `CompressedRistretto` point. 
    /// This proces is done following the guidelines and steps
    /// mentioned in: https://ristretto.group/formulas/decoding.html
    /// 
    /// # Returns
    /// - If the decompression/decoding succeeds -> `Some(RistrettoPoint)`. 
    /// - If the decompression/decoding fails -> `None`.
    pub fn decompress(&self) -> Option<RistrettoPoint> {
        // Step 1: Check that the byte-string is a valid FieldElement. 

        // As Ristretto paper says: "If the implementation's field element
        // encoding function produces canonical outputs, one way to check 
        // that s_bytes is a canonical encoding (in step 1) is to decode 
        // s_bytes into sss, then re-encode sss into s_bytes_check, and ensure 
        // that s_bytes == s_bytes_check.
        let s: FieldElement = FieldElement::from_bytes(&self.as_bytes());
        let s_check = s.to_bytes();
        let s_correct_enc = s_check.ct_eq(&self.as_bytes());
        let s_is_positive = s.is_positive();

        // If the byte-encoding was incorrect or the representation is
        // a negative `FieldElement` (according to the definition of 
        // positive found on Decaf paper), return `None`. 
        if s_is_positive.unwrap_u8() == 0u8 || s_correct_enc.unwrap_u8() == 0u8 {
            return None
        };

        // Step 2: Attempt to decompress the CompressedRistretto. 
        let one = FieldElement::one();

        // u1 = 1 + as² with a = -1. 
        let u1 = one - s.square();
        // u2 = 1 - as² with a = -1. +
        let u2 = one + s.square();
        let u2_sq = u2.square();

        // v = a*d*u1² - u2²
        let v = -(constants::EDWARDS_D * u1.square()) - u2_sq;
        // I = 1/sqrt(v*u2²), returns `None` if the sqrt does not exist. 
        let I = match (v*u2_sq).inv_sqrt() {
            None => return None,
            Some(x) => x
        };
        // Compute the Extended Point Coordinates Y & T
        let Dx = I*u2;
        let Dy = I*Dx*v;

        // Compute ABS(2*s*Dx) and negate if it is negative.
        let mut x = (&s + &s) * Dx;
        let x_is_pos = x.is_positive();
        x.conditional_negate(!x_is_pos);
        // Compute Y and T coordinates. 
        let y = u1 * Dy;
        let t = x * y;

        if t.is_positive().unwrap_u8() == 0u8 || y == FieldElement::zero() {
            return None
        };

        Some(RistrettoPoint(EdwardsPoint {
            X: x,
            Y: y,
            Z: one,
            T: t
        }))
    }
}


pub struct RistrettoPoint (pub(crate) EdwardsPoint);

impl Debug for RistrettoPoint {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "{:?}", &self.0)
    }
}

impl ConstantTimeEq for RistrettoPoint {
    /// As specified on the Ristretto protocol docs: 
    /// https://ristretto.group/formulas/equality.html
    /// and we are on the twisted case, we compare 
    /// X1*Y2 == Y1*X2. 
    fn ct_eq(&self, other: &RistrettoPoint) -> Choice {
        (self.0.X * other.0.Y).to_bytes().ct_eq(
            &(self.0.Y * other.0.X).to_bytes())
    }
}

impl PartialEq for RistrettoPoint {
    fn eq(&self, other: &RistrettoPoint) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8      
    }
}

impl Eq for RistrettoPoint {}

impl Identity for RistrettoPoint {
    /// Gives back the Identity point for the Extended Edwards Coordinates
    /// which is endoded as a `RistrettoPoint` with coordinates: 
    /// `(X, Y, Z, T)` = `(0, 1, 1, 0)`.
    fn identity() -> RistrettoPoint {
        RistrettoPoint(EdwardsPoint::identity())
    }
}

impl Default for RistrettoPoint {
    /// Gives back the Identity point for the Extended Edwards Coordinates
    /// which is endoded as a `RistrettoPoint` with coordinates: 
    /// `(X, Y, Z, T)` = `(0, 1, 1, 0)`.
    fn default() -> RistrettoPoint {
        RistrettoPoint::identity()
    }
}



impl RistrettoPoint {
    /// Encode a Ristretto point represented by the point `(X:Y:Z:T)`
    /// in extended coordinates. 
    #[allow(non_snake_case)]
    pub fn compress(&self) -> CompressedRistretto {
        let u1 = (self.0.Z + self.0.Y) * (self.0.Z - self.0.Y);
        let u2 = self.0.X * self.0.Y;
        let I = (u1 * u2).inv_sqrt().expect("This is not a valid point representative");
        let D1 = u1 * I;
        let D2 = u2 * I;
        let Zinv = D1 * D2 * self.0.T;
        let mut xy;
        let D;
        if (self.0.T * Zinv).is_positive().unwrap_u8() == 0u8 {
            xy = 
                ((self.0.Y * constants::INV_SQRT_A),
                self.0.X * constants::MINUS_SQRT_A);
            D = D1 * constants::INV_SQRT_A_MINUS_D;
        } else {
            xy = (self.0.X, self.0.Y);
            D = D2;
        };

        xy.1.conditional_negate(!(xy.0*Zinv).is_positive());
        // We are on the Twisted case, so a = -1. 
        // Then s = ABS((Z-Y) * D)
        let mut s = (self.0.Z - xy.1) * D;
        s.conditional_negate(!s.is_positive());

        CompressedRistretto(s.to_bytes()) 
    }

    
}

mod tests {
    use hex;
    use super::*;


    #[test]
    fn bad_encodings() {
        // The following are invalid encodings, which should all be rejected.
        // These are designed to test each of the checks during decoding.

        let bad_encodings = [
            // These are all bad because they're non-canonical field encodings.
            "00ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "f3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "edffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            // These are all bad because they're negative field elements.
            "0100000000000000000000000000000000000000000000000000000000000000",
            "01ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "ed57ffd8c914fb201471d1c3d245ce3c746fcbe63a3679d51b6a516ebebe0e20",
            "c34c4e1826e5d403b78e246e88aa051c36ccf0aafebffe137d148a2bf9104562",
            "c940e5a4404157cfb1628b108db051a8d439e1a421394ec4ebccb9ec92a8ac78",
            "47cfc5497c53dc8e61c91d17fd626ffb1c49e2bca94eed052281b510b1117a24",
            "f1c6165d33367351b0da8f6e4511010c68174a03b6581212c71c0e1d026c3c72",
            "87260f7a2f12495118360f02c26a470f450dadf34a413d21042b43b9d93e1309",
            // These are all bad because they give a nonsquare x^2.
            "26948d35ca62e643e26a83177332e6b6afeb9d08e4268b650f1f5bbd8d81d371",
            "4eac077a713c57b4f4397629a4145982c661f48044dd3f96427d40b147d9742f",
            "de6a7b00deadc788eb6b6c8d20c0ae96c2f2019078fa604fee5b87d6e989ad7b",
            "bcab477be20861e01e4a0e295284146a510150d9817763caf1a6f4b422d67042",
            "2a292df7e32cababbd9de088d1d1abec9fc0440f637ed2fba145094dc14bea08",
            "f4a9e534fc0d216c44b218fa0c42d99635a0127ee2e53c712f70609649fdff22",
            "8268436f8c4126196cf64b3c7ddbda90746a378625f9813dd9b8457077256731",
            "2810e5cbc2cc4d4eece54f61c6f69758e289aa7ab440b3cbeaa21995c2f4232b",
            // These are all bad because they give a negative xy value.
            "3eb858e78f5a7254d8c9731174a94f76755fd3941c0ac93735c07ba14579630e",
            "a45fdc55c76448c049a1ab33f17023edfb2be3581e9c7aade8a6125215e04220",
            "d483fe813c6ba647ebbfd3ec41adca1c6130c2beeee9d9bf065c8d151c5f396e",
            "8a2e1d30050198c65a54483123960ccc38aef6848e1ec8f5f780e8523769ba32",
            "32888462f8b486c68ad7dd9610be5192bbeaf3b443951ac1a8118419d9fa097b",
            "227142501b9d4355ccba290404bde41575b037693cef1f438c47f8fbf35d1165",
            "5c37cc491da847cfeb9281d407efc41e15144c876e0170b499a96a22ed31e01e",
            "445425117cb8c90edcbc7c1cc0e74f747f2c1efa5630a967c64f287792a48a4b",
            // This is s = -1, which causes y = 0.
            "ecffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
        ];

        let mut bad_bytes = [0u8; 32];
        for bad_encoding in &bad_encodings {
            bad_bytes.copy_from_slice(&hex::decode(bad_encoding).unwrap());
            println!("Try: {:?}", CompressedRistretto(bad_bytes));
            assert!(CompressedRistretto(bad_bytes).decompress().is_none());
        }
    }
}

