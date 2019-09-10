//! Implementation of the Ristretto Protocol over the
//! Doppio curve.
//! 
//! Notes extracted from: https://ristretto.group/ristretto.html. 
//! Go there for the full lecture or check the paper here: 
//! https://tools.ietf.org/pdf/draft-hdevalence-cfrg-ristretto-00.pdf
//! 
//! # What's Ristretto?
//! Ristretto is a construction of a prime-order group using a non-prime-order Edwards curve.
//! The Decaf paper suggests using a non-prime-order curve E\mathcal EE to implement a prime-order
//! group by constructing a quotient group. Ristretto uses the same idea, but with different formulas,
//! in order to allow the use of cofactor-888 curves such as Curve25519.
//! 
//! Internally, a Ristretto point is represented by an Edwards point. 
//! Two Edwards points P,QP, QP,Q may represent the same Ristretto point, in the same way that 
//! different projective (X,Y,Z) coordinates may represent the same Edwards point. 
//! 
//! Group operations on Ristretto points are carried out with no overhead by performing the 
//! operations on the representative Edwards points.



use crate::constants;
use crate::edwards::{EdwardsPoint, double_and_add};
use crate::field::FieldElement;
use crate::scalar::Scalar;
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
        println!("Is positive? {:?} Is correct encoded?{:?}", s_is_positive.unwrap_u8(), s_correct_enc.unwrap_u8());
        if s_is_positive.unwrap_u8() == 0u8 || s_correct_enc.unwrap_u8() == 0u8 {
            return None
        };
        // Step 2: Attempt to decompress the CompressedRistretto. 
        let one = FieldElement::one();

        // u1 = 1 + as² with a = -1. 
        let u1 = one - s.square();
        // u2 = 1 - as² with a = -1. 
        let u2 = one + s.square();
        let u2_sq = u2.square();

        // v = a*d*u1² - u2²
        let v = -(constants::EDWARDS_D * u1.square()) - u2_sq;
        // I = 1/sqrt(v*u2²), returns `None` if the sqrt does not exist. 
        let (ok, I) = (v*u2_sq).inv_sqrt();
        if ok.unwrap_u8() == 0 {return None};
        
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

        println!("Is QR? {:?}  Is positive? {:?} Y is zero?  {:?}",ok.unwrap_u8(), t.is_positive().unwrap_u8(), y);
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

#[derive(Clone, Copy)]
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
        let X1Y2 = &self.0.X * &other.0.Y;
        let Y1X2 = &self.0.Y * &other.0.X;
        let X1X2 = &self.0.X * &other.0.X;
        let Y1Y2 = &self.0.Y * &other.0.Y;

        X1Y2.ct_eq(&Y1X2) | X1X2.ct_eq(&Y1Y2)
    }
}

impl PartialEq for RistrettoPoint {
    fn eq(&self, other: &RistrettoPoint) -> bool {
        self.compress().ct_eq(&other.compress()).unwrap_u8() == 1u8      
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

impl<'a> Neg for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    /// Negates a `RistrettoPoint` giving it's negated representation
    /// as a result.
    /// 
    /// Since the negative of a point is (-X:Y:Z:-T), it
    /// gives as a result: `(-X:Y:Z:-T)`.
    fn neg(self) -> RistrettoPoint {
        RistrettoPoint(-self.0)
    }
}

impl Neg for RistrettoPoint {
    type Output = RistrettoPoint;
    /// Negates a `RistrettoPoint` giving it's negated representation
    /// as a result.
    /// 
    /// Since the negative of a point is (-X:Y:Z:-T), it
    /// gives as a result: `(-X:Y:Z:-T)`.
    fn neg(self) -> RistrettoPoint {
        RistrettoPoint(-&self.0)
    }
}

impl<'a, 'b> Add<&'a RistrettoPoint> for &'b RistrettoPoint {
    type Output = RistrettoPoint;
    /// Performs the addition of two RistrettoPoints following the 
    /// Twisted Edwards Extended Coordinates formulae but using as
    /// `d` parameter the `RISTRETTO_D` instead of the `EDWARDS_D`.
    /// 
    /// This implementation is specific for curves with `a = -1` as 
    /// the isomorphic twist is for Doppio.
    /// 
    /// [Source: 2008 Hisil–Wong–Carter–Dawson], 
    /// (http://eprint.iacr.org/2008/522), Section 3.1.
    fn add(self, other: &'a RistrettoPoint) -> RistrettoPoint {
        unimplemented!()
    }    
}

impl Add<RistrettoPoint> for RistrettoPoint {
    type Output = RistrettoPoint;
    /// Performs the addition of two RistrettoPoints following the 
    /// Twisted Edwards Extended Coordinates formulae but using as
    /// `d` parameter the `RISTRETTO_D` instead of the `EDWARDS_D`.
    /// 
    /// This implementation is specific for curves with `a = -1` as 
    /// the isomorphic twist is for Doppio.
    /// 
    /// [Source: 2008 Hisil–Wong–Carter–Dawson], 
    /// (http://eprint.iacr.org/2008/522), Section 3.1.
    fn add(self, other: RistrettoPoint) -> RistrettoPoint {
        &self + &other
    }
}

impl<'a> Double for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    /// Performs the point doubling operation
    /// ie. `2*P` over the Twisted Edwards Extended
    /// Coordinates.
    /// 
    /// This implementation is specific for curves with `a = -1` as 
    /// the isomorphic twist is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson, 
    /// http://eprint.iacr.org/2008/522, Section 3.1.
    /// Cost: 4M+ 4S+ 1D
    fn double(self) -> RistrettoPoint {
        RistrettoPoint(self.0.double())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    /// Scalar multiplication: compute `self * Scalar`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    /// 
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
    /// Guide to Elliptic Curve Cryptography. 
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: &'b Scalar) -> RistrettoPoint {
        double_and_add(self, scalar)
    }
}

impl Mul<Scalar> for RistrettoPoint {
    type Output = RistrettoPoint;
    /// Scalar multiplication: compute `self * Scalar`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    /// 
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
    /// Guide to Elliptic Curve Cryptography. 
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: Scalar) -> RistrettoPoint {
        &self * &scalar
    }
}


impl RistrettoPoint {
    /// Encode a Ristretto point represented by the point `(X:Y:Z:T)`
    /// in extended coordinates. 
    #[allow(non_snake_case)]
    pub fn compress(&self) -> CompressedRistretto {
        let u1 = (self.0.Z + self.0.Y) * (self.0.Z - self.0.Y);
        let u2 = self.0.X * self.0.Y;
        println!("{:?}", (u1 * u2.square()).inv_sqrt());
        let (_, I) = (u1 * u2.square()).inv_sqrt();
        let D1 = u1 * I;
        let D2 = u2 * I;
        let Zinv = D1 * D2 * self.0.T;
        let mut xy;
        let D;
        if (self.0.T * Zinv).is_positive().unwrap_u8() == 0u8 {
            xy = 
                (constants::SQRT_MINUS_ONE * self.0.Y,
                constants::SQRT_MINUS_ONE * self.0.X);
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

    /// Computes the Ristretto Elligator map.
    ///
    /// # Note
    ///
    /// This method is not public because it's just used for hashing
    /// to a point -- proper elligator support is deferred for now.
    #[allow(non_snake_case)]
    pub(crate) fn elligator_ristretto_flavor(r_0: &FieldElement) -> RistrettoPoint {
        let (i, d) = (&constants::MINUS_SQRT_A, &constants::EDWARDS_D);
        let one = FieldElement::one();
        let one_minus_d_sq = &one - &d.square();
        let d_minus_one_sq = (d - &one).square();

        let r = i * &r_0.square();
        let N_s = &(&r + &one) * &one_minus_d_sq;
        let mut c = -&one;
        let D = &(&c - &(d * &r)) * &(&r + d);

        let (Ns_D_is_sq, mut s) = N_s.sqrt_ratio_i(&D);

        let mut s_prime = &s * r_0;
        let s_prime_is_pos = s_prime.is_positive();
        s_prime.conditional_negate(s_prime_is_pos);

        s.conditional_assign(&s_prime, !Ns_D_is_sq);
        c.conditional_assign(&r, !Ns_D_is_sq);

        let N_t = &(&(&c * &(&r - &one)) * &d_minus_one_sq) - &D;
        let s_sq = s.square();


        // The conversion from W_i is exactly the conversion from P1xP1.
        RistrettoPoint(EdwardsPoint{
            X: &(&s + &s) * &D,
            Z: &N_t * &constants::SQRT_AD_MINUS_ONE,
            Y: &FieldElement::one() - &s_sq,
            T: &FieldElement::one() + &s_sq,
        })
    }
}

mod tests {
    use hex;
    use super::*;

    #[test]
    fn point_comp_decomp_equivalence() {
        // This should give the RISTRETTO_BASEPOINT. 
        let bc = &constants::RISTRETTO_BASEPOINT.compress();
        println!("{:?}", bc);
        let bcd = bc.decompress().unwrap();
        println!("{:?}", bcd);
        let bcdc = bcd.compress();
        println!("{:?}", bcdc);
        let bcdcd = bcdc.decompress().unwrap();
        println!("{:?}", bcdcd);

        for point in &bcdcd.0.coset4() {
            println!("{:?}", point);
        };

        assert!(bcdcd == constants::RISTRETTO_BASEPOINT);
    }

    #[test]
    fn gqerhqerhreq() {
        let a = EdwardsPoint::new_random_point();
        println!("{:?}", a);
        for point in &constants::FOUR_COSET_GROUP {
            let p = &a + point;
            println!("{:?}", p);
            println!("Y pos? {:?}", p.Y.is_positive());
            println!("X*Y pos? {:?}", (p.Y * p.X).is_positive());
        };
    }


    #[test]
    fn rhqwehrwth() {
        let (ok, res) = (constants::EDWARDS_A - constants::EDWARDS_D).inv_sqrt();

        assert!(res == constants::INV_SQRT_A_MINUS_D);
    }

    #[test]
    fn ehrthrt() {
        let basepoint = EdwardsPoint {
            X: FieldElement([3461346045645288, 972387018214097, 4435206378704739, 3440261531857766, 17448879694677]),
            Y: FieldElement([4083500257631147, 3305309464452868, 1629709588575767, 4306635512831061, 7662705179761]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([212017923156762, 2250920437320259, 2648273971913259, 353963282208220, 3122143411880])
        };

        println!("{:?}", basepoint * Scalar::from(&8u8));
    }
}

