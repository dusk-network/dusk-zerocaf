#![allow(non_snake_case)]
//! Edwards Point operation implementations and definitions.
//! Encoding/decoding processes implementations
//! and support for all kind of interactions with them.

use crate::field::FieldElement;
use crate::scalar::Scalar;
use crate::montgomery::MontgomeryPoint;
use crate::constants;
use crate::traits::Identity;
use crate::traits::ops::*;


use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use rand::{Rng, thread_rng};

use std::default::Default;
use std::fmt::Debug;

use core::ops::{Index, IndexMut};
use std::ops::{Add, Sub, Mul, Neg};



// ------------- Common Point fn declarations ------------- //

/// Implementation of the standard algorithm of `double_and_add`.
/// This is a function implemented for Generic points that have
/// implemented `Add`, `Double`, `Identity` and `Clone`.
/// 
/// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
/// Guide to Elliptic Curve Cryptography. 
/// Springer Professional Computing. New York: Springer-Verlag. 
/// 
/// We implement this and not `windowing algorithms` because we 
/// prioritize less constraints on R1CS over the computational 
/// costs of the algorithm.
pub fn double_and_add<'b, 'a, T>(point: &'a T, scalar: &'b Scalar) -> T 
    where for<'c> &'c T: Add<Output = T>  + Double<Output = T>,
    T: Identity + Clone {

    let mut N = point.clone();
    let mut n = scalar.clone();
    let mut Q = T::identity();

    while n != Scalar::zero() {
        if !n.is_even() {
            Q = &Q + &N;
        };

        N = N.double();
        n = n.inner_half();
    }  
    Q
}

 /// Multiply by the cofactor: return (8 P).
pub fn mul_by_cofactor<'a, T>(point: &'a T) -> T 
    where for<'c> &'c T: Mul<&'c Scalar, Output = T> {
    point * &Scalar::from(&8u8)
}

/// Compute ([2^k] * P (mod l)).
/// 
/// Note: The maximum pow allowed is 249 since otherways
/// we will be able to get results greater than the
/// prime of the sub-group.
pub fn mul_by_pow_2<'a, 'b, T>(point: &'a T, _k: &'b u64) -> T 
    where for<'c> &'c T: Mul<&'c Scalar, Output = T> {
    point * &Scalar::two_pow_k(_k)
}



// ---------------- Point Structs ---------------- //

/// The first 255 bits of a `CompressedEdwardsY` represent the
/// (y)-coordinate.  The high bit of the 32nd byte gives the sign of (x).
#[derive(Copy, Clone)]
pub struct CompressedEdwardsY(pub [u8; 32]);

impl ConstantTimeEq for CompressedEdwardsY {
    fn ct_eq(&self, other: &CompressedEdwardsY) -> Choice {
        self.to_bytes().ct_eq(&other.to_bytes())
    }
}

impl PartialEq for CompressedEdwardsY {
    fn eq(&self, other: &CompressedEdwardsY) -> bool {
        self.ct_eq(&other).unwrap_u8() == 1u8
    }
}

impl Eq for CompressedEdwardsY {}

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
/// over the Twisted Edwards Extended Coordinates eg. (X, Y, Z, T).
/// 
/// Extended coordinates represent x y as`(X Y Z T)` satisfying the following equations:
/// x=X/Z
/// y=Y/Z
/// x*y=T/Z
#[derive(Copy, Clone)] 
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

impl ConstantTimeEq for EdwardsPoint {
    fn ct_eq(&self, other: &EdwardsPoint) -> Choice {
        AffinePoint::from(self).ct_eq(&AffinePoint::from(other))
    }
}

impl PartialEq for EdwardsPoint {
    fn eq(&self, other: &EdwardsPoint) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl Eq for EdwardsPoint {}

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

impl<'a> From<&'a ProjectivePoint> for EdwardsPoint {
    /// Given (X:Y:Z) in ε passing to εε can beperformed 
    /// in 3M+ 1S by computing (X*Z, Y*Z, X*Y, Z^2). 
    /// 
    /// Twisted Edwards Curves Revisited - 
    /// Huseyin Hisil, Kenneth Koon-Ho Wong, Gary Carter, 
    /// and Ed Dawson, Section 3.
    fn from(point: &'a ProjectivePoint) -> EdwardsPoint {
        EdwardsPoint {
            X: &point.X * &point.Z,
            Y: &point.Y * &point.Z,
            Z: point.Z.square(),
            T: &point.X * &point.Y
        }
    }
}

impl<'a> From<&'a AffinePoint> for EdwardsPoint {
    /// In affine form, each elliptic curve point has 2 coordinates, 
    /// like (x,y). In the new projective form, 
    /// each point will have 3 coordinates, like (X,Y,Z), 
    /// with the restriction that Z is never zero.
    ///  
    /// The forward mapping is given by (x,y)→(xz,yz,z), 
    /// for any non-zero z (usually chosen to be 1 for convenience).
    /// 
    /// After this is done, we move from Projective to Extended by
    /// setting the new coordinate `T = X * Y`.
    fn from(point: &'a AffinePoint) -> EdwardsPoint {
        EdwardsPoint {
            X: point.X,
            Y: point.Y,
            Z: FieldElement::one(),
            T: point.X * point.Y
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
    /// 
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
        let H = &B + &A;

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
    /// Scalar multiplication: compute `self * Scalar`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    /// 
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
    /// Guide to Elliptic Curve Cryptography. 
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: &'b Scalar) -> EdwardsPoint {
        double_and_add(self, scalar)
    }
}

impl Mul<Scalar> for EdwardsPoint {
    type Output = EdwardsPoint;
    /// Scalar multiplication: compute `Scalar * self`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    /// 
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
    /// Guide to Elliptic Curve Cryptography. 
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: Scalar) -> EdwardsPoint {
        double_and_add(&self, &scalar)
    }
}

impl<'a> Double for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Performs the point doubling operation
    /// ie. `2*P` over the Twisted Edwards Extended
    /// Coordinates.
    /// 
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson, 
    /// http://eprint.iacr.org/2008/522, Section 3.1.
    /// Cost: 4M+ 4S+ 1D
    fn double(self) -> EdwardsPoint {
        let two: FieldElement = FieldElement::from(&2u8);
        
        let A = self.X.square();
        let B = self.Y.square();
        let C = two * self.Z.square();
        let D = -A;
        let E = (self.X + self.Y) * (self.X + self.Y) -A -B;
        let G = D + B;
        let F = G - C;
        let H = D - B;

        EdwardsPoint {
            X: E * F,
            Y: G * H,
            Z: F * G,
            T: E * H
        }
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

    /// This function tries to build a Point over the Doppio Curve from
    /// a `Y` coordinate and a Choice that determines the Sign o the `X`
    /// coordinate that the user wants to use.
    /// 
    /// The function gets `X` by solving: 
    /// `+-X = mod_sqrt((y^2 -1)/(dy^2 - a))`. 
    /// 
    /// The sign of `x` is choosen with a `Choice` parameter. 
    /// 
    /// For Choice(0) -> Negative result. 
    /// For Choice(1) -> Positive result.
    /// 
    /// Then Z is always equal to `1`.
    /// 
    /// # Returns
    /// `Some(EdwardsPoint)` if there exists a result for the `mod_sqrt`
    /// and `None` if the resulting `x^2` isn't a QR modulo `FIELD_L`. 
    pub fn new_from_y_coord(y: &FieldElement, sign: Choice) -> Option<EdwardsPoint> {
        match ProjectivePoint::new_from_y_coord(&y, sign) {
            None => return None,
            Some(point) => return Some(EdwardsPoint::from(&point)),
        }
    } 

    /// This function tries to build a Point over the Doppio Curve from
    /// a random `Y` coordinate and a random Choice that determines the 
    /// Sign o the `X` coordinate.
    pub fn new_random_point() -> EdwardsPoint {
        // Simply generate a random `ProjectivePoint`
        // and once we get one that is valid, switch 
        // it to Extended Coordinates.
        EdwardsPoint::from(&ProjectivePoint::new_random_point())   
    }
}

/// A `ProjectivePoint` represents a point on the Doppio Curve expressed
/// over the Twisted Edwards Projective Coordinates eg. (X, Y, Z).
///  
/// For Z1≠0 the point (X1:Y1:Z1) represents the affine point (x1= X1/Z1, y1= Y1/Z1)
/// on EE,a,d.
/// Projective coordinates represent `x` `y` as `(X Y Z`) satisfying the following equations:
/// 
/// x=X/Z
/// 
/// y=Y/Z
/// 
/// Expressing an elliptic curve in twisted Edwards form saves time in arithmetic, 
/// even when the same curve can be expressed in the Edwards form. 
#[derive(Copy, Clone)]
pub struct ProjectivePoint {
    pub X: FieldElement,
    pub Y: FieldElement,
    pub Z: FieldElement
}

impl Debug for ProjectivePoint {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "
        ProjectivePoint {{
            X: {:?},
            Y: {:?},
            Z: {:?}
        }};", self.X, self.Y, self.Z)
    }
}

impl ConstantTimeEq for ProjectivePoint {
    fn ct_eq(&self, other: &ProjectivePoint) -> Choice {
        AffinePoint::from(self).ct_eq(&AffinePoint::from(other))
    }
}

impl PartialEq for ProjectivePoint {
    fn eq(&self, other: &ProjectivePoint) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl Eq for ProjectivePoint {}

impl Default for ProjectivePoint {
    /// Returns the default ProjectivePoint Extended Coordinates: (0, 1, 1). 
    fn default() -> ProjectivePoint {
        ProjectivePoint::identity()
    }
}

impl Identity for ProjectivePoint {
    /// Returns the Edwards Point identity value = `(0, 1, 1)`.
    fn identity() -> ProjectivePoint {
        ProjectivePoint {
            X: FieldElement::zero(),
            Y: FieldElement::one(),
            Z: FieldElement::one()
        }
    }
}

impl<'a> From<&'a EdwardsPoint> for ProjectivePoint {
    /// Given (X:Y:T:Z) in εε, passing to ε is cost-free by 
    /// simply ignoring `T`.
    /// 
    /// Twisted Edwards Curves Revisited - 
    /// Huseyin Hisil, Kenneth Koon-Ho Wong, Gary Carter, 
    /// and Ed Dawson, Section 3.
    fn from(point: &'a EdwardsPoint) -> ProjectivePoint {
        ProjectivePoint{
            X: point.X,
            Y: point.Y,
            Z: point.Z
        }
    }
}

impl<'a> From<&'a AffinePoint> for ProjectivePoint {
    /// The key idea of projective coordinates is that instead of 
    /// performing every division immediately, we defer the 
    /// divisions by multiplying them into a denominator.
    /// 
    /// In affine form, each elliptic curve point has 2 coordinates, 
    /// like (x,y). In the new projective form, 
    /// each point will have 3 coordinates, like (X,Y,Z), 
    /// with the restriction that Z is never zero.
    ///  
    /// The forward mapping is given by (x,y)→(xz,yz,z), 
    /// for any non-zero z (usually chosen to be 1 for convenience).
    fn from(point: &'a AffinePoint) -> ProjectivePoint {
        ProjectivePoint {
            X: point.X,
            Y: point.Y,
            Z: FieldElement::one()
        }
    }
}

impl<'a> Neg for &'a ProjectivePoint {
    type Output = ProjectivePoint;
    /// Negates an `ProjectivePoint` giving it as a result.
    /// Since the negative of a point is (-X:Y:Z:-T), it
    /// gives as a result: `(-X, Y, Z, -T)`.
    fn neg(self) -> ProjectivePoint {
       ProjectivePoint{
           X: -&self.X,
           Y:   self.Y,
           Z:   self.Z
       }
    }
}

impl Neg for ProjectivePoint {
    type Output = ProjectivePoint;
    /// Negates an `ProjectivePoint` giving it as a result
    fn neg(self) -> ProjectivePoint {
        -&self
    }
}

impl<'a, 'b> Add<&'b ProjectivePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;
    /// Add two ProjectivePoints and give the resulting `ProjectivePoint`.
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// 
    /// Bernstein D.J., Birkner P., Joye M., Lange T., Peters C. 
    /// (2008) Twisted Edwards Curves. In: Vaudenay S. (eds) 
    /// Progress in Cryptology – AFRICACRYPT 2008. AFRICACRYPT 2008. 
    /// Lecture Notes in Computer Science, vol 5023. Springer, Berlin, Heidelberg.
    /// See: https://eprint.iacr.org/2008/013.pdf - Section 6.
    #[inline]
    fn add(self, other: &'b ProjectivePoint) -> ProjectivePoint {
        let A = self.Z * other.Z;
        let B = A.square();
        let C = self.X * other.X;
        let D = self.Y * other.Y;
        let E = constants::EDWARDS_D * C * D;
        let F = B - E;
        let G = B + E;

        ProjectivePoint {
            X: A * (F * (((self.X + self.Y) * (other.X + other.Y) - C) - D)),
            Y: A * G * (D + C),
            Z: F * G
        }
    }
}

impl Add<ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;
    /// Add two ProjectivePoints and give the resulting `ProjectivePoint`.
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// 
    /// Bernstein D.J., Birkner P., Joye M., Lange T., Peters C. 
    /// (2008) Twisted Edwards Curves. In: Vaudenay S. (eds) 
    /// Progress in Cryptology – AFRICACRYPT 2008. AFRICACRYPT 2008. 
    /// Lecture Notes in Computer Science, vol 5023. Springer, Berlin, Heidelberg.
    /// See: https://eprint.iacr.org/2008/013.pdf - Section 6.
    #[inline]
    fn add(self, other: ProjectivePoint) -> ProjectivePoint {
        &self + &other
    }
}

impl<'a, 'b> Sub<&'b ProjectivePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;
    /// Add two ProjectivePoints, negating the second one,
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// 
    /// Bernstein D.J., Birkner P., Joye M., Lange T., Peters C. 
    /// (2008) Twisted Edwards Curves. In: Vaudenay S. (eds) 
    /// Progress in Cryptology – AFRICACRYPT 2008. AFRICACRYPT 2008. 
    /// Lecture Notes in Computer Science, vol 5023. Springer, Berlin, Heidelberg.
    /// See: https://eprint.iacr.org/2008/013.pdf - Section 6.
    fn sub(self, other: &'b ProjectivePoint) -> ProjectivePoint {
        self + &(-other)
    }
}

impl Sub<ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;
    /// Add two ProjectivePoints, negating the second one,
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// 
    /// Bernstein D.J., Birkner P., Joye M., Lange T., Peters C. 
    /// (2008) Twisted Edwards Curves. In: Vaudenay S. (eds) 
    /// Progress in Cryptology – AFRICACRYPT 2008. AFRICACRYPT 2008. 
    /// Lecture Notes in Computer Science, vol 5023. Springer, Berlin, Heidelberg.
    /// See: https://eprint.iacr.org/2008/013.pdf - Section 6.
    #[inline]
    fn sub(self, other: ProjectivePoint) -> ProjectivePoint {
        &self - &other
    }
}

impl<'a, 'b> Mul<&'a Scalar> for &'b ProjectivePoint {
    type Output = ProjectivePoint;

    /// Scalar multiplication: compute `Scalar * self`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    /// 
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
    /// Guide to Elliptic Curve Cryptography. 
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: &'a Scalar) -> ProjectivePoint {
        double_and_add(self, scalar)
    }
}

impl Mul<Scalar> for ProjectivePoint {
    type Output = ProjectivePoint;

    /// Scalar multiplication: compute `Scalar * self`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    /// 
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004). 
    /// Guide to Elliptic Curve Cryptography. 
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: Scalar) -> ProjectivePoint {
        &self * &scalar
    }
}

impl<'a> Double for &'a ProjectivePoint {
    type Output = ProjectivePoint;
    /// Double the given point following:
    /// This implementation is specific for curves with `a = -1` as Doppio is.
    /// 
    /// /// Bernstein D.J., Birkner P., Joye M., Lange T., Peters C. 
    /// (2008) Twisted Edwards Curves. In: Vaudenay S. (eds) 
    /// Progress in Cryptology – AFRICACRYPT 2008. AFRICACRYPT 2008. 
    /// Lecture Notes in Computer Science, vol 5023. Springer, Berlin, Heidelberg.
    /// See: https://eprint.iacr.org/2008/013.pdf - Section 6.
    /// 
    /// Cost: 3M+ 4S+ +7a + 1D.
    fn double(self) -> ProjectivePoint {
        let B = (&self.X + &self.Y).square();
        let C = self.X.square();
        let D = self.Y.square();
        let E = &constants::EDWARDS_A * &C;
        let F = &E + &D;
        let H = self.Z.square();
        let J = &F - &(&FieldElement::from(&2u8) * &H);

        ProjectivePoint {
            X: &(&(&B - &C) -&D) * &J,
            Y: &F * &(&E - &D),
            Z: &F *&J
        }
    }
}

impl ProjectivePoint {

    /// This function tries to build a Point over the Doppio Curve from
    /// a `Y` coordinate and a Choice that determines the Sign o the `X`
    /// coordinate that the user wants to use.
    /// 
    /// The function gets `X` by solving: 
    /// `+-X = mod_sqrt((y^2 -1)/(dy^2 - a))`. 
    /// 
    /// The sign of `x` is choosen with a `Choice` parameter. 
    /// 
    /// For Choice(0) -> Negative result. 
    /// For Choice(1) -> Positive result.
    /// 
    /// Then Z is always equal to `1`.
    /// 
    /// # Returns
    /// `Some(ProjectivePoint)` if there exists a result for the `mod_sqrt`
    /// and `None` if the resulting `x^2` isn't a QR modulo `FIELD_L`. 
    pub fn new_from_y_coord(y: &FieldElement, sign: Choice) -> Option<ProjectivePoint> {
        let x;
        let xx = (y.square() - FieldElement::one()) / ((constants::EDWARDS_D * y.square()) - constants::EDWARDS_A);
        // Match the sqrt(x) Option.
        match xx.mod_sqrt(sign) {
            // If we get `None`, we know that 
            None => return None,
            Some(_x) => x = _x,
        };

        Some(ProjectivePoint {
            X: x,
            Y: y.clone(),
            Z: FieldElement::one()
        })
    }

    /// This function tries to build a Point over the Doppio Curve from
    /// a random `Y` coordinate and a random Choice that determines the 
    /// Sign o the `X` coordinate.
    pub fn new_random_point() -> ProjectivePoint {
        // Gen a random `Y` coordinate value.
        let y = FieldElement::generate_random();
        // Gen a random sign choice. 
        let sign = Choice::from(thread_rng().gen_range(0u8, 1u8));

        // Until we don't get a valid `Y` value, call the
        // function recursively.
        match ProjectivePoint::new_from_y_coord(&y, sign) {
            None => ProjectivePoint::new_random_point(),
            Some(point) => return point,
        }    
    }
}

/// A `AffinePoint` represents a point on the Doppio Curve expressed
/// over the Twisted Edwards Affine Coordinates also known as 
/// cartesian coordinates: (X, Y). 
/// 
/// This Twisted Edwards coordinates are ONLY implemented for
/// equalty testing, since all of the Point operations
/// defined over them are much slower than the same ones 
/// defined over Extended or Projective coordinates. 
pub struct AffinePoint {
    pub X: FieldElement,
    pub Y: FieldElement
}

impl Debug for AffinePoint {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "
        AffinePoint {{
            X: {:?},
            Y: {:?}
        }};", self.X, self.Y)
    }
}

impl Default for AffinePoint {
    /// Returns the default Twisted Edwards AffinePoint Coordinates: (0, 1). 
    fn default() -> AffinePoint {
        AffinePoint::identity()
    }
}

impl Identity for AffinePoint {
    /// Returns the Edwards Point identity value = `(0, 1)`.
    fn identity() -> AffinePoint {
        AffinePoint {
            X: FieldElement::zero(),
            Y: FieldElement::one()
        }
    }
}

impl ConstantTimeEq for AffinePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.X.ct_eq(&other.X) & self.Y.ct_eq(&other.Y)
    }
}

impl PartialEq for AffinePoint {
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(&other))
    }
}

impl Eq for AffinePoint {}

impl<'a> From<&'a EdwardsPoint> for AffinePoint {
    /// Given (X:Y:Z:T) in εε, passing to affine can be 
    /// performed in 3M+ 1I by computing:
    /// 
    /// First, move to Projective Coordinates by removing `T`.
    ///  
    /// Then, reduce the point from Projective to Affine
    /// coordinates computing: (X*Zinv, Y*Zinv, Z*Zinv).
    /// 
    /// And once Z coord = `1` we can simply remove it. 
    /// 
    /// Twisted Edwards Curves Revisited - 
    /// Huseyin Hisil, Kenneth Koon-Ho Wong, Gary Carter, 
    /// and Ed Dawson.
    fn from(point: &'a EdwardsPoint) -> AffinePoint {
        let A = point.Z.inverse();
        AffinePoint {
            X: point.X * A,
            Y: point.Y * A
        }
    }
}

impl<'a> From<&'a ProjectivePoint> for AffinePoint {
    /// Reduce the point from Projective to Affine
    /// coordinates computing: (X*Zinv, Y*Zinv, Z*Zinv).
    /// 
    /// And once Z coord = `1` we can simply remove it. 
    /// 
    /// Twisted Edwards Curves Revisited - 
    /// Huseyin Hisil, Kenneth Koon-Ho Wong, Gary Carter, 
    /// and Ed Dawson.
    fn from(point: &'a ProjectivePoint) -> AffinePoint {
        AffinePoint {
            X: point.X * point.Z.inverse(),
            Y: point.Y * point.Z.inverse()
        }
    }
}

impl<'a> Neg for &'a AffinePoint {
    type Output = AffinePoint;
    /// Negates an `AffinePoint` giving it as a result.
    /// Since the negative of a point is (-X:Y), it
    /// gives as a result: `(-X, Y)`.
    fn neg(self) -> AffinePoint {
       AffinePoint{
           X: -&self.X,
           Y:   self.Y
       }
    }
}

impl Neg for AffinePoint {
    type Output = AffinePoint;
    /// Negates an `AffinePoint` giving it as a result.
    /// Since the negative of a point is (-X:Y), it
    /// gives as a result: `(-X, Y)`.
    fn neg(self) -> AffinePoint {
        -&self
    }
}



#[allow(dead_code)]
#[cfg(test)]
/// Module used for tesing `EdwardsPoint` operations and implementations.
/// Also used to check the correctness of the transformation functions.
pub mod tests {
    use super::*;

    pub(self) static P1_AFFINE: AffinePoint = AffinePoint {
        X: FieldElement([23, 0, 0, 0, 0]),
        Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998])
    };

    pub(self) static P2_AFFINE: AffinePoint = AffinePoint {
        X: FieldElement([68, 0, 0, 0, 0]),
        Y: FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]),
    };

    pub(self) static P1_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([23, 0, 0, 0, 0]),
        Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([4351986304670635, 4020128726404030, 674192131526433, 1158854437106827, 6468984742885])
    };

    pub(self) static P2_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([68, 0, 0, 0, 0]),
        Y: FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([3505259403500377, 292342788271022, 2608000066641474, 796697979921534, 2995435405555])
    };

    /// `P4_EXTENDED = P1_EXTENDED + P2_EXTENDED` over Twisted Edwards Extended Coordinates. 
    pub(self) static P4_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([1933054591726350, 4403792816774408, 3029093546253310, 1134491999944368, 8146384875494]),
        Y: FieldElement([3369514960042642, 4355698098047571, 1547650124195635, 1314697306673062, 12051634278308]),
        Z: FieldElement([576719056868307, 1763329757922533, 3184642959683715, 2550235128581121, 11094626825862]),
        T: FieldElement([4345938036071968, 1280347559053553, 3286762790776823, 3577757860131876, 6505793015434])
    };

    /// `P3_EXTENDED = 2P1_EXTENDED` over Twisted Edwards Extended Coordinates. 
    pub(self) static P3_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([3851124475403222, 3539758816612178, 1146717153316815, 2152796892260637, 5956037993247]),
        Y: FieldElement([980361497621373, 1671502808757874, 2143986549518967, 1109176323830729, 9039277193734]),
        Z: FieldElement([2942534902618579, 3556685217095302, 1974617438797742, 1657071371119364, 16635295697052]),
        T: FieldElement([2487305805734419, 681684275336734, 499518740758148, 156812857292600, 3978688323434])
    };

    pub(self) static P1_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([23, 0, 0, 0, 0]),
        Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
        Z: FieldElement([1, 0, 0, 0, 0])
    };

    pub(self) static P2_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([68, 0, 0, 0, 0]),
        Y: FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]),
        Z: FieldElement([1, 0, 0, 0, 0])
    };

    /// `P4_PROJECTIVE = P1_PROJECTIVE + P2_PROJECTIVE` over Twisted Edwards Projective Coordinates. 
    pub(self) static P4_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([1933054591726350, 4403792816774408, 3029093546253310, 1134491999944368, 8146384875494]),
        Y: FieldElement([3369514960042642, 4355698098047571, 1547650124195635, 1314697306673062, 12051634278308]),
        Z: FieldElement([576719056868307, 1763329757922533, 3184642959683715, 2550235128581121, 11094626825862])
    };

    /// `P3_PROJECTIVE = 2P1_PROJECTIVE` over Twisted Edwards Projective Coordinates. 
    pub(self) static P3_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([3851124475403222, 3539758816612178, 1146717153316815, 2152796892260637, 5956037993247]),
        Y: FieldElement([980361497621373, 1671502808757874, 2143986549518967, 1109176323830729, 9039277193734]),
        Z: FieldElement([2942534902618579, 3556685217095302, 1974617438797742, 1657071371119364, 16635295697052])
    };


    //------------------ Tests ------------------//

    #[test]
    fn from_projective_to_extended() {
        let p3_extended_proj = ProjectivePoint {
            X: FieldElement([23, 0, 0, 0, 0]),
            Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
            Z: FieldElement([1, 0, 0, 0, 0])
        };

        assert!(EdwardsPoint::from(&p3_extended_proj) == P1_EXTENDED);
    }

    #[test]
    fn from_extended_to_projective() {
        let p3_extended_proj = ProjectivePoint {
            X: FieldElement([23, 0, 0, 0, 0]),
            Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
            Z: FieldElement([1, 0, 0, 0, 0])
        };

        assert!(p3_extended_proj == ProjectivePoint::from(&P1_EXTENDED));
    }


    #[test]
    fn extended_point_neg() {
        let a = EdwardsPoint::default();
        let inv_a = EdwardsPoint {
            X: FieldElement::zero(),
            Y: FieldElement::one(),
            Z: FieldElement::one(),
            T: FieldElement::zero()
        };

        let res = -a;
        assert!(res == inv_a);
    }

    #[test]
    fn extended_coords_neg_identity() {
        let res = - &EdwardsPoint::identity();

        assert!(res == EdwardsPoint::identity())
    }

    #[test]
    fn extended_point_addition() {
        let res = &P1_EXTENDED + &P2_EXTENDED;

        assert!(res == P4_EXTENDED);
    }

    #[test]
    fn extended_point_doubling_by_addition() {
        let res: EdwardsPoint = &P1_EXTENDED + &P1_EXTENDED;
        
        assert!(res == P3_EXTENDED);
    }

    #[test]
    fn extended_point_doubling() {
        let res = P1_EXTENDED.double();
        
        assert!(res == P3_EXTENDED);
    }

    #[test]
    fn extended_double_and_add() {
        // Since we compute ( 8* P ) we don't need
        // to test `mul_by_cofactor` if this passes.
        let expect = P1_EXTENDED.double().double().double();
        let res = P1_EXTENDED * Scalar::from(&8u8);

        assert!(expect == res);
    }

    #[test]
    fn extended_point_generation() {
        let y2 = FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]);
        let p2 = EdwardsPoint::new_from_y_coord(&y2, Choice::from(0u8)).unwrap();

        assert!(p2 == P2_EXTENDED); 

        let y1 = FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]);
        let p1 = EdwardsPoint::new_from_y_coord(&y1, Choice::from(1u8)).unwrap();

        assert!(p1 == P1_EXTENDED);

        // `A = 15` which is not a QR over the prime of the field.
        let y_failure = FieldElement::from(&15u8);
        let p_fail = EdwardsPoint::new_from_y_coord(&y_failure, Choice::from(0u8));
        assert!(p_fail.is_none())
    }

    #[test]
    fn projective_point_neg() {
        let a = ProjectivePoint::default();
        let inv_a = ProjectivePoint {
            X: FieldElement::zero(),
            Y: FieldElement::one(),
            Z: FieldElement::one()
        };

        let res = -a;
        assert!(res == inv_a);
    }

    #[test]
    fn projective_coords_neg_identity() {
        let res = - &ProjectivePoint::identity();

        assert!(res == ProjectivePoint::identity())
    }

    #[test]
    fn projective_point_addition() {
        let res = P1_PROJECTIVE + P2_PROJECTIVE;
        
        assert!(res == P4_PROJECTIVE);
    }

    #[test]
    fn projective_point_doubling() {
        let res = P1_PROJECTIVE.double();

        assert!(res == P3_PROJECTIVE);
    }

    #[test]
    fn projective_double_and_add() {
        // Since we compute ( 8* P ) we don't need
        // to test `mul_by_cofactor` if this passes.
        let expect = P1_PROJECTIVE.double().double().double();
        let res = P1_PROJECTIVE * Scalar::from(&8u8);

        assert!(expect == res);
    }

    #[test]
    fn projective_point_generation() {
        let y2 = FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]);
        let p2 = ProjectivePoint::new_from_y_coord(&y2, Choice::from(0u8)).unwrap();

        assert!(p2 == P2_PROJECTIVE); 

        let y1 = FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]);
        let p1 = ProjectivePoint::new_from_y_coord(&y1, Choice::from(1u8)).unwrap();

        assert!(p1 == P1_PROJECTIVE);

        // `A = 15` which is not a QR over the prime of the field.
        let y_failure = FieldElement::from(&15u8);
        let p_fail = ProjectivePoint::new_from_y_coord(&y_failure, Choice::from(0u8));
        assert!(p_fail.is_none())
    } 

    #[test]
    fn affine_point_ct_eq() {
        // Test true case for equalty.
        assert!(P1_AFFINE.ct_eq(&P2_AFFINE).unwrap_u8() == 0u8);

        // Test false case for equalty.
        assert!(P1_AFFINE.ct_eq(&P1_AFFINE).unwrap_u8() == 1u8);

        // After Point doubling case. Z's are different but points are equivalent.
        assert!(AffinePoint::from(&P3_PROJECTIVE)
            .ct_eq(&AffinePoint::from(&P1_PROJECTIVE.double())).
            unwrap_u8() 
        == 1u8);
    }  

    #[test]
    fn affine_point_eq() {
        assert!(P1_AFFINE == P1_AFFINE);
        assert!(P1_AFFINE != P2_AFFINE);

        assert!(AffinePoint::from(&P3_PROJECTIVE) == AffinePoint::from(&P1_PROJECTIVE.double()));
    }
}
