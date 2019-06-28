//! This is the definition and implementation
//! of a Weierstrass point.

use core::ops::Add;
use core::ops::Mul;

use crate::field::FieldElement;
use crate::scalar::Scalar;

/// Holds the \\(u\\)-coordinate of a point on the Montgomery form of
/// Curve25519 or its twist.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct WeierstrassPoint(pub [u8; 32]);

impl Default for WeierstrassPoint {
    fn default() -> WeierstrassPoint {
        WeierstrassPoint([0u8; 32])
    }
}

impl Eq for WeierstrassPoint {}

impl WeierstrassPoint {
    /// View this `WeierstrassPoint` as an array of bytes.
    pub fn as_bytes<'a>(&'a self) -> &'a [u8; 32] {
        &self.0
    }

    /// Convert this `WeierstrassPoint` to an array of bytes.
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }
}

/// A `ProjectivePoint` represents a point on the Weierstrass form of the elliptic curve.
#[derive(Debug, PartialEq)]
pub struct ProjectivePoint {
    pub y: FieldElement,
}

impl Default for ProjectivePoint {
    fn default() -> ProjectivePoint {
        let field_elem = ProjectivePoint {
            y: FieldElement::one(),
        };
        field_elem
    }
}

impl<'a, 'b> Add<&'b ProjectivePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;
    fn add(self, _rhs: &'b ProjectivePoint) -> ProjectivePoint {
        let result = ProjectivePoint {
            y: self.y.add(&_rhs.y),
        };
        result
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a ProjectivePoint {
    type Output = ProjectivePoint;
    /// Multiply this `WeiestrassPoint` by a `Scalar`.
    fn mul(self, _rhs: &'b Scalar) -> ProjectivePoint {
        unimplemented!()
    }
}
