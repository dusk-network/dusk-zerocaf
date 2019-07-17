//! Implementation that provides support for Montgomery Points
//! over the Doppio curve.
//! 
//! A `MontgomeryPoint` is represented as the `u-coordinate`
//! of itself in LE bytes-format.

use crate::edwards::EdwardsPoint;
use crate::field::FieldElement;

use subtle::Choice;
use subtle::ConstantTimeEq;

/// Holds the u-coordinate of a point on the Montgomery form of
/// Doppio-curve or its twist.
#[derive(Copy, Clone, Debug)]
pub struct MontgomeryPoint(pub [u8; 32]);

/// Equality of `MontgomeryPoint`s is defined mod p.
impl ConstantTimeEq for MontgomeryPoint {
    fn ct_eq(&self, other: &MontgomeryPoint) -> Choice {
        let self_fe = FieldElement::from_bytes(&self.0);
        let other_fe = FieldElement::from_bytes(&other.0);

        self_fe.ct_eq(&other_fe)
    }
}

impl Default for MontgomeryPoint {
    fn default() -> MontgomeryPoint {
        MontgomeryPoint([0u8; 32])
    }
}

impl PartialEq for MontgomeryPoint {
    fn eq(&self, other: &MontgomeryPoint) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl Eq for MontgomeryPoint {}

impl MontgomeryPoint {
    /// View this `MontgomeryPoint` as an array of bytes.
    pub fn as_bytes<'a>(&'a self) -> &'a [u8; 32] {
        &self.0
    }

    /// Convert this `MontgomeryPoint` to an array of bytes.
    pub fn to_bytes(&self) -> [u8; 32] {
        self.0
    }

    /// Attempt to convert to an `EdwardsPoint`, using the supplied
    /// choice of sign for the `EdwardsPoint`.
    pub fn to_edwards(&self, _sign: u8) -> Option<EdwardsPoint> {
        unimplemented!()
    }
}