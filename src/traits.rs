//! Module for Public Trait implementations.


/// Trait for getting the identity element of a point type.
pub trait Identity {
    /// Returns the identity element of the curve.
    fn identity() -> Self;
}