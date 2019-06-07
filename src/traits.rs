//! Module for public Trait implementations.



// Public Traits


/// Trait for getting the identity element of a point type.
pub trait Identity {
    /// Returns the identity element of the curve.
    /// Can be used as a constructor.
    fn identity() -> Self;
}