//! Contains the curve-constants needed by different algorithm implementations.

use crate::field::FieldElement;
use crate::edwards::{EdwardsPoint, CompressedEdwardsY};
use crate::ristretto::{CompressedRistretto, RistrettoPoint};

#[cfg(feature = "u64_backend")]
pub use crate::backend::u64::constants::*;

/// Holds the value of the Curve basepoint, which has been constructed
/// from taking `y-coodrinate = 3/5 (mod l)`.
/// The positive sign is choosen for it, so we leave it on it's cannonical bytes
/// encoding. 
pub const BASEPOINT_COMPRESSED: CompressedEdwardsY = 
    CompressedEdwardsY([194, 24, 45, 158, 220, 161, 164, 1, 
                        231, 42, 46, 200, 184, 98, 31, 166, 
                        153, 153, 153, 153, 153, 153, 153, 153, 
                        153, 153, 153, 153, 153, 153, 153, 9]);


/// Ristretto Basepoint on compressed format. 
pub const RISTRETTO_BASEPOINT_COMPRESSED: CompressedRistretto =
    CompressedRistretto([2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);