//! This module contains curve-specific constant values, 
//! needed to either testing or implementing algorithms.

use crate::scalar::Scalar;
use crate::field::FieldElement;

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub (crate) static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);