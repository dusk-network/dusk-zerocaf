//! This module contains backend-specific constant values as the 64-bit limbs of curve constants.

use crate::backend::u64::field::FieldElement;
use crate::backend::u64::scalar::Scalar;

/// `L` is the order of base point for Doppio, i.e. `2^249 - 15145038707218910765482344729778085401`.
pub const L: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551]);

/// Scalar-LFACTOR is the value that satisfies the equation: `L * LFACTOR = -1 (mod 2^52)`
/// In this case, `LFACTOR` is the one used for the Montgomery Reduction algorithm,
/// implemented on Scalar Arithmetics module.
pub const LFACTOR: u64 = 547593343082025;

/// Montgomery modulus defined for Scalar arithmetics, `R^2 = (2^260)^2 % L`
pub const RR: Scalar = Scalar([1682248870925813, 4078880264703668, 2289123149127681, 4169238435752846, 2104335664921]);

/// `FIELD_L` is the order of the Prime field for Doppio, i.e. 2^252 + 27742317777372353535851937790883648493`
pub const FIELD_L: FieldElement = FieldElement([671914833335277, 3916664325105025, 1367801, 0, 17592186044416]);

/// Montgomery modulus defined for FieldElement arithmetics, `R^2 = (2^260)^2 % FIELD_L`
pub const RR_FIELD: FieldElement = FieldElement([2764609938444603, 3768881411696287, 1616719297148420, 1087343033131391, 10175238647962]);

/// FieldElement-LFACTOR is the value that satisfies the equation: `L * LFACTOR = -1 (mod 2^52)`
/// In this case, `LFACTOR` is the one used for the Montgomery Reduction algorithm,
/// implemented on FieldElement Arithmetics module.
pub const LFACTOR_FIELD: u64 = 1439961107955227;

/// Montgomery modulus defined for FieldElements on `inverse()` function's scope. 
/// It's used for the Montgomery Mul operation that takes place on the `Inversion
/// operation`. It's defined as: `R^2 = (2^253)^2 % L`.
pub const INV_RR: FieldElement = FieldElement([2210115751650724, 3809421927348411, 2357176729341513, 3420097284349172, 7483527818736]);
