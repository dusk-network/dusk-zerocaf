//! This module contains backend-specific constant values as the 64-bit limbs of curve constants.

use crate::backend::u64::field::FieldElement;
use crate::backend::u64::scalar::Scalar;


/// `L` is the order of base point for Doppio, in this case it is equivalent to 2^249 + 14490550575682688738086195780655237219
pub const L: Scalar = Scalar([1129677152307299, 1363544697812651, 714439, 0, 2199023255552]);

/// `(L - 1) / 2` used to check positiveness of a `FieldElement` on the Decaf paper. 
pub(crate) const POS_RANGE: FieldElement = FieldElement([2587757230352886, 4210131976237760, 683900, 0, 8796093022208]);

/// `(L - 1) / 2` used to check positiveness of a `FieldElement` on the Decaf paper. 
pub(crate) const POS_RANGE: FieldElement = FieldElement([2587757230352886, 4210131976237760, 683900, 0, 8796093022208]);

/// Scalar-LFACTOR is the value that satisfies the equation: `L * LFACTOR = -1 (mod 2^52)`
/// In this case, `LFACTOR` is the one used for the Montgomery Reduction algorithm,
/// implemented on Scalar Arithmetics module.
pub const LFACTOR: u64 = 1331240223835829;

/// Montgomery modulus defined for Scalar arithmetics, `R^2 = (2^260)^2 % L`
pub const RR: Scalar = Scalar([137682194168839, 3209056245311277, 1480926248458276, 2533620989757837, 1314911199310]);

/// `FIELD_L` is the order of the Prime field for Doppio, n this case it is equivalent to 2^252 + 27742317777372353535851937790883648493`
pub const FIELD_L: FieldElement = FieldElement([671914833335277, 3916664325105025, 1367801, 0, 17592186044416]);

/// Montgomery modulus defined for FieldElement arithmetics, `R^2 = (2^260)^2 % FIELD_L`
pub const RR_FIELD: FieldElement = FieldElement([2764609938444603, 3768881411696287, 1616719297148420, 1087343033131391, 10175238647962]);

/// FieldElement-LFACTOR is the value that satisfies the equation: `L * LFACTOR = -1 (mod 2^52)`
/// In this case, `LFACTOR` is the one used for the Montgomery Reduction algorithm,
/// implemented on FieldElement Arithmetics module.
pub const LFACTOR_FIELD: u64 = 1439961107955227;

/// Montgomery modulus defined for FieldElements on `inverse()` functions scope. 
/// It is used for the Montgomery Mul operation that takes place on the `Inversion
// operation`. It's defined as: `R^2 = (2^253)^2 % L`
pub const INV_RR: FieldElement = FieldElement([2210115751650724, 3809421927348411, 2357176729341513, 3420097284349172, 7483527818736]);
