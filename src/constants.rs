//! Contains the curve-constants needed by different algorithm implementations.

use crate::field::FieldElement;
use crate::ristretto::CompressedRistretto;

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);

/// Edwards `d` variable value = `86649/86650 (mod l)` equals:
/// `1201935917638644956968126114584555454358623906841733991436515590915937358637`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_D: FieldElement = FieldElement([939392471225133, 587442007554368, 4497154776428662, 4184267646867733, 2921744366591]);

/// Ristretto `d = -121665/121666 (mod l) = 3939827596983564740704173820156434201364153519171569545177282234949263247860`.
pub const RISTRETTO_D: FieldElement = FieldElement([3932549999839732, 251630550069643, 1998572741534892, 2238362972951308, 9577190362565]);

/// `1/SQRT(a) (mod l)` equals: ``.
pub static INV_SQRT_A: FieldElement = FieldElement([3075585030474777, 2451921961843096, 1194333869305507, 2218299809671669, 7376823328646]);

/// `-SQRT(-1) (mod l)` equals: ``. 
pub static MINUS_SQRT_A: FieldElement = FieldElement([3075585030474777, 2451921961843096, 1194333869305507, 2218299809671669, 7376823328646]);

/// INV_SQRT_A_MINUS_D
pub const INV_SQRT_A_MINUS_D: FieldElement = FieldElement([1, 0, 0, 0, 0]);

/// SQRT_AD_MINUS_ONE
pub const SQRT_AD_MINUS_ONE : FieldElement = FieldElement([1, 0 , 0, 0, 0]);

/// The Ristretto basepoint, in `CompressedRistretto` format.
pub const RISTRETTO_BASEPOINT_COMPRESSED: CompressedRistretto =
    CompressedRistretto([0xe2, 0xf2, 0xae, 0x0a, 0x6a, 0xbc, 0x4e, 0x71,
                         0xa8, 0x84, 0xa9, 0x61, 0xc5, 0x00, 0x51, 0x5f,
                         0x58, 0xe3, 0x0b, 0x6a, 0xa5, 0x82, 0xdd, 0x8d,
                         0xb6, 0xa6, 0x59, 0x45, 0xe0, 0x8d, 0x2d, 0x76]);
