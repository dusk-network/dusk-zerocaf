//! Contains the curve-constants needed by different algorithm implementations.

use crate::field::FieldElement;
use crate::edwards::{EdwardsPoint, CompressedEdwardsY};
use crate::ristretto::{CompressedRistretto, RistrettoPoint};

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);

/// Edwards `d` variable value = `86649/86650 (mod l)` equals:
/// `1201935917638644956968126114584555454358623906841733991436515590915937358637`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_D: FieldElement = FieldElement([939392471225133, 587442007554368, 4497154776428662, 4184267646867733, 2921744366591]);

/// Ristretto `d = -86649 (mod l)`.
pub const RISTRETTO_D: FieldElement = FieldElement([671914833248628, 3916664325105025, 1367801, 0, 17592186044416]);

/// Holds the value of the Doppio basepoint, which has been choosen as the `y`
/// coordinate `100171752`. 
/// The positive sign is choosen for it, so we leave it on it's cannonical bytes
/// encoding. 
/// UNINPLEMENTED
pub const DOPPIO_BASEPOINT_COMPRESSED: CompressedEdwardsY = 
        CompressedEdwardsY([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]); 


/// Comes from taking the `y` Twisted Edwards coordinate as: `100171752`.
pub const DOPPIO_BASEPOINT: EdwardsPoint = EdwardsPoint {
            X: FieldElement([3265320031919788, 2344618945868358, 1522956767231782, 3674566506878787, 1422874481139]),
            Y: FieldElement([100171752, 0, 0, 0, 0]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([247433686587364, 3650682313482504, 3458624897327137, 1443086282535945, 8688752063094])
        };

/// Holds the value of one of both `sqrt(-1 (mod p)) values. 
/// `SQRT_MINUS_ONE = 3034649101460298094273452163494570791663566989388331537498831373842135895065`. 
pub const SQRT_MINUS_ONE: FieldElement = FieldElement([3075585030474777, 2451921961843096, 1194333869305507, 2218299809671669, 7376823328646]); 

/// `(+)1/SQRT(a) (mod l)` equals: `4202356475871964119699734399548423449193549369991576068503119564443318355924`.
pub static INV_SQRT_A: FieldElement = FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]);

/// `(-)SQRT(a) (mod l)` equals: `4202356475871964119699734399548423449193549369991576068503119564443318355924`. 
pub static MINUS_SQRT_A: FieldElement = FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]);

/// `INV_SQRT_A_MINUS_D = 6597598851246620382811217884446123466325309731460906167807735083660259066907`.
pub const INV_SQRT_A_MINUS_D: FieldElement = FieldElement([2513811247091666, 4155030717723932, 2858009955324939, 76218730049292, 1554311650468]);

/// `SQRT_AD_MINUS_ONE = 639406726085641831161968678596870774531806627919001438194215854625195184082`.
pub const SQRT_AD_MINUS_ONE : FieldElement = FieldElement([2661703213614107, 4265233234751588, 1645589673413357, 4427380897321203, 16037874393947]);

/// The Ristretto basepoint, in `CompressedRistretto` format.
pub const RISTRETTO_BASEPOINT_COMPRESSED: CompressedRistretto =
    CompressedRistretto([45, 101, 136, 106, 139, 202, 154, 178, 57, 176, 95, 56, 189, 31, 96, 8, 216, 220, 5, 29, 234, 82, 195, 238, 188, 74, 48, 243, 219, 91, 136, 5]);

/// The Ristretto Basepoint is the same as the Curve Basepoint. 
pub const RISTRETTO_BASEPOINT: RistrettoPoint = RistrettoPoint(DOPPIO_BASEPOINT);
