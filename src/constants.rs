//! Contains the curve-constants needed by different algorithm implementations.

use crate::field::FieldElement;
use crate::edwards::{EdwardsPoint, CompressedEdwardsY};
use crate::ristretto::{CompressedRistretto, RistrettoPoint};

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);

/// Edwards `d` variable value = `-126296/126297 (mod l)` equals:
/// `951605751702391019481481818669129158712512026257330939079110344917983315091`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_D: FieldElement = FieldElement([3304133203739795, 2446467598308289, 1534112949566882, 2032729967918914, 2313225441931]);

/// Holds the value of the Doppio basepoint, which has been choosen as the `y`
/// coordinate `100171752`. 
/// The positive sign is choosen for it, so we leave it on it's cannonical bytes
/// encoding. 
/// UNINPLEMENTED
pub const DOPPIO_BASEPOINT_COMPRESSED: CompressedEdwardsY = 
        CompressedEdwardsY([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]); 


/// Comes from taking the `x` Twisted Edwards coordinate as: `4`.
pub const DOPPIO_BASEPOINT: EdwardsPoint = EdwardsPoint {
            X: FieldElement([3461346045645288, 972387018214097, 4435206378704739, 3440261531857766, 17448879694677]),
            Y: FieldElement([4083500257631147, 3305309464452868, 1629709588575767, 4306635512831061, 7662705179761]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([212017923156762, 2250920437320259, 2648273971913259, 353963282208220, 3122143411880])
        };

/// 4Coset of a RistrettoPoint. 
pub(crate) const FOUR_COSET_GROUP: [EdwardsPoint; 3] = 
    [
        EdwardsPoint {
            X: FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]),
            Y: FieldElement([0, 0, 0, 0, 0]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([0, 0, 0, 0, 0])
        }, 

        EdwardsPoint {
            X: FieldElement([0, 0, 0, 0, 0]),
            Y: FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([0, 0, 0, 0, 0])
        },
        
        EdwardsPoint {
            X: FieldElement([3075585030474777, 2451921961843096, 1194333869305507, 2218299809671669, 7376823328646]),
            Y: FieldElement([0, 0, 0, 0, 0]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([0, 0, 0, 0, 0])
        },
    ];

/// Holds the value of one of both `sqrt(-1 (mod p)) values. 
/// `SQRT_MINUS_ONE = 3034649101460298094273452163494570791663566989388331537498831373842135895065`. 
pub const SQRT_MINUS_ONE: FieldElement = FieldElement([3075585030474777, 2451921961843096, 1194333869305507, 2218299809671669, 7376823328646]); 

/// `(+)1/SQRT(a) (mod l)` equals: `4202356475871964119699734399548423449193549369991576068503119564443318355924`.
pub static INV_SQRT_A: FieldElement = FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]);

/// `(-)SQRT(a) (mod l)` equals: `4202356475871964119699734399548423449193549369991576068503119564443318355924`. 
pub static MINUS_SQRT_A: FieldElement = FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]);

/// `INV_SQRT_A_MINUS_D = 482283834104289360917429750399313974390948281833312135312952165682596457149`.
pub const INV_SQRT_A_MINUS_D: FieldElement = FieldElement([550050132044477, 3953042081665262, 2971403105229349, 212915494370164, 1172367057772]);

/// `SQRT_AD_MINUS_ONE = `.
pub const SQRT_AD_MINUS_ONE : FieldElement = FieldElement([3601277882726560, 1817821323014817, 1726005090908779, 2111284621343800, 648674458156]);

/// The Ristretto basepoint, in `CompressedRistretto` format.
pub const RISTRETTO_BASEPOINT_COMPRESSED: CompressedRistretto =
    CompressedRistretto([45, 101, 136, 106, 139, 202, 154, 178, 57, 176, 95, 56, 189, 31, 96, 8, 216, 220, 5, 29, 234, 82, 195, 238, 188, 74, 48, 243, 219, 91, 136, 5]);

/// The Ristretto Basepoint is the same as the Curve Basepoint. 
pub const RISTRETTO_BASEPOINT: RistrettoPoint = RistrettoPoint(DOPPIO_BASEPOINT);
