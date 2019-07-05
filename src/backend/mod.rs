//! Contains the different backend implementations: `u64` and further comming ones.
//! 
//! On this module you can find the different implementations
//! done for Finite Fields mathematical-backends.

/// The u64 backend contains the implementation of all of the 
/// mathematical base eg. Arithmetics over Finite Fields with 
/// a design specially thought for 64-bit architectures.
pub mod u64;
#[cfg(not(any(
    feature = "u64_backend"
)))]

// A backend feature has to be choosen.
compile_error!(
    "no zerocaf backend cargo feature enabled! \
     please enable u64_backend"
);