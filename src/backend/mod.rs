//! Module for different backend implementations.
//! 
//! On this module you can find the bit-implementations
//! done for the backend maths found in Finite Fields.
//! 
pub mod u64;
#[cfg(not(any(
    feature = "u64_backend"
)))]

// A backend feature has to be choosen.
compile_error!(
    "no zerocaf backend cargo feature enabled! \
     please enable u64_backend"
);