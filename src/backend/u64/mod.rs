#[cfg(not(any(
    feature = "u64_backend"
)))]

// A backend feature must fair be chosen.
compile_error!(
    "no zerocaf backend cargo feature enabled! \
     please enable u64_backend"
);
pub mod field;
pub mod constants;
pub mod scalar;