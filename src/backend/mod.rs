pub mod u64;
#[cfg(not(any(
    feature = "u64_backend"
)))]

// A backend feature has to be choosen.
compile_error!(
    "no doppio-curve backend cargo feature enabled! \
     please enable u64_backend"
);