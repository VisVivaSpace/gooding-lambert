fn main() {
    #[cfg(feature = "gooding-ffi")]
    cc::Build::new().file("csrc/lamb.c").compile("lamb");
}
