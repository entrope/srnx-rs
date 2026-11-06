#![doc = include_str!("../README.md")]
//!
//! ## API stability
//!
//! The public API is [`reader`], [`writer`] and [`model`].  The hidden modules are
//! format internals, public only for `srnx-diff` and the integration tests,
//! and are not covered by semver.

pub mod model;
pub mod reader;
pub mod writer;

mod chunk;
mod compat;
mod evtf;
mod leb128;
mod pack;
mod plan;
mod sate;
mod sdir;
mod soc2;
mod srnx_header;

#[doc(hidden)]
pub mod crinex;
#[doc(hidden)]
pub mod diff;
#[doc(hidden)]
pub mod digest;
#[doc(hidden)]
pub mod epoc;
#[doc(hidden)]
pub mod file;
#[doc(hidden)]
pub mod rinex2_obs;
#[doc(hidden)]
pub mod rinex3_obs;
#[doc(hidden)]
pub mod rinex_header;
#[doc(hidden)]
pub mod socd;
#[doc(hidden)]
pub mod srnx_reader;

/// Read a fixture from `tests/data`, panicking if it is missing.
#[cfg(test)]
pub(crate) fn test_data(name: &str) -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The RINEX text of fixture `stem` in `tests/data`: `stem.rnx`, or for
/// fixtures stored only in compact form, `stem.crx` decompressed.
#[cfg(test)]
pub(crate) fn test_rinex(stem: &str) -> Vec<u8> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    std::fs::read(dir.join(format!("{stem}.rnx")))
        .unwrap_or_else(|_| crinex::decompress(&test_data(&format!("{stem}.crx"))).unwrap())
}
