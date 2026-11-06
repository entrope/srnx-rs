//! Convert an observation file (SRNX, RINEX or Compact RINEX) to SRNX.
//!
//! `cargo run --release --example to_srnx -- [--socd] <input> <output.srnx>`

use std::process::ExitCode;

use srnx::writer::{self, ObsChunk, WriteOptions};

fn main() -> ExitCode {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut options = WriteOptions::default();
    if args.first().is_some_and(|a| a == "--socd") {
        options.obs_chunk = ObsChunk::Socd;
        args.remove(0);
    }
    let [input, output] = args.as_slice() else {
        eprintln!("usage: to_srnx [--socd] <input> <output.srnx>");
        return ExitCode::from(2);
    };
    let result = (|| -> Result<(), Box<dyn std::error::Error>> {
        let file = srnx::reader::read(&std::fs::read(input)?)?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(output)?);
        writer::write_with(&file, &options, &mut out)?;
        std::io::Write::flush(&mut out)?;
        Ok(())
    })();
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{input}: {e}");
            ExitCode::FAILURE
        }
    }
}
