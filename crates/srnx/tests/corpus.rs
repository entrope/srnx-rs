//! Bring-your-own-corpus tests, run only when `SRNX_CORPUS` names a
//! directory of observation files.  Each file is checked against the
//! reference tools, which must be on `PATH`: `rnx2srnx`, `srnx2rnx` and
//! `srnx_scan` (SRNX), and `crx2rnx` (RNXCMP; `$CRX2RNX` overrides it).
//!
//! - RINEX and Compact RINEX files are converted with `rnx2srnx`, both to
//!   SOC2 (revision 1.1, or 1.2 for picosecond epochs) and with `--socd`
//!   (revision 1.0).  Our reading of each SRNX file must match our reading
//!   of the source, and of `srnx2rnx`'s decoding of the SRNX file.
//! - Our Compact RINEX decompression must match `crx2rnx` byte for byte.
//! - SRNX files must match `srnx2rnx`'s decoding and `srnx_scan`'s
//!   statistics.
//! - Every file is written as SRNX with our writer, both with `SOC2` and
//!   `SOCD` chunks.  Our reading of the result must match the source, and
//!   so must `rnxcmp`'s comparison of the two files (`$RNXCMP` overrides
//!   it; files `rnxcmp` cannot read are not compared), and `srnx2rnx` must
//!   decode it.
//!
//! Files in other formats are skipped, as are files the reference tools
//! reject (if we reject them too).  Every file is checked, and the test
//! fails at the end if any check failed.  Large files take a while;
//! `--release` helps.

mod common;

use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

use srnx::crinex;
use srnx::model::ObsFile;
use srnx::reader::{self, Format};
use srnx::writer::{self, ObsChunk, WriteOptions};

/// A temporary file, removed when dropped.
struct TempFile(PathBuf);

impl TempFile {
    fn new(extension: &str) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let n = NEXT.fetch_add(1, Ordering::Relaxed);
        Self(std::env::temp_dir().join(format!(
            "srnx-corpus-{}-{n}.{extension}",
            std::process::id()
        )))
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Run a reference tool.  `Ok(None)` if it ran and failed (for example by
/// rejecting the file); `Err` if it could not be run at all.
fn tool(
    program: impl AsRef<std::ffi::OsStr>,
    args: &[&Path],
    stdin: Option<&Path>,
) -> Result<Option<Vec<u8>>, String> {
    let program = program.as_ref();
    let mut command = Command::new(program);
    command.args(args).stderr(Stdio::piped());
    if let Some(path) = stdin {
        command.stdin(std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?);
    }
    let output = command
        .output()
        .map_err(|e| format!("cannot run {}: {e}", program.to_string_lossy()))?;
    Ok(output.status.success().then_some(output.stdout))
}

fn read(bytes: &[u8], label: &str) -> ObsFile {
    reader::read(bytes).unwrap_or_else(|e| panic!("{label}: {e}"))
}

/// Our reading of `srnx2rnx`'s decoding of `srnx` must match `expected`.
fn check_srnx2rnx(srnx: &Path, expected: &ObsFile, label: &str) -> Result<(), String> {
    let decoded = TempFile::new("rnx");
    if tool("srnx2rnx", &[srnx, &decoded.0], None)?.is_none() {
        return Err(format!("{label}: srnx2rnx failed"));
    }
    let from_tool = read(&std::fs::read(&decoded.0).unwrap(), "srnx2rnx output");
    common::assert_same(expected, &from_tool, label, "srnx2rnx output");
    Ok(())
}

/// Convert a RINEX or Compact RINEX file with `rnx2srnx` and compare.
fn check_rinex(path: &Path, source: &ObsFile) -> Result<(), String> {
    for options in [&[][..], &["--socd"][..]] {
        let srnx = TempFile::new("srnx");
        let mut args: Vec<&Path> = options.iter().map(Path::new).collect();
        args.extend([path, srnx.0.as_path()]);
        if tool("rnx2srnx", &args, None)?.is_none() {
            eprintln!("skip {}: rnx2srnx {options:?} rejected it", path.display());
            return Ok(());
        }
        let label = format!("rnx2srnx {options:?} output");
        let converted = read(&std::fs::read(&srnx.0).unwrap(), &label);
        common::assert_same(source, &converted, &path.display().to_string(), &label);
        check_srnx2rnx(&srnx.0, &converted, &label)?;
    }
    Ok(())
}

/// Write `source` (read from `path`) with our writer and compare.
fn check_writer(path: &Path, source: &ObsFile) -> Result<(), String> {
    let rnxcmp = std::env::var_os("RNXCMP").unwrap_or_else(|| "rnxcmp".into());
    // rnxcmp exits nonzero for files it cannot read, as for differences.
    let comparable = tool(&rnxcmp, &[path, path], None)?.is_some();
    if !comparable {
        eprintln!(
            "{}: rnxcmp cannot read it; not comparing with rnxcmp",
            path.display()
        );
    }
    for obs_chunk in [ObsChunk::Soc2, ObsChunk::Socd] {
        let label = format!("our {obs_chunk:?} output");
        let mut options = WriteOptions::default();
        options.obs_chunk = obs_chunk;
        let mut bytes = Vec::new();
        writer::write_with(source, &options, &mut bytes)
            .map_err(|e| format!("{}: {label}: {e}", path.display()))?;
        let written = read(&bytes, &label);
        common::assert_same(source, &written, &path.display().to_string(), &label);
        let srnx = TempFile::new("srnx");
        std::fs::write(&srnx.0, &bytes).unwrap();
        if comparable && tool(&rnxcmp, &[path, &srnx.0], None)?.is_none() {
            return Err(format!(
                "{}: rnxcmp finds {label} different",
                path.display()
            ));
        }
        check_srnx2rnx(&srnx.0, &written, &label)?;
    }
    Ok(())
}

fn check_crinex(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let program = std::env::var_os("CRX2RNX").unwrap_or_else(|| "crx2rnx".into());
    let Some(want) = tool(program, &[Path::new("-")], Some(path))? else {
        eprintln!("skip {}: crx2rnx rejected it", path.display());
        return Ok(());
    };
    let got = crinex::decompress(bytes).map_err(|e| format!("{}: {e}", path.display()))?;
    if got != want {
        return Err(format!(
            "{}: decompression differs from crx2rnx",
            path.display()
        ));
    }
    Ok(())
}

fn check_srnx(path: &Path, bytes: &[u8], file: &ObsFile) -> Result<(), String> {
    check_srnx2rnx(path, file, &path.display().to_string())?;
    let Some(scan) = tool("srnx_scan", &[path], None)? else {
        return Err(format!("{}: srnx_scan failed", path.display()));
    };
    common::assert_matches_scan(
        bytes,
        &String::from_utf8_lossy(&scan),
        &path.display().to_string(),
    );
    Ok(())
}

/// Whether the reference tool for `format` accepts the file at `path`.
fn reference_accepts(path: &Path, format: Format) -> Result<bool, String> {
    let output = TempFile::new("out");
    let program = if format == Format::Srnx {
        "srnx2rnx"
    } else {
        "rnx2srnx"
    };
    Ok(tool(program, &[path, &output.0], None)?.is_some())
}

fn check_file(path: &Path) -> Result<bool, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let Some(format) = Format::detect(&bytes) else {
        eprintln!("skip {}: not an observation file", path.display());
        return Ok(false);
    };
    let file = match reader::read(&bytes) {
        Ok(file) => file,
        Err(e) if !reference_accepts(path, format)? => {
            eprintln!(
                "skip {}: we and the reference tool both reject it ({e})",
                path.display()
            );
            return Ok(false);
        }
        Err(e) => {
            return Err(format!(
                "{}: {e} (the reference tool accepts it)",
                path.display()
            ))
        }
    };
    check_reader(path, &bytes, format, &file)?;
    check_writer(path, &file)?;
    Ok(true)
}

fn check_reader(path: &Path, bytes: &[u8], format: Format, file: &ObsFile) -> Result<(), String> {
    match format {
        Format::Srnx => check_srnx(path, bytes, file)?,
        Format::Crinex => {
            check_crinex(path, bytes)?;
            check_rinex(path, file)?;
        }
        _ => check_rinex(path, file)?,
    }
    Ok(())
}

#[test]
fn corpus() {
    let Some(dir) = std::env::var_os("SRNX_CORPUS") else {
        eprintln!("skip corpus test: set SRNX_CORPUS to a directory of observation files");
        return;
    };
    let mut paths: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", Path::new(&dir).display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_file())
        .collect();
    paths.sort();

    let mut checked = 0;
    let mut failures = Vec::new();
    for path in &paths {
        match panic::catch_unwind(AssertUnwindSafe(|| check_file(path))) {
            Ok(Ok(true)) => checked += 1,
            Ok(Ok(false)) => {}
            Ok(Err(message)) => failures.push(message),
            Err(payload) => failures.push(
                payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| payload.downcast_ref::<&str>().map(|s| (*s).to_owned()))
                    .unwrap_or_else(|| format!("{}: panicked", path.display())),
            ),
        }
    }
    eprintln!("checked {checked} files, {} failed", failures.len());
    assert!(failures.is_empty(), "{}", failures.join("\n\n"));
}
