# Test fixtures

The integration tests in `crates/srnx/tests/` and the unit tests in
`crates/srnx/src/` read only the files in this directory. They need no
external tools.

The repository's licence (the Unlicense) covers the hand-written fixtures and
the scripts here. It does not cover the excerpts of third-party data listed
below; those are under their sources' terms.

## Hand-written fixtures

| Stem | Purpose |
| --- | --- |
| `crinex_v211_features` | RINEX 2.11 and CRINEX 1.0 corner cases |
| `crinex_v304_features` | RINEX 3.04 and CRINEX 3.0 corner cases |
| `rinex211_flag1_flag6` | Power-failure and cycle-slip records, zero values |
| `rinex402_picoseconds` | RINEX 4.02 picosecond epochs, CRINEX 3.1, zero values |
| `rinex402_table_a4_variants` | RINEX 4.02 epoch-record variants |

These were written for this project. Each `.crx` was made from the `.rnx` by
`RNX2CRX` (RNXCMP 4.2.0), and each `.crx.rnx` is the output of `CRX2RNX`
4.2.0 on the `.crx`. `make_fixtures.py` makes `.srnx` and `.scan` files for
`crinex_v211_features`, `crinex_v304_features`, `rinex211_flag1_flag6` and
`rinex402_picoseconds`.  `rinex402_table_a4_variants` has none: it has a
truncated epoch and a BeiDou satellite without `SYS / # / OBS TYPES`, which
`rnx2srnx` rejects, and it has no `.crx` either.

## Excerpts of station data

| Stem | Source file | Station, date | Excerpt |
| --- | --- | --- | --- |
| `n302_v211_30s` | `n3021830.22o` | N302, 2022-07-02 | First 120 epochs (1 hour at 30 s) |
| `mty2_v211_15s` | `mty22000.20o` | MTY2, 2020-07-18 | Records 200–279 (00:50 to 01:09:30 at 15 s), including a flag-4 event |
| `at01_v211_1hz` | `at012000.20o` | AT01, 2020-07-18 | First 120 epochs (2 minutes at 1 s) |
| `balm_v304_1hz` | `BALM00AUS_S_20262810000_15M_01S_MO.crx` | BALM, 2026-10-08 | First 120 epochs (2 minutes at 1 s) |
| `balm_v402_1hz_ps` | as `balm_v304_1hz` | as `balm_v304_1hz` | Relabelled RINEX 4.02, with synthetic picosecond digits |

Changes from the source files:

- Every excerpt keeps the full header but only a window of epochs, so
  header records such as `TIME OF LAST OBS` describe the whole source file.
- `n302_v211_30s` has `# OF SATELLITES` and `PRN / # OF OBS` records added,
  counting the excerpt.
- `balm_v402_1hz_ps` gives epoch *k* (from 1) `k × 37123 mod 100000`
  picoseconds. These times are made up; they exist to exercise picosecond
  handling.
- `at01_v211_1hz`, `balm_v304_1hz` and `balm_v402_1hz_ps` are stored only as
  Compact RINEX (`.crx`), to save space; tests decompress them. The
  decompressed text follows `CRX2RNX`'s formatting (for example, no leading
  zero before a decimal point), so it differs slightly from the source.

### NOAA CORS Network (NCN)

`n302_v211_30s`, `mty2_v211_15s` and `at01_v211_1hz` are derived from data
of the NOAA Continuously Operating Reference Stations (CORS) Network (NCN),
managed by NOAA's National Geodetic Survey.
NOAA states that its data are open to the public and can be used as
desired, and requests attribution.

These files are modified (truncated and reformatted), so they are not
original, unaltered NOAA data. Their use here does not imply endorsement by
or affiliation with NOAA.

### Geoscience Australia

`balm_v304_1hz` and `balm_v402_1hz_ps` are derived from data in the
Geoscience Australia GNSS data archive, which Geoscience Australia makes
available under the [Creative Commons Attribution 3.0 Australia
licence](https://creativecommons.org/licenses/by/3.0/au/).

These files are adaptations: excerpted, reformatted, and in the case of
`balm_v402_1hz_ps`, relabelled and given synthetic epoch times.

## Generated files

Every `.srnx` file was made by `rnx2srnx` from
[SRNX](https://github.com/entrope/srnx/), and each `.scan` file is
`srnx_scan` output for the `.srnx` of the same name.
For each excerpt, `NAME.srnx` uses the default options (revision 1.1, or 1.2
for picosecond epochs) and `NAME.socd.srnx` uses `--socd` (revision 1.0).
`n302_v211_30s.null.srnx` and `n302_v211_30s.blake3.srnx` use non-default
chunk and file digests.

To regenerate the excerpts and every generated file, put the source files
named above in one directory, put `crx2rnx`, `rnx2crx`, `rnx2srnx` and
`srnx_scan` on `PATH`, and run:

```sh
python3 make_fixtures.py /path/to/source/files
```

## Checking against other data

`crates/srnx/tests/corpus.rs` checks a directory of your own observation
files against the reference tools:

```sh
SRNX_CORPUS=/path/to/files cargo test --release --test corpus -- --nocapture
```
