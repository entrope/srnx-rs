#!/usr/bin/env python3
"""Regenerate the corpus-derived test fixtures in this directory.

Usage: make_fixtures.py CORPUS_DIR

CORPUS_DIR holds the source files named in FIXTURES below.  The
script cuts a window of epochs out of each one, then runs the reference
tools, which must be on PATH: crx2rnx and rnx2crx (RNXCMP), and rnx2srnx and
srnx_scan.  See README.md for where the source files come from.

Outputs, for each fixture NAME:
  NAME.rnx            the excerpt (RINEX 2 or 3/4), only for small fixtures
  NAME.crx            rnx2crx of the excerpt; tests decompress it for larger
                      fixtures, which saves space
  NAME.srnx           rnx2srnx of the excerpt (SOC2, default digests)
  NAME.socd.srnx      rnx2srnx --socd (revision 1.0 SOCD chunks)
  NAME.*.scan         srnx_scan of each .srnx
plus the extra variants listed per fixture, and NAME.srnx and NAME.scan for
the hand-written fixtures in SMALL.
"""

import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def label(line):
    return line[60:80].rstrip()


def split(lines):
    """Split an observation file into (version, header, records)."""
    end = next(i for i, l in enumerate(lines) if label(l) == "END OF HEADER")
    header, body = lines[: end + 1], lines[end + 1 :]
    version = float(header[0][:9])
    records = []
    i = 0
    if version >= 3:
        while i < len(body):
            if not body[i].startswith(">"):
                raise ValueError(f"line {end + 2 + i}: expected an epoch record")
            n = int(body[i][32:35])
            records.append(body[i : i + 1 + n])
            i += 1 + n
        return version, header, records

    nobs = None
    for l in header:
        if label(l) == "# / TYPES OF OBSERV" and l[:6].strip():
            nobs = int(l[:6])
    while i < len(body):
        flag = int(body[i][28:29] or "0")
        n = int(body[i][29:32])
        if flag in (0, 1, 6):
            length = math.ceil(n / 12) + n * math.ceil(nobs / 5)
        else:
            length = 1 + n
            for l in body[i + 1 : i + length]:
                if label(l) == "# / TYPES OF OBSERV" and l[:6].strip():
                    nobs = int(l[:6])
        records.append(body[i : i + length])
        i += length
    return version, header, records


def v2_prn_obs(header, records):
    """PRN / # OF OBS and # OF SATELLITES records counting `records`."""
    nobs = next(
        int(l[:6]) for l in header if label(l) == "# / TYPES OF OBSERV" and l[:6].strip()
    )
    counts = {}
    for rec in records:
        if int(rec[0][28:29] or "0") != 0:
            continue
        n = int(rec[0][29:32])
        sat_text = "".join(l[32:68] for l in rec[: math.ceil(n / 12)])
        sats = [sat_text[3 * k : 3 * k + 3].replace(" ", "0") for k in range(n)]
        obs_lines = rec[math.ceil(n / 12) :]
        per_sat = math.ceil(nobs / 5)
        for k, sat in enumerate(sats):
            text = "".join(l.ljust(80) for l in obs_lines[k * per_sat : (k + 1) * per_sat])
            fields = [text[80 * (j // 5) + 16 * (j % 5) : 80 * (j // 5) + 16 * (j % 5) + 14] for j in range(nobs)]
            row = counts.setdefault(sat, [0] * nobs)
            for j, f in enumerate(fields):
                if f.strip() and float(f) != 0.0:
                    row[j] += 1
    out = [f"{len(counts):6d}".ljust(60) + "# OF SATELLITES"]
    for sat in sorted(counts):
        row = counts[sat]
        for start in range(0, nobs, 9):
            lead = f"   {sat}" if start == 0 else "      "
            out.append((lead + "".join(f"{c:6d}" for c in row[start : start + 9])).ljust(60) + "PRN / # OF OBS")
    return out


def v402_picoseconds(header, records):
    """Relabel a RINEX 3 excerpt as 4.02 and give every epoch picoseconds."""
    header = ["     4.02" + header[0][9:]] + header[1:]
    out = []
    for k, rec in enumerate(records, start=1):
        ps = k * 37_123 % 100_000
        out.append([f"{rec[0][:56]:<56} {ps:05}"] + rec[1:])
    return header, out


def run(args, stdin=None, stdout=None):
    subprocess.run(args, stdin=stdin, stdout=stdout, cwd=HERE, check=True)


def read_lines(path):
    if path.suffix == ".crx" or path.suffix.endswith("d"):
        with open(path, "rb") as f:
            text = subprocess.run(["crx2rnx", "-"], stdin=f, capture_output=True, check=True).stdout
    else:
        text = path.read_bytes()
    lines = text.decode("ascii").split("\n")
    if lines[-1] == "":
        lines.pop()
    return [l.rstrip("\r") for l in lines]


def write_crx(rnx, out):
    """Compress with rnx2crx, fixing the date in its CRINEX PROG / DATE record
    so that regenerating the fixtures gives identical files."""
    with open(rnx, "rb") as f:
        text = subprocess.run(["rnx2crx", "-"], stdin=f, capture_output=True, check=True).stdout
    lines = text.split(b"\n")
    assert lines[1].endswith(b"CRINEX PROG / DATE"), lines[1]
    lines[1] = lines[1][:40] + b"01-Jan-26 00:00".ljust(20) + lines[1][60:]
    out.write_bytes(b"\n".join(lines))


def write_srnx(rnx, out, *options):
    run(["rnx2srnx", *options, str(rnx), str(out)])
    with open(out.with_suffix(".scan"), "wb") as f:
        run(["srnx_scan", out.name], stdout=f)


# (name, source file, first record, record count, keep .rnx?, PRN / # OF OBS?,
#  extra SRNX variants)
FIXTURES = [
    # NOAA CORS, RINEX 2.11 at 30 s, with PRN / # OF OBS records added.
    ("n302_v211_30s", "n3021830.22o", 0, 120, True, True,
     [("null", ["--chunk-digest=0", "--file-digest=0"]),
      ("blake3", ["--chunk-digest=5", "--file-digest=4"])]),
    # NOAA CORS, RINEX 2.11 at 15 s, spanning a flag-4 event (record 240).
    ("mty2_v211_15s", "mty22000.20o", 200, 80, True, False, []),
    # NOAA CORS, RINEX 2.11 at 1 Hz with 20 observation types.
    ("at01_v211_1hz", "at012000.20o", 0, 120, False, False, []),
    # Geoscience Australia, RINEX 3.04 at 1 Hz, six GNSS.
    ("balm_v304_1hz", "BALM00AUS_S_20262810000_15M_01S_MO.crx", 0, 120, False, False, []),
]


# Hand-written fixtures (not excerpts) that also get an SRNX version.
# rnx2srnx rejects rinex402_table_a4_variants, which has a truncated epoch
# and a satellite of a system without SYS / # / OBS TYPES.
SMALL = [
    "crinex_v211_features",
    "crinex_v304_features",
    "rinex211_flag1_flag6",
    "rinex402_picoseconds",
]


def main():
    corpus = Path(sys.argv[1])
    for name, source, first, count, keep_rnx, prn_obs, variants in FIXTURES:
        version, header, records = split(read_lines(corpus / source))
        records = records[first : first + count]
        if prn_obs:
            header = header[:-1] + v2_prn_obs(header, records) + header[-1:]
        rnx = HERE / f"{name}.rnx"
        rnx.write_text("\n".join(header + [l for r in records for l in r]) + "\n")
        write_crx(rnx, HERE / f"{name}.crx")
        write_srnx(rnx, HERE / f"{name}.srnx")
        write_srnx(rnx, HERE / f"{name}.socd.srnx", "--socd")
        for suffix, options in variants:
            write_srnx(rnx, HERE / f"{name}.{suffix}.srnx", *options)
        if version >= 3:
            # RINEX 4.02 picosecond epochs exercise EPO2 and SOC2 rate prediction
            # in picosecond units, and CRINEX 3.1.
            h402, r402 = v402_picoseconds(header, records)
            rnx402 = HERE / f"{name.replace('v304', 'v402')}_ps.rnx"
            rnx402.write_text("\n".join(h402 + [l for r in r402 for l in r]) + "\n")
            write_crx(rnx402, rnx402.with_suffix(".crx"))
            write_srnx(rnx402, rnx402.with_suffix(".srnx"))
            rnx402.unlink()
        if not keep_rnx:
            rnx.unlink()
    # SRNX versions of the hand-written fixtures.
    for stem in SMALL:
        write_srnx(HERE / f"{stem}.rnx", HERE / f"{stem}.srnx")


if __name__ == "__main__":
    main()
