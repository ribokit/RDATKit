"""Command-line entry points for rdat_kit.

Two subcommands:

  rdat_kit validate  <file.rdat> [file2.rdat ...]
      Parse an RDAT file and run RDATFile.validate(). Exit code:
        0 — file parsed and validate() returned no warnings
        1 — parse failure
        2 — validate() returned at least one warning (file parsed OK
            but has shape/annotation/seqpos issues worth fixing)

  rdat_kit to_md     <file.rdat>  [--rmdb-id ID]
      Emit a Jekyll/RMDB front-matter `.md` stub on stdout, populated
      from the RDAT (sequence, structure, offset, annotations,
      comments). Optionally accept --rmdb-id to override the default
      derived from the filename. Designed to feed
      https://github.com/DasLab/rmdb.github.io contributions.
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from datetime import date as _date

from .handler import RDATFile


# ── validate ────────────────────────────────────────────────────────── #

def cmd_validate(args: argparse.Namespace) -> int:
    """Returns process exit code."""
    rc = 0
    for path in args.files:
        rdat = RDATFile()
        try:
            with open(path) as f:
                rdat.load(f)
        except Exception as e:
            print(f"{path}: PARSE FAIL — {e}", file=sys.stderr)
            rc = max(rc, 1)
            continue

        messages = rdat.validate() or []
        if messages:
            print(f"{path}: {len(messages)} warning(s)")
            for m in messages:
                print(f"  {m}")
            rc = max(rc, 2)
        else:
            n_constructs = len(rdat.constructs)
            n_data = sum(len(c.data) for c in rdat.constructs.values())
            print(f"{path}: OK   ({n_constructs} construct(s), "
                  f"{n_data} data row(s), RDAT_VERSION={rdat.version})")
    return rc


# ── to_md ───────────────────────────────────────────────────────────── #

# Canonical RMDB ID looks like  <PREFIX>_<CHEM>_<NNNN>  e.g. 4WY16S_DCP_0000.
# Files may also carry a version suffix: 4WY16S_DCP_0000_2.rdat → strip the _2.
_RMDB_ID_RE = re.compile(r"^([A-Z0-9]+_[A-Z0-9]+_\d{4})(_\d+)?$")


def _slug_from_path(path: str) -> str:
    """Infer the RMDB_ID from a filename like 'OK7BLIB_2A3_0000.rdat'
    or '4WY16S_DCP_0000_2.rdat' (versioned)."""
    base = os.path.basename(path)
    base = re.sub(r"\.rdat$", "", base, flags=re.IGNORECASE)
    m = _RMDB_ID_RE.match(base)
    return m.group(1) if m else base


def _yaml_quote(s) -> str:
    if s is None or s == "":
        return '""'
    s = str(s).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{s}"'


def _yaml_list(items) -> str:
    return "[" + ", ".join(_yaml_quote(i) for i in items) + "]"


def _yaml_block(s) -> str:
    if not s:
        return '""'
    text = str(s).replace("\r\n", "\n").replace("\t", "    ")
    text = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]", "", text)
    lines = [ln.strip() for ln in text.split("\n")]
    while lines and not lines[-1]:
        lines.pop()
    if not lines:
        return '""'
    return "|\n" + "\n".join("    " + ln for ln in lines)


def cmd_to_md(args: argparse.Namespace) -> int:
    rdat = RDATFile()
    try:
        with open(args.file) as f:
            rdat.load(f)
    except Exception as e:
        print(f"{args.file}: PARSE FAIL — {e}", file=sys.stderr)
        return 1

    if not rdat.constructs:
        print(f"{args.file}: no constructs found", file=sys.stderr)
        return 1

    # Use first construct as the canonical view.
    name, c = next(iter(rdat.constructs.items()))
    rmdb_id = args.rmdb_id or _slug_from_path(args.file)

    # Aggregate annotations (top-level + construct-level), de-duped.
    annot = {}
    for src in (rdat.annotations or {}, getattr(c, "annotations", {}) or {}):
        for k, v in src.items():
            if isinstance(v, list):
                annot.setdefault(k, [])
                annot[k] = list(dict.fromkeys(annot[k] + v))
            else:
                annot.setdefault(k, v)

    construct_count = len(rdat.constructs)
    data_points = sum(len(getattr(d, "values", [])) for c2 in rdat.constructs.values()
                                                    for d in c2.data)

    fm = ["---"]
    fm.append(f"rmdb_id: {_yaml_quote(rmdb_id)}")
    fm.append(f"permalink: /detail/{rmdb_id}/")
    fm.append(f"name: {_yaml_quote(name)}")
    fm.append(f"category: \"General\"     # General | RNA_Puzzles | Eterna")
    fm.append(f"date: {_date.today().isoformat()}        # update with the original deposition date if known")
    fm.append(f"sequence: {_yaml_quote(c.sequence)}")
    fm.append(f"structure: {_yaml_quote(getattr(c, 'structure', '') or '.' * len(c.sequence))}")
    fm.append(f"offset: {getattr(c, 'offset', 0)}")
    fm.append(f"construct_count: {construct_count}")
    fm.append(f"data_points: {data_points}")
    if rdat.comments:
        fm.append(f"comments: " + _yaml_block(rdat.comments))
    if annot:
        fm.append("annotation:")
        for k, v in sorted(annot.items()):
            if isinstance(v, list):
                fm.append(f"  {k}: {_yaml_list(v)}")
            else:
                fm.append(f"  {k}: {_yaml_quote(v)}")
    # Citation skeleton — contributor fills in
    fm.append("citation:")
    fm.append('  authors: ""')
    fm.append('  title:   ""')
    fm.append('  journal: ""')
    fm.append('  year:    ""')
    fm.append('  doi:     ""')
    fm.append('  pubmed:  ""')
    fm.append(f"thumbnail: /assets/thumbnails/{rmdb_id}.png")
    fm.append(f"rdat:      https://github.com/DasLab/rmdb.github.io/releases/"
              f"download/data-general/{rmdb_id}.rdat   "
              f"# update release tag (data-eterna|puzzle|riboswitches|rna-structures|general)")
    fm.append("---")
    fm.append("")
    print("\n".join(fm))
    return 0


# ── main dispatcher ────────────────────────────────────────────────── #

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="rdat_kit",
        description="RDATKit CLI — validate RDAT files and emit RMDB entry stubs.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    v = sub.add_parser("validate", help="Parse + schema-check one or more RDAT files.")
    v.add_argument("files", nargs="+", help="*.rdat path(s)")
    v.set_defaults(func=cmd_validate)

    t = sub.add_parser(
        "to_md",
        help="Emit a Jekyll/RMDB front-matter .md stub for the given RDAT.",
    )
    t.add_argument("file", help="*.rdat path")
    t.add_argument("--rmdb-id", default=None,
                   help="Override the RMDB_ID (default: derived from filename)")
    t.set_defaults(func=cmd_to_md)

    th = sub.add_parser(
        "thumbnail",
        help="Render a reactivity heatmap PNG from an RDAT.",
    )
    th.add_argument("file", help="*.rdat path")
    th.add_argument("--out", default=".",
                    help="Output directory (default: .)")
    th.add_argument("--rmdb-id", default=None,
                    help="Override the RMDB_ID (default: derived from filename)")
    th.add_argument("--max-rows", type=int, default=1000,
                    help="Truncate to first N data rows (default 1000)")
    th.add_argument("--short-px", type=int, default=440,
                    help="Short-axis target pixel count (default 440)")
    th.add_argument("--long-max-px", type=int, default=2400,
                    help="Long-axis cap in pixels (default 2400)")
    th.set_defaults(func=cmd_thumbnail)

    return p


# ── thumbnail ──────────────────────────────────────────────────────── #
#
# Render a reactivity heatmap from an RDAT. matplotlib + numpy are
# imported lazily so users who only need validate / to_md don't pay
# for the heavy dependencies. If they're missing we emit a clear
# install hint.

def cmd_thumbnail(args: argparse.Namespace) -> int:
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        matplotlib.rcParams['font.family'] = 'Helvetica'
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(
            "rdat_kit thumbnail requires matplotlib + numpy.\n"
            "  pip install matplotlib numpy\n"
            f"(missing: {e})",
            file=sys.stderr,
        )
        return 1

    rdat = RDATFile()
    try:
        with open(args.file) as f:
            rdat.load(f)
    except Exception as e:
        print(f"{args.file}: PARSE FAIL — {e}", file=sys.stderr)
        return 1
    if not rdat.constructs:
        print(f"{args.file}: no constructs found", file=sys.stderr)
        return 1

    name, c = next(iter(rdat.constructs.items()))
    rmdb_id = args.rmdb_id or _slug_from_path(args.file)

    total_rows = len(c.data)
    use_rows = min(total_rows, args.max_rows)
    values = np.array([d.values for d in c.data[:use_rows]], dtype=float)
    n_rows, n_cols = values.shape

    print(f"{rmdb_id}: {use_rows}/{total_rows} rows × {n_cols} cols, "
          f"value range [{np.nanmin(values):.2f}…{np.nanmax(values):.2f}]"
          + (f"; truncated to first {use_rows:,}"
             if total_rows > use_rows else ""))

    # Robust contrast: 0 to 90th percentile of finite values
    finite = values[np.isfinite(values)]
    vmax = max(0.05, float(np.percentile(finite, 90))) if finite.size else 1.0
    arr = np.where(np.isfinite(values), values, 0.0)

    # Legacy aspect logic (matches Server_RMDB media.py):
    # - eterna libraries OR <3 rows → 'auto' (stretch into default figure)
    # - otherwise → 'equal' (square cells)
    is_eterna = "ETERNA" in rmdb_id.upper()
    legacy_aspect = "auto" if (is_eterna or n_rows < 3) else "equal"

    dpi = 100
    if legacy_aspect == "auto":
        fig, ax = plt.subplots(dpi=dpi, facecolor="white")
    else:
        short_cells = min(n_rows, n_cols)
        long_cells = max(n_rows, n_cols)
        px = min(args.short_px / short_cells,
                 args.long_max_px / long_cells)
        w = max(int(round(n_cols * px)), 1)
        h = max(int(round(n_rows * px)), 1)
        fig, ax = plt.subplots(figsize=(w / dpi, h / dpi),
                               dpi=dpi, facecolor="white")

    ax.imshow(arr, cmap="Greys", vmin=0, vmax=vmax,
              aspect=legacy_aspect, interpolation="nearest")
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{rmdb_id}.png")
    plt.savefig(out_path, dpi=dpi, facecolor="white",
                bbox_inches="tight",
                pil_kwargs={"optimize": True, "compress_level": 9})
    plt.close(fig)
    print(f"  wrote {out_path} ({os.path.getsize(out_path)/1024:.1f} KB)")
    return 0


# ── main ────────────────────────────────────────────────────────────── #

def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
