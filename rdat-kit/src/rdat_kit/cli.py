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

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
