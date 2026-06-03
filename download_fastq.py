#!/usr/bin/env python3
"""
Download FASTQ files for CCLE cell lines with RNA fusion data.

Pipeline:
  OmicsFusionFilteredSupplementary.csv  (ModelID)
      → Model.csv                        (ModelID → StrippedCellLineName)
      → filereport_read_run_PRJNA523380.tsv (StrippedCellLineName → run_accession + FTP URLs)
      → wget download
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Download FASTQ files for CCLE cell lines with fusion data."
    )
    p.add_argument(
        "--fusions",
        default="OmicsFusionFilteredSupplementary.csv",
        help="Path to OmicsFusionFilteredSupplementary.csv",
    )
    p.add_argument(
        "--model",
        default="Model.csv",
        help="Path to Model.csv",
    )
    p.add_argument(
        "--filereport",
        default="filereport_read_run_PRJNA523380.tsv",
        help="Path to filereport_read_run_PRJNA523380.tsv",
    )
    p.add_argument(
        "--outdir",
        default="fastq_downloads",
        help="Output directory for downloaded FASTQ files (default: fastq_downloads/)",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="Only download the first N matched cell lines",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be downloaded without downloading",
    )
    return p.parse_args()


def load_model_ids(fusions_path):
    df = pd.read_csv(fusions_path, usecols=["ModelID"])
    ids = df["ModelID"].dropna().unique().tolist()
    print(f"[1/4] Loaded {len(ids)} unique ModelIDs from {fusions_path}")
    return ids


def load_model_map(model_path):
    df = pd.read_csv(model_path, usecols=["ModelID", "CellLineName", "StrippedCellLineName"])
    model_map = df.set_index("ModelID")[["CellLineName", "StrippedCellLineName"]].to_dict("index")
    print(f"[2/4] Loaded {len(model_map)} entries from {model_path}")
    return model_map


def load_sra_report(report_path):
    """Build lookup: stripped_cell_line_name -> list of (run_accession, [ftp_url, ...])"""
    df = pd.read_csv(
        report_path,
        sep="\t",
        usecols=["run_accession", "sample_alias", "fastq_ftp", "library_strategy"],
    )
    rna_df = df[df["library_strategy"] == "RNA-Seq"]
    n_excluded = len(df) - len(rna_df)
    if n_excluded:
        excluded_counts = dict(df[df["library_strategy"] != "RNA-Seq"]["library_strategy"].value_counts())
        print(f"         (excluded {n_excluded} non-RNA-Seq rows: {excluded_counts})")
    lookup = {}
    for _, row in rna_df.iterrows():
        alias = str(row["sample_alias"]).strip()
        run_acc = str(row["run_accession"]).strip()
        ftp_raw = str(row["fastq_ftp"]).strip() if pd.notna(row["fastq_ftp"]) else ""
        if not ftp_raw:
            continue
        urls = [f"ftp://{u.strip()}" for u in ftp_raw.split(";") if u.strip()]
        # The alias is StrippedCellLineName + "_" + TISSUE (e.g. HCC56_LARGE_INTESTINE)
        # Extract the cell-line portion by splitting on the first underscore that separates
        # it from an all-caps tissue suffix. We store keyed by the full alias for exact
        # matching; the matching function below does the prefix lookup.
        lookup.setdefault(alias, []).append((run_acc, urls))
    print(f"[3/4] Loaded {len(lookup)} sample aliases from {report_path}")
    return lookup


def match_cell_line(stripped_name, sra_lookup):
    """Return all (run_accession, urls) for aliases that equal or start with stripped_name_."""
    matches = []
    for alias, runs in sra_lookup.items():
        if alias == stripped_name or alias.startswith(stripped_name + "_"):
            matches.extend(runs)
    return matches


def download_file(url, outdir):
    """Download a single file with wget --continue. Returns 'done', 'skipped', or 'failed'."""
    filename = url.split("/")[-1]
    dest = Path(outdir) / filename

    if dest.exists() and dest.stat().st_size > 0:
        return "skipped"

    result = subprocess.run(
        ["wget", "--continue", "--quiet", "-P", outdir, url],
        check=False,
    )
    return "done" if result.returncode == 0 else "failed"


def main():
    args = parse_args()

    for path in [args.fusions, args.model, args.filereport]:
        if not os.path.exists(path):
            sys.exit(f"Error: file not found: {path}")

    if not args.dry_run:
        os.makedirs(args.outdir, exist_ok=True)

    model_ids = load_model_ids(args.fusions)
    model_map = load_model_map(args.model)
    sra_lookup = load_sra_report(args.filereport)

    # Resolve each ModelID to matched SRA runs
    matched = []       # list of (model_id, cell_line_name, stripped, run_acc, urls)
    unmatched_ids = [] # ModelIDs with no Model.csv entry
    unmatched_cells = []  # cell lines found in Model.csv but absent from SRA report

    for mid in model_ids:
        if mid not in model_map:
            unmatched_ids.append(mid)
            continue
        info = model_map[mid]
        stripped = str(info["StrippedCellLineName"]).strip()
        cell_name = str(info["CellLineName"]).strip()
        runs = match_cell_line(stripped, sra_lookup)
        if not runs:
            unmatched_cells.append((mid, cell_name, stripped))
        else:
            for run_acc, urls in runs:
                matched.append((mid, cell_name, stripped, run_acc, urls))

    print(f"[4/4] Matched {len({r[0] for r in matched})} ModelIDs → "
          f"{len(matched)} run(s) total")

    # Write unmatched logs
    if unmatched_ids:
        with open("unmatched_model_ids.txt", "w") as f:
            f.write("\n".join(unmatched_ids) + "\n")
        print(f"  Warning: {len(unmatched_ids)} ModelIDs not found in Model.csv "
              f"→ unmatched_model_ids.txt")

    if unmatched_cells:
        with open("unmatched_cells.txt", "w") as f:
            for mid, name, stripped in unmatched_cells:
                f.write(f"{mid}\t{name}\t{stripped}\n")
        print(f"  Warning: {len(unmatched_cells)} cell lines not found in SRA report "
              f"→ unmatched_cells.txt")

    # Apply --limit: restrict to the first N unique cell lines
    if args.limit is not None:
        seen_stripped = []
        limited = []
        for row in matched:
            stripped = row[2]
            if stripped not in seen_stripped:
                seen_stripped.append(stripped)
            if len(seen_stripped) > args.limit:
                break
            limited.append(row)
        print(f"  --limit {args.limit}: restricted to {len(limited)} run(s) "
              f"across {len(seen_stripped[:args.limit])} cell line(s)")
        matched = limited

    if not matched:
        print("No runs to download. Exiting.")
        return

    # Write metadata manifest
    os.makedirs(args.outdir, exist_ok=True)
    metadata_path = os.path.join(args.outdir, "metadata.csv")
    with open(metadata_path, "w") as f:
        f.write("ModelID,CellLineName,run_accession\n")
        for mid, cell_name, stripped, run_acc, urls in matched:
            f.write(f"{mid},{cell_name},{run_acc}\n")
    print(f"  Metadata written to {metadata_path}")

    # Flatten all download tasks: (label, run_acc, url)
    tasks = []
    for mid, cell_name, stripped, run_acc, urls in matched:
        for url in urls:
            tasks.append((f"{stripped}/{run_acc}", run_acc, url))

    print(f"\n  {len(tasks)} file(s) to download across {len(matched)} run(s)\n")

    if args.dry_run:
        for label, run_acc, url in tasks:
            print(f"  [dry-run] wget -c {url} -P {args.outdir}  ({label})")
        print(f"\n{'=' * 60}")
        print(f"Done. (dry-run mode — no files were downloaded)")
        return

    errors = []
    for i, (label, run_acc, url) in enumerate(tasks, 1):
        filename = url.split("/")[-1]
        print(f"  [{i}/{len(tasks)}] {filename} ({label})")
        status = download_file(url, args.outdir)
        print(f"    → {status}")
        if status == "failed":
            errors.append((run_acc, url))

    # Summary
    print()
    print("=" * 60)
    print(f"Done.")
    print(f"  Files processed: {len(tasks)}")
    print(f"  Errors         : {len(errors)}")
    if errors:
        error_log = "download_errors.txt"
        with open(error_log, "w") as f:
            for run_acc, url in errors:
                f.write(f"{run_acc}\t{url}\n")
        print(f"  Failed URLs written to {error_log}")


if __name__ == "__main__":
    main()
