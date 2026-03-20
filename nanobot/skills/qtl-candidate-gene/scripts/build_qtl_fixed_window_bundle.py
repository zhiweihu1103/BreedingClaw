#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


GENOMEWIDE_P_DEFAULT = 0.05 / 4033947


@dataclass(frozen=True)
class Gene:
    gene_id: str
    chrom: int
    start: int
    end: int
    strand: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build LD-supported QTL summaries and plotting tables with a fixed-window fallback."
    )
    parser.add_argument("--gwas-table", required=True)
    parser.add_argument("--gene-annotation", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--p-threshold", type=float, required=True)
    parser.add_argument("--window-bp", type=int, required=True, help="Initial clustering radius and fallback half-window in bp.")
    parser.add_argument("--trait-label", default="")
    parser.add_argument("--reference-fai", default="")
    parser.add_argument("--max-loci", type=int, default=3)
    parser.add_argument("--genomewide-p", type=float, default=GENOMEWIDE_P_DEFAULT)
    parser.add_argument("--ld-table", default="")
    parser.add_argument("--geno-tfile-prefix", default="")
    parser.add_argument("--geno-bfile-prefix", default="")
    parser.add_argument("--bed-cache-prefix", default="")
    parser.add_argument("--plink-bin", default="plink")
    parser.add_argument("--plink-threads", type=int, default=1)
    parser.add_argument("--ld-r2", type=float, default=0.2)
    parser.add_argument("--local-flank-bp", type=int, default=200000)
    parser.add_argument("--min-qtl-pad-bp", type=int, default=50000)
    return parser.parse_args()


def chrom_num(raw: str) -> int:
    match = re.search(r"(\d+)$", str(raw).strip())
    if not match:
        raise ValueError(f"Cannot parse chromosome number from: {raw}")
    return int(match.group(1))


def read_gff_genes(path: Path) -> tuple[dict[int, list[Gene]], dict[int, list[int]]]:
    genes_by_chr: dict[int, list[Gene]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            try:
                chrom = chrom_num(fields[0])
            except ValueError:
                continue
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = {}
            for part in fields[8].split(";"):
                if "=" in part:
                    key, value = part.split("=", 1)
                    attrs[key] = value
            gene_id = attrs.get("ID", f"gene_{chrom}_{start}_{end}")
            genes_by_chr[chrom].append(Gene(gene_id=gene_id, chrom=chrom, start=start, end=end, strand=strand))
    starts_by_chr: dict[int, list[int]] = {}
    for chrom, genes in genes_by_chr.items():
        genes.sort(key=lambda item: (item.start, item.end, item.gene_id))
        starts_by_chr[chrom] = [gene.start for gene in genes]
    return genes_by_chr, starts_by_chr


def genes_in_region(
    genes_by_chr: dict[int, list[Gene]],
    starts_by_chr: dict[int, list[int]],
    chrom: int,
    start: int,
    end: int,
) -> list[Gene]:
    genes = genes_by_chr.get(chrom, [])
    starts = starts_by_chr.get(chrom, [])
    if not genes:
        return []
    index = bisect_left(starts, start)
    while index > 0 and genes[index - 1].end >= start:
        index -= 1
    out = []
    while index < len(genes) and genes[index].start <= end:
        gene = genes[index]
        if gene.end >= start:
            out.append(gene)
        index += 1
    return out


def nearest_genes(
    genes_by_chr: dict[int, list[Gene]],
    starts_by_chr: dict[int, list[int]],
    chrom: int,
    pos: int,
    limit: int = 3,
) -> list[Gene]:
    genes = genes_by_chr.get(chrom, [])
    starts = starts_by_chr.get(chrom, [])
    if not genes:
        return []
    index = bisect_left(starts, pos)
    candidates: list[tuple[int, Gene]] = []
    for j in range(max(0, index - 10), min(len(genes), index + 10)):
        gene = genes[j]
        distance = 0 if gene.start <= pos <= gene.end else min(abs(pos - gene.start), abs(pos - gene.end))
        candidates.append((distance, gene))
    candidates.sort(key=lambda item: (item[0], item[1].start, item[1].gene_id))
    return [gene for _, gene in candidates[:limit]]


def load_chrom_sizes(path: str) -> dict[int, int]:
    if not path:
        return {}
    sizes: dict[int, int] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                chrom = chrom_num(fields[0])
            except ValueError:
                continue
            sizes[chrom] = int(fields[1])
    return sizes


def detect_schema(fieldnames: list[str]) -> str:
    names = set(fieldnames)
    if {"trait", "chrom", "pos", "snp_id", "pvalue"}.issubset(names):
        return "gwas_summary"
    if {"SNP", "CHR", "BP", "P"}.issubset(names):
        return "pmap"
    raise ValueError(f"Unsupported GWAS table schema: {fieldnames}")


def standardize_row(row: dict[str, str], schema: str, trait_label: str) -> dict[str, object]:
    if schema == "gwas_summary":
        return {
            "trait": row["trait"] or trait_label,
            "chrom": chrom_num(row["chrom"]),
            "bp": int(float(row["pos"])),
            "snp_id": row["snp_id"],
            "pvalue": float(row["pvalue"]),
            "beta": row.get("beta", ""),
            "stderr": row.get("stderr", ""),
        }
    return {
        "trait": trait_label,
        "chrom": chrom_num(row["CHR"]),
        "bp": int(float(row["BP"])),
        "snp_id": row["SNP"],
        "pvalue": float(row["P"]),
        "beta": "",
        "stderr": "",
    }


def build_loci(
    sig_hits: list[dict[str, object]],
    cluster_bp: int,
    genomewide_p: float,
    max_loci: int,
) -> list[dict[str, object]]:
    loci: list[dict[str, object]] = []
    by_chrom: dict[int, list[dict[str, object]]] = defaultdict(list)
    for hit in sig_hits:
        by_chrom[int(hit["chrom"])].append(hit)
    for chrom, hits in by_chrom.items():
        work = sorted(hits, key=lambda item: (float(item["pvalue"]), int(item["bp"]), str(item["snp_id"])))
        while work:
            lead = work[0]
            lo = max(1, int(lead["bp"]) - cluster_bp)
            hi = int(lead["bp"]) + cluster_bp
            members = [row for row in work if lo <= int(row["bp"]) <= hi]
            work = [row for row in work if not (lo <= int(row["bp"]) <= hi)]
            loci.append(
                {
                    "chrom": chrom,
                    "lead_snp": str(lead["snp_id"]),
                    "lead_bp": int(lead["bp"]),
                    "lead_p": float(lead["pvalue"]),
                    "coarse_start": min(int(row["bp"]) for row in members),
                    "coarse_end": max(int(row["bp"]) for row in members),
                    "member_count": len(members),
                    "gw_member_count": sum(1 for row in members if float(row["pvalue"]) <= genomewide_p),
                    "members": members,
                }
            )
    loci.sort(key=lambda item: (float(item["lead_p"]), int(item["chrom"]), int(item["lead_bp"])))
    for index, locus in enumerate(loci[:max_loci], start=1):
        locus["locus_id"] = f"L{index:03d}"
        locus["is_genomewide"] = bool(float(locus["lead_p"]) <= genomewide_p)
    return loci[:max_loci]


def run(cmd: list[str]) -> None:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError("Command failed:\n" + " ".join(cmd) + "\n" + proc.stdout)


def prefix_has_bed(prefix: str) -> bool:
    if not prefix:
        return False
    base = Path(prefix)
    return all(Path(str(base) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def ensure_bed_source(args: argparse.Namespace, outdir: Path) -> str:
    if prefix_has_bed(args.geno_bfile_prefix):
        return args.geno_bfile_prefix

    if not args.geno_tfile_prefix:
        return ""

    if args.bed_cache_prefix:
        cache_prefix = Path(args.bed_cache_prefix)
    else:
        cache_prefix = outdir.parent.parent / "_ld_cache" / Path(args.geno_tfile_prefix).name

    cache_prefix.parent.mkdir(parents=True, exist_ok=True)
    if prefix_has_bed(str(cache_prefix)):
        return str(cache_prefix)

    run(
        [
            args.plink_bin,
            "--tfile",
            args.geno_tfile_prefix,
            "--allow-extra-chr",
            "--allow-no-sex",
            "--make-bed",
            "--out",
            str(cache_prefix),
        ]
    )
    return str(cache_prefix)


def sanitize_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def parse_ld_lines(path: Path) -> list[dict[str, object]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = None
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if header is None:
                header = fields
                continue
            record = dict(zip(header, fields))
            snp = record.get("SNP_B") or record.get("SNP")
            bp = record.get("BP_B") or record.get("BP")
            chrom = record.get("CHR_B") or record.get("CHR")
            r2 = record.get("R2")
            if not snp or bp is None or chrom is None or r2 is None:
                continue
            try:
                rows.append(
                    {
                        "snp": snp,
                        "chrom": chrom_num(chrom),
                        "bp": int(float(bp)),
                        "r2": float(r2),
                    }
                )
            except ValueError:
                continue
    return rows


def compute_ld_rows(
    args: argparse.Namespace,
    bed_prefix: str,
    chrom: int,
    lead_snp: str,
    start: int,
    end: int,
    tmp_dir: Path,
) -> list[dict[str, object]]:
    if not bed_prefix:
        return []
    win_kb = max(1, math.ceil((end - start) / 1000))
    out_prefix = tmp_dir / f"{sanitize_name(lead_snp)}.{chrom}.{start}.{end}"
    run(
        [
            args.plink_bin,
            "--bfile",
            bed_prefix,
            "--allow-extra-chr",
            "--allow-no-sex",
            "--r2",
            "--ld-snp",
            lead_snp,
            "--chr",
            str(chrom),
            "--from-bp",
            str(start),
            "--to-bp",
            str(end),
            "--ld-window",
            "99999",
            "--ld-window-kb",
            str(win_kb),
            "--ld-window-r2",
            "0",
            "--threads",
            str(args.plink_threads),
            "--out",
            str(out_prefix),
        ]
    )
    return parse_ld_lines(Path(str(out_prefix) + ".ld"))


def load_ld_table(path: str) -> dict[str, list[dict[str, object]]]:
    if not path or not Path(path).exists():
        return {}
    table: dict[str, list[dict[str, object]]] = defaultdict(list)
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}
        for row in reader:
            lead = row.get("lead_snp") or row.get("lead_snp_id") or row.get("SNP_A")
            snp = row.get("snp") or row.get("snp_id") or row.get("SNP_B")
            chrom = row.get("chrom") or row.get("CHR_B") or row.get("CHR")
            bp = row.get("bp") or row.get("BP_B") or row.get("BP")
            r2 = row.get("r2") or row.get("R2")
            if not lead or not snp or chrom is None or bp is None or r2 is None:
                continue
            try:
                table[str(lead)].append(
                    {
                        "snp": str(snp),
                        "chrom": chrom_num(str(chrom)),
                        "bp": int(float(bp)),
                        "r2": float(r2),
                    }
                )
            except ValueError:
                continue
    return dict(table)


def apply_ld_boundaries(
    locus: dict[str, object],
    ld_rows: list[dict[str, object]],
    ld_r2: float,
    local_flank_bp: int,
    min_qtl_pad_bp: int,
) -> tuple[dict[str, object], dict[str, float]]:
    q1 = int(locus["coarse_start"])
    q2 = int(locus["coarse_end"])
    ld_lookup: dict[str, float] = {str(row["snp"]): float(row["r2"]) for row in ld_rows}
    ld_lookup[str(locus["lead_snp"])] = 1.0

    high = [
        row
        for row in ld_rows
        if str(row["snp"]) != str(locus["lead_snp"]) and float(row["r2"]) >= ld_r2
    ]
    if high:
        q1 = min(q1, min(int(row["bp"]) for row in high))
        q2 = max(q2, max(int(row["bp"]) for row in high))
        qtl_type = "ld_supported"
    else:
        qtl_type = "fixed_window_fallback"

    if q1 == q2:
        q1 = max(1, q1 - min_qtl_pad_bp)
        q2 = q2 + min_qtl_pad_bp
    elif q2 - q1 < min_qtl_pad_bp:
        mid = (q1 + q2) // 2
        half = max(1, min_qtl_pad_bp // 2)
        q1 = max(1, mid - half)
        q2 = mid + half

    locus["qtl_start"] = q1
    locus["qtl_end"] = q2
    locus["panel_start"] = max(1, q1 - local_flank_bp)
    locus["panel_end"] = q2 + local_flank_bp
    locus["qtl_type"] = qtl_type
    return locus, ld_lookup


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def slice_rows_by_region(
    global_by_chr: dict[int, list[dict[str, object]]],
    bp_index_by_chr: dict[int, list[int]],
    chrom: int,
    start: int,
    end: int,
) -> list[dict[str, object]]:
    rows = global_by_chr.get(chrom, [])
    bp_index = bp_index_by_chr.get(chrom, [])
    if not rows:
        return []
    left = bisect_left(bp_index, start)
    right = bisect_right(bp_index, end)
    return rows[left:right]


def main() -> None:
    args = parse_args()
    gwas_path = Path(args.gwas_table)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    trait_label = args.trait_label or gwas_path.name.replace(".gwas_summary.tsv", "").replace(".all.QQ.pmap.txt", "")
    genes_by_chr, starts_by_chr = read_gff_genes(Path(args.gene_annotation))
    chrom_sizes = load_chrom_sizes(args.reference_fai)
    ld_table_by_lead = load_ld_table(args.ld_table)
    bed_prefix = ensure_bed_source(args, outdir)
    tmp_dir = outdir / "_ld_tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    significant_hits: list[dict[str, object]] = []
    global_rows: list[dict[str, object]] = []
    global_by_chr: dict[int, list[dict[str, object]]] = defaultdict(list)

    with gwas_path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Empty GWAS table: {gwas_path}")
        schema = detect_schema(reader.fieldnames)
        for row in reader:
            std = standardize_row(row, schema=schema, trait_label=trait_label)
            chrom = int(std["chrom"])
            bp = int(std["bp"])
            pvalue = max(float(std["pvalue"]), float.fromhex("0x1.0p-1022"))
            chrom_sizes[chrom] = max(chrom_sizes.get(chrom, 0), bp)
            item = {
                "snp_id": std["snp_id"],
                "chrom": chrom,
                "bp": bp,
                "p_value": pvalue,
                "mlog10_p": -math.log10(pvalue),
            }
            global_rows.append(item)
            global_by_chr[chrom].append(item)
            if pvalue <= args.p_threshold:
                significant_hits.append(std)

    if not global_rows:
        raise SystemExit(f"No rows found in GWAS table: {gwas_path}")

    for chrom in list(global_by_chr):
        global_by_chr[chrom].sort(key=lambda row: (int(row["bp"]), str(row["snp_id"])))
    bp_index_by_chr = {
        chrom: [int(row["bp"]) for row in rows]
        for chrom, rows in global_by_chr.items()
    }

    loci = build_loci(significant_hits, args.window_bp, args.genomewide_p, args.max_loci)
    if not loci:
        raise SystemExit(f"No significant loci found at p <= {args.p_threshold}")

    qtl_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    qtl_plot_rows: list[dict[str, object]] = []
    local_rows: list[dict[str, object]] = []
    highlight_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    added_genes: set[str] = set()
    status_rows: list[dict[str, object]] = []

    try:
        for locus in loci:
            ld_start = max(1, int(locus["lead_bp"]) - args.window_bp)
            ld_end = int(locus["lead_bp"]) + args.window_bp
            ld_rows = ld_table_by_lead.get(str(locus["lead_snp"]), [])
            ld_source = "precomputed" if ld_rows else "none"
            if not ld_rows and bed_prefix:
                try:
                    ld_rows = compute_ld_rows(
                        args=args,
                        bed_prefix=bed_prefix,
                        chrom=int(locus["chrom"]),
                        lead_snp=str(locus["lead_snp"]),
                        start=ld_start,
                        end=ld_end,
                        tmp_dir=tmp_dir,
                    )
                    ld_source = "plink"
                except Exception as exc:
                    ld_rows = []
                    ld_source = f"failed:{exc}"

            locus, ld_lookup = apply_ld_boundaries(
                locus=locus,
                ld_rows=ld_rows,
                ld_r2=args.ld_r2,
                local_flank_bp=args.local_flank_bp,
                min_qtl_pad_bp=args.min_qtl_pad_bp,
            )

            genes = genes_in_region(
                genes_by_chr,
                starts_by_chr,
                int(locus["chrom"]),
                int(locus["qtl_start"]),
                int(locus["qtl_end"]),
            )
            if not genes:
                genes = nearest_genes(
                    genes_by_chr,
                    starts_by_chr,
                    int(locus["chrom"]),
                    int(locus["lead_bp"]),
                    limit=10,
                )
            genes = genes[:10]
            top_gene_ids = [gene.gene_id for gene in genes[:3]]

            qtl_rows.append(
                {
                    "trait": trait_label,
                    "locus_id": locus["locus_id"],
                    "chrom": locus["chrom"],
                    "lead_snp": locus["lead_snp"],
                    "lead_bp": locus["lead_bp"],
                    "lead_p": locus["lead_p"],
                    "qtl_type": locus["qtl_type"],
                    "coarse_start": locus["coarse_start"],
                    "coarse_end": locus["coarse_end"],
                    "qtl_start": locus["qtl_start"],
                    "qtl_end": locus["qtl_end"],
                    "qtl_span_bp": int(locus["qtl_end"]) - int(locus["qtl_start"]) + 1,
                    "member_count": locus["member_count"],
                    "gw_member_count": locus["gw_member_count"],
                    "top_candidate_genes": "; ".join(top_gene_ids),
                }
            )

            qtl_plot_rows.append(
                {
                    "locus_id": locus["locus_id"],
                    "chrom": locus["chrom"],
                    "lead_snp_id": locus["lead_snp"],
                    "lead_bp": locus["lead_bp"],
                    "lead_p": locus["lead_p"],
                    "qtl_start": locus["qtl_start"],
                    "qtl_end": locus["qtl_end"],
                    "panel_start": locus["panel_start"],
                    "panel_end": locus["panel_end"],
                    "qtl_type": locus["qtl_type"],
                }
            )

            for rank, gene in enumerate(genes, start=1):
                distance = 0 if gene.start <= int(locus["lead_bp"]) <= gene.end else min(abs(int(locus["lead_bp"]) - gene.start), abs(int(locus["lead_bp"]) - gene.end))
                candidate_rows.append(
                    {
                        "trait": trait_label,
                        "locus_id": locus["locus_id"],
                        "rank": rank,
                        "gene_id": gene.gene_id,
                        "chrom": gene.chrom,
                        "gene_start": gene.start,
                        "gene_end": gene.end,
                        "lead_snp": locus["lead_snp"],
                        "lead_bp": locus["lead_bp"],
                        "lead_p": locus["lead_p"],
                        "distance_bp": distance,
                        "in_qtl": int(gene.start <= int(locus["qtl_end"]) and gene.end >= int(locus["qtl_start"])),
                    }
                )
                if rank <= 3:
                    highlight_rows.append(
                        {
                            "locus_id": locus["locus_id"],
                            "gene_id": gene.gene_id,
                            "rank": rank,
                            "label": gene.gene_id,
                        }
                    )
                if gene.gene_id not in added_genes:
                    gene_rows.append(
                        {
                            "gene_id": gene.gene_id,
                            "chrom": gene.chrom,
                            "start_bp": gene.start,
                            "end_bp": gene.end,
                            "strand": gene.strand,
                        }
                    )
                    added_genes.add(gene.gene_id)

            local_subset = slice_rows_by_region(
                global_by_chr=global_by_chr,
                bp_index_by_chr=bp_index_by_chr,
                chrom=int(locus["chrom"]),
                start=int(locus["panel_start"]),
                end=int(locus["panel_end"]),
            )
            for row in local_subset:
                local_rows.append(
                    {
                        "locus_id": locus["locus_id"],
                        "snp_id": row["snp_id"],
                        "chrom": row["chrom"],
                        "bp": row["bp"],
                        "p_value": row["p_value"],
                        "mlog10_p": row["mlog10_p"],
                        "r2": ld_lookup.get(str(row["snp_id"]), 0.0),
                        "is_lead": int(str(row["snp_id"]) == str(locus["lead_snp"])),
                    }
                )

            status_rows.append(
                {
                    "locus_id": locus["locus_id"],
                    "lead_snp": locus["lead_snp"],
                    "qtl_type": locus["qtl_type"],
                    "ld_source": ld_source,
                    "ld_row_count": len(ld_rows),
                }
            )
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    plot_dir = outdir / "plot_tables"
    write_tsv(
        outdir / f"{trait_label}.qtl_summary.tsv",
        [
            "trait",
            "locus_id",
            "chrom",
            "lead_snp",
            "lead_bp",
            "lead_p",
            "qtl_type",
            "coarse_start",
            "coarse_end",
            "qtl_start",
            "qtl_end",
            "qtl_span_bp",
            "member_count",
            "gw_member_count",
            "top_candidate_genes",
        ],
        qtl_rows,
    )
    write_tsv(
        outdir / f"{trait_label}.top10_candidate_genes.tsv",
        [
            "trait",
            "locus_id",
            "rank",
            "gene_id",
            "chrom",
            "gene_start",
            "gene_end",
            "lead_snp",
            "lead_bp",
            "lead_p",
            "distance_bp",
            "in_qtl",
        ],
        candidate_rows,
    )
    write_tsv(
        plot_dir / "chrom_sizes.tsv",
        ["chrom", "length_bp"],
        [{"chrom": chrom, "length_bp": chrom_sizes[chrom]} for chrom in sorted(chrom_sizes)],
    )
    write_tsv(
        plot_dir / "global_manhattan.tsv",
        ["snp_id", "chrom", "bp", "p_value", "mlog10_p"],
        global_rows,
    )
    write_tsv(
        plot_dir / "qtl_regions.tsv",
        [
            "locus_id",
            "chrom",
            "lead_snp_id",
            "lead_bp",
            "lead_p",
            "qtl_start",
            "qtl_end",
            "panel_start",
            "panel_end",
            "qtl_type",
        ],
        qtl_plot_rows,
    )
    write_tsv(
        plot_dir / "local_manhattan.tsv",
        ["locus_id", "snp_id", "chrom", "bp", "p_value", "mlog10_p", "r2", "is_lead"],
        local_rows,
    )
    write_tsv(
        plot_dir / "gene_models.tsv",
        ["gene_id", "chrom", "start_bp", "end_bp", "strand"],
        gene_rows,
    )
    write_tsv(
        plot_dir / "highlight_genes.tsv",
        ["locus_id", "gene_id", "rank", "label"],
        highlight_rows,
    )
    (outdir / f"{trait_label}.status.json").write_text(
        json.dumps(
            {
                "trait": trait_label,
                "selected_loci": len(qtl_rows),
                "candidate_genes": len(candidate_rows),
                "ld_mode": "plink" if bed_prefix else ("precomputed" if ld_table_by_lead else "fixed_window_fallback"),
                "loci": status_rows,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
