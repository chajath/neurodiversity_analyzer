#!/usr/bin/env python3
"""
23andMe Raw Data Analyzer — Neurodiversity-Associated SNP Checker
(Quantified Edition with Odds Ratios, Confidence Intervals, and Polygenic Scores)

Parses a downloaded 23andMe raw data file and checks for SNPs
associated with various neurodevelopmental/neuropsychiatric conditions
from published genome-wide association studies (GWAS).

Supported conditions:
    - Autism Spectrum Disorder (ASD)
    - Attention Deficit Hyperactivity Disorder (ADHD)
    - Dyslexia
    - Obsessive-Compulsive Disorder (OCD)
    - Tourette Syndrome
    - Bipolar Disorder
    - Schizophrenia (shared neurodevelopmental overlap)
    - Cognitive Ability / Giftedness
    - Sensory Processing Sensitivity / Overexcitability
    - Psychomotor Intensity

IMPORTANT DISCLAIMER:
    This tool is for EDUCATIONAL and INFORMATIONAL purposes only.
    It is NOT a diagnostic tool. These conditions are highly polygenic
    and multifactorial. No set of SNPs can predict or diagnose any of
    them. Always consult a qualified genetic counselor or medical
    professional for interpretation.

Usage:
    python neurodiversity_analyzer.py <23andme_raw_data.txt> [--condition CONDITION]

    Without --condition, all panels are analyzed.
"""

import argparse
import math
import os
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class SNPMarker:
    """A SNP associated with a condition from published research."""

    rsid: str
    gene: str
    chromosome: str
    risk_allele: str
    odds_ratio: float  # Published OR for risk allele
    ci_lower: float  # 95% CI lower bound
    ci_upper: float  # 95% CI upper bound
    risk_allele_freq: float  # Risk allele frequency in general population
    description: str
    source: str


@dataclass
class ConditionPanel:
    """A set of SNPs associated with a neurodevelopmental condition."""

    name: str
    short_name: str
    prevalence: float  # Population prevalence (e.g. 0.01 = 1%)
    snps: list[SNPMarker]
    notes: list[str]


# ---------------------------------------------------------------------------
# Data loading from YAML files
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).parent / "data"


def _load_panels() -> dict[str, ConditionPanel]:
    """Load all condition panels from data/panels.yaml."""
    with open(DATA_DIR / "panels.yaml", "r") as f:
        raw = yaml.safe_load(f)
    panels = {}
    for key, pdata in raw.items():
        snps = [
            SNPMarker(
                rsid=s["rsid"],
                gene=s["gene"],
                chromosome=s["chromosome"],
                risk_allele=s["risk_allele"],
                odds_ratio=s["odds_ratio"],
                ci_lower=s["ci_lower"],
                ci_upper=s["ci_upper"],
                risk_allele_freq=s["risk_allele_freq"],
                description=s["description"],
                source=s["source"],
            )
            for s in pdata["snps"]
        ]
        panels[key] = ConditionPanel(
            name=pdata["name"],
            short_name=pdata["short_name"],
            prevalence=pdata["prevalence"],
            snps=snps,
            notes=pdata["notes"],
        )
    return panels


def _load_interactions() -> list[tuple[str, str, str, str, float, str, str]]:
    """Load gene-gene interactions from data/interactions.yaml."""
    with open(DATA_DIR / "interactions.yaml", "r") as f:
        raw = yaml.safe_load(f)
    return [
        (
            i["rsid1"],
            i["rsid2"],
            i["gene1"],
            i["gene2"],
            i["interaction_or"],
            i["description"],
            i["source"],
        )
        for i in raw
    ]


def _load_books() -> dict[str, dict]:
    """Load book composite indices from data/books.yaml."""
    with open(DATA_DIR / "books.yaml", "r") as f:
        raw = yaml.safe_load(f)
    books = {}
    for name, bdata in raw.items():
        components = [
            (c["panel"], c["weight"], c["dimension"])
            for c in bdata.get("components", [])
        ]
        cross_rg = {}
        for pair_str, val in bdata.get("cross_rg", {}).items():
            k1, k2 = pair_str.split(",", 1)
            cross_rg[(k1.strip(), k2.strip())] = val
        books[name] = {
            "author": bdata.get("author", ""),
            "description": bdata.get("description", ""),
            "components": components,
            "cross_rg": cross_rg,
        }
    return books


ALL_PANELS = _load_panels()
GENE_INTERACTIONS = _load_interactions()
COMPOSITE_INDICES = _load_books()


# ---------------------------------------------------------------------------
# Shared genes across panels
# ---------------------------------------------------------------------------
def find_shared_snps(panels: list[ConditionPanel]) -> dict[str, list[str]]:
    """Find SNP rsids that appear in multiple panels -> list of condition names."""
    rsid_conditions: dict[str, list[str]] = {}
    for panel in panels:
        for snp in panel.snps:
            rsid_conditions.setdefault(snp.rsid, []).append(panel.short_name)
    return {rsid: conds for rsid, conds in rsid_conditions.items() if len(conds) > 1}


# ---------------------------------------------------------------------------
# File parsing
# ---------------------------------------------------------------------------
def _parse_lines(lines, snp_data: dict[str, tuple[str, str, str]]) -> None:
    """Parse lines from a 23andMe raw data file into snp_data dict."""
    for line in lines:
        if isinstance(line, bytes):
            line = line.decode("utf-8", errors="replace")
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            parts = line.split(",")
        if len(parts) >= 4:
            rsid, chrom, position, genotype = (
                parts[0].strip(),
                parts[1].strip(),
                parts[2].strip(),
                parts[3].strip(),
            )
            if rsid.startswith("rs") or rsid.startswith("i"):
                snp_data[rsid] = (chrom, position, genotype)


def parse_23andme_file(filepath: str) -> dict[str, tuple[str, str, str]]:
    """Parse a 23andMe raw data file (.txt or .zip). Returns rsid -> (chromosome, position, genotype)."""
    snp_data: dict[str, tuple[str, str, str]] = {}

    if not os.path.isfile(filepath):
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    if zipfile.is_zipfile(filepath):
        with zipfile.ZipFile(filepath, "r") as zf:
            txt_files = [
                n
                for n in zf.namelist()
                if n.endswith(".txt") and not n.startswith("__MACOSX")
            ]
            if not txt_files:
                print(
                    "Error: No .txt file found inside the zip archive.", file=sys.stderr
                )
                sys.exit(1)
            # Use the first (or largest) .txt file
            txt_name = txt_files[0]
            if len(txt_files) > 1:
                txt_name = max(txt_files, key=lambda n: zf.getinfo(n).file_size)
            print(f"  Extracting from zip: {txt_name}")
            with zf.open(txt_name) as f:
                _parse_lines(f, snp_data)
    else:
        with open(filepath, "r", encoding="utf-8") as f:
            _parse_lines(f, snp_data)

    return snp_data


def count_risk_alleles(genotype: str, risk_allele: str) -> int:
    """Count copies of the risk allele in the genotype. Returns -1 for no-call."""
    if not genotype or genotype in ("--", ".."):
        return -1
    return sum(1 for a in genotype if a == risk_allele)


# ---------------------------------------------------------------------------
# Polygenic score calculations
# ---------------------------------------------------------------------------
def compute_per_snp_or(snp: SNPMarker, n_risk: int) -> float:
    """
    Compute the genotype-level odds ratio for a given number of risk alleles.
    Assumes a multiplicative (log-additive) model:
        0 copies -> OR = 1.0
        1 copy   -> OR = published OR
        2 copies -> OR = published OR^2
    """
    if n_risk <= 0:
        return 1.0
    return snp.odds_ratio**n_risk


def compute_polygenic_score(
    snps_results: list[tuple[SNPMarker, int]],
) -> tuple[float, float, float, float, float]:
    """
    Compute a log-additive polygenic score from observed genotypes.

    Returns (combined_or, ci_lower, ci_upper, log_or_score, combined_se).

    The combined OR is the product of per-SNP ORs (multiplicative model).
    The 95% CI is computed by propagating the per-SNP log(OR) standard errors:
        SE(log OR) ~ (log(CI_upper) - log(CI_lower)) / (2 * 1.96)
        Combined SE = sqrt(sum(SE_i^2 * n_risk_i^2))
    """
    total_log_or = 0.0
    total_var = 0.0

    for snp, n_risk in snps_results:
        if n_risk <= 0:
            continue
        log_or = math.log(snp.odds_ratio) * n_risk
        se = (math.log(snp.ci_upper) - math.log(snp.ci_lower)) / (2 * 1.96)
        total_log_or += log_or
        total_var += (se * n_risk) ** 2

    combined_se = math.sqrt(total_var) if total_var > 0 else 0
    combined_or = math.exp(total_log_or)
    ci_low = math.exp(total_log_or - 1.96 * combined_se)
    ci_high = math.exp(total_log_or + 1.96 * combined_se)

    return combined_or, ci_low, ci_high, total_log_or, combined_se


def compute_expected_score_stats(
    snps: list[SNPMarker],
) -> tuple[float, float]:
    """
    Compute the expected (mean) and std dev of the polygenic log-OR score
    in the general population, assuming Hardy-Weinberg equilibrium.

    For each SNP with risk allele freq p:
        E[n_risk] = 2p
        Var[n_risk] = 2p(1-p)
        E[log_or_contribution] = 2p * log(OR)
        Var[log_or_contribution] = 2p(1-p) * log(OR)^2
    """
    mean_score = 0.0
    var_score = 0.0

    for snp in snps:
        p = snp.risk_allele_freq
        lor = math.log(snp.odds_ratio)
        mean_score += 2 * p * lor
        var_score += 2 * p * (1 - p) * lor * lor

    return mean_score, math.sqrt(var_score) if var_score > 0 else 0.0


# ---------------------------------------------------------------------------
# Known gene-gene interaction effects (epistasis)
# When two risk genes co-occur, the combined effect may be > multiplicative.
# Format: (rsid1, rsid2, gene1, gene2, interaction_or, description, source)
# ---------------------------------------------------------------------------
def check_interactions(
    snp_data: dict[str, tuple[str, str, str]],
    panel_snps: list[SNPMarker],
) -> list[tuple[str, str, str, str, float, int, int, float, str]]:
    """
    Check for known gene-gene interactions among genotyped SNPs.
    Returns list of (rsid1, rsid2, gene1, gene2, interaction_or,
    n_risk1, n_risk2, effective_or_boost, description).
    Only reports interactions where BOTH SNPs carry >= 1 risk allele.
    """
    # Build set of panel rsids for relevance
    panel_rsids = {s.rsid for s in panel_snps}
    found = []

    for rsid1, rsid2, gene1, gene2, inter_or, desc, _source in GENE_INTERACTIONS:
        # At least one SNP should be in the panel
        if rsid1 not in panel_rsids and rsid2 not in panel_rsids:
            continue
        if rsid1 not in snp_data or rsid2 not in snp_data:
            continue
        _, _, geno1 = snp_data[rsid1]
        _, _, geno2 = snp_data[rsid2]

        # Find risk alleles from panel or interaction definition
        risk1 = risk2 = None
        for s in panel_snps:
            if s.rsid == rsid1:
                risk1 = s.risk_allele
            if s.rsid == rsid2:
                risk2 = s.risk_allele
        # For SNPs not in this panel, check ALL_PANELS
        if risk1 is None:
            for p in ALL_PANELS.values():
                for s in p.snps:
                    if s.rsid == rsid1:
                        risk1 = s.risk_allele
                        break
        if risk2 is None:
            for p in ALL_PANELS.values():
                for s in p.snps:
                    if s.rsid == rsid2:
                        risk2 = s.risk_allele
                        break
        if risk1 is None or risk2 is None:
            continue

        n1 = count_risk_alleles(geno1, risk1)
        n2 = count_risk_alleles(geno2, risk2)
        if n1 > 0 and n2 > 0:
            # Effective boost = interaction OR scaled by copies
            boost = inter_or ** (min(n1, n2) / 2)
            found.append((rsid1, rsid2, gene1, gene2, inter_or, n1, n2, boost, desc))

    return found


def score_to_percentile(z: float) -> float:
    """Convert a z-score to a percentile using the standard normal CDF approximation."""
    if z < -8:
        return 0.0
    if z > 8:
        return 100.0
    a1, a2, a3, a4, a5 = (
        0.254829592,
        -0.284496736,
        1.421413741,
        -1.453152027,
        1.061405429,
    )
    p_const = 0.3275911
    sign = 1 if z >= 0 else -1
    x = abs(z) / math.sqrt(2)
    t = 1.0 / (1.0 + p_const * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)
    return round((0.5 * (1.0 + sign * y)) * 100, 1)


# ---------------------------------------------------------------------------
# Analysis & reporting
# ---------------------------------------------------------------------------
def analyze_panel(
    panel: ConditionPanel,
    snp_data: dict[str, tuple[str, str, str]],
    shared_snps: dict[str, list[str]],
) -> None:
    """Analyze and print results for one condition panel."""

    print()
    print("=" * 100)
    print(f"  {panel.name.upper()} — SNP ANALYSIS")
    print(f"  Population prevalence: ~{panel.prevalence * 100:.1f}%")
    print("=" * 100)
    print()

    found_count = 0
    risk_count = 0
    total_risk_alleles = 0
    total_possible = 0
    results: list[tuple[SNPMarker, str, int, float]] = []
    score_inputs: list[tuple[SNPMarker, int]] = []

    for snp in panel.snps:
        if snp.rsid in snp_data:
            found_count += 1
            _, _, genotype = snp_data[snp.rsid]
            n_risk = count_risk_alleles(genotype, snp.risk_allele)

            if n_risk >= 0:
                geno_or = compute_per_snp_or(snp, n_risk)
                score_inputs.append((snp, n_risk))
            else:
                geno_or = 1.0

            if n_risk < 0:
                pass
            elif n_risk == 0:
                total_possible += 2
            else:
                risk_count += 1
                total_risk_alleles += n_risk
                total_possible += 2

            results.append((snp, genotype, n_risk, geno_or))

    if found_count == 0:
        print(f"  No {panel.short_name}-associated SNPs were found in your data.")
        print("  This may be due to chip version differences.\n")
        return

    # Compute column widths dynamically
    gene_w = max(len("Gene"), max(len(s.gene) for s, _, _, _ in results)) + 2

    # Print SNP table
    print(
        f"  {'SNP':<14} {'Gene':<{gene_w}} {'Geno':<6} {'Risk':<5} {'#':<4} "
        f"{'Pub OR':<9} {'95% CI':<16} {'Your OR':<10} {'RAF':<6}"
    )
    sep_w = 14 + gene_w + 6 + 5 + 4 + 9 + 16 + 10 + 6
    print("  " + "-" * sep_w)

    for snp, genotype, n_risk, geno_or in results:
        if n_risk < 0:
            n_str, yor_str = "--", "--"
        elif n_risk == 0:
            n_str, yor_str = "0", "1.000"
        else:
            n_str = str(n_risk)
            yor_str = f"{geno_or:.3f}"

        ci_str = f"({snp.ci_lower:.2f}-{snp.ci_upper:.2f})"
        shared_tag = ""
        if snp.rsid in shared_snps:
            others = [c for c in shared_snps[snp.rsid] if c != panel.short_name]
            if others:
                shared_tag = f"  [{','.join(others)}]"

        print(
            f"  {snp.rsid:<14} {snp.gene:<{gene_w}} {genotype:<6} {snp.risk_allele:<5} "
            f"{n_str:<4} {snp.odds_ratio:<9.2f} {ci_str:<16} {yor_str:<10} "
            f"{snp.risk_allele_freq:<6.2f}{shared_tag}"
        )

    # Compute polygenic score
    combined_or, ci_low, ci_high, log_or_score, combined_se = compute_polygenic_score(
        score_inputs
    )

    # Z-score against population distribution (for SNPs that were genotyped)
    genotyped_snps = [snp for snp, _ in score_inputs]
    pop_mean_geno, pop_sd_geno = compute_expected_score_stats(genotyped_snps)
    z_score = (log_or_score - pop_mean_geno) / pop_sd_geno if pop_sd_geno > 0 else 0.0
    percentile = score_to_percentile(z_score)

    # Z-score 95% CI
    log_or_low = log_or_score - 1.96 * combined_se
    log_or_high = log_or_score + 1.96 * combined_se
    z_low = (log_or_low - pop_mean_geno) / pop_sd_geno if pop_sd_geno > 0 else 0.0
    z_high = (log_or_high - pop_mean_geno) / pop_sd_geno if pop_sd_geno > 0 else 0.0
    pctile_low = score_to_percentile(z_low)
    pctile_high = score_to_percentile(z_high)

    # Interpretation (based on population Z)
    if z_score < -1.0:
        interp = "WELL BELOW population average"
    elif z_score < -0.5:
        interp = "BELOW population average"
    elif z_score < 0.5:
        interp = "AVERAGE genetic load"
    elif z_score < 1.0:
        interp = "ABOVE population average"
    elif z_score < 2.0:
        interp = "NOTABLY ABOVE population average"
    else:
        interp = "WELL ABOVE population average"

    pct_str = (
        f"{(total_risk_alleles / total_possible * 100):.1f}%"
        if total_possible > 0
        else "N/A"
    )

    # Summary section
    print()
    print(f"  --- QUANTIFIED SUMMARY: {panel.short_name} ---")
    print()
    print(f"    SNPs in panel:             {len(panel.snps)}")
    print(f"    SNPs genotyped:            {found_count}")
    print(f"    SNPs with risk allele(s):  {risk_count}")
    print(
        f"    Total risk alleles:        {total_risk_alleles} / {total_possible} possible"
    )
    print(f"    Risk allele proportion:    {pct_str}")
    print()
    print("  --- POLYGENIC SCORE (population Z) ---")
    print()
    print(f"    Combined OR:         {combined_or:.3f}")
    print(f"    95% CI:              {ci_low:.3f} - {ci_high:.3f}")
    print(f"    Log-OR score:        {log_or_score:.4f} (yours)")
    print(f"    Population expected: {pop_mean_geno:.4f} +/- {pop_sd_geno:.4f}")
    print(
        f"    Z-score:             {z_score:+.2f}  (95% CI: {z_low:+.2f} to {z_high:+.2f})"
    )
    print(
        f"    Percentile:          {percentile:.1f}th  (95% CI: {pctile_low:.1f}th to {pctile_high:.1f}th)"
    )
    print()
    print(f"    >> {interp}")

    # Gene interaction effects
    interactions = check_interactions(snp_data, panel.snps)
    if interactions:
        print()
        print("  --- GENE INTERACTIONS DETECTED ---")
        for rsid1, rsid2, gene1, gene2, inter_or, n1, n2, boost, desc in interactions:
            print()
            print(f"    {gene1} x {gene2}  (interaction OR: {inter_or:.2f})")
            print(f"      {rsid1} ({n1}x risk) + {rsid2} ({n2}x risk)")
            print(f"      Effective boost: {boost:.3f}x beyond multiplicative")
            print(f"      {desc}")

    # Flagged variants
    flagged = [(s, g, n, o) for s, g, n, o in results if n > 0]
    if flagged:
        print()
        print("  --- FLAGGED VARIANTS (>= 1 risk allele) ---")
        for snp, genotype, n_risk, geno_or in flagged:
            label = "ELEVATED" if n_risk == 2 else "MODERATE"
            print()
            print(f"    {snp.rsid} ({snp.gene}) -- {label}")
            print(f"      Genotype:  {genotype}  ({n_risk}x {snp.risk_allele})")
            print(
                f"      Your OR:   {geno_or:.3f}  (published: {snp.odds_ratio:.2f}, 95% CI: {snp.ci_lower:.2f}-{snp.ci_upper:.2f})"
            )
            print(f"      RAF:       {snp.risk_allele_freq:.0%}")
            print(f"      {snp.description}")
            print(f"      Ref: {snp.source}")

    # Notes
    print()
    print(f"  Notes for {panel.short_name}:")
    for i, note in enumerate(panel.notes, 1):
        print(f"    {i}. {note}")
    print()


def print_cross_condition_summary(
    panels: list[ConditionPanel],
    snp_data: dict[str, tuple[str, str, str]],
) -> None:
    """Print a summary comparing risk across all analyzed conditions."""

    # Collect all row data first to compute column widths
    row_data = []
    for panel in panels:
        score_inputs = []
        found = risk = total_risk = total_poss = 0
        for snp in panel.snps:
            if snp.rsid in snp_data:
                found += 1
                _, _, genotype = snp_data[snp.rsid]
                n = count_risk_alleles(genotype, snp.risk_allele)
                if n >= 0:
                    score_inputs.append((snp, n))
                if n > 0:
                    risk += 1
                    total_risk += n
                    total_poss += 2
                elif n == 0:
                    total_poss += 2

        pct = f"{(total_risk / total_poss * 100):.1f}%" if total_poss > 0 else "N/A"
        allele_str = f"{total_risk}/{total_poss}" if total_poss > 0 else "--"

        combined_or, ci_low, ci_high, log_or_score, comb_se = compute_polygenic_score(
            score_inputs
        )
        genotyped_snps = [s for s, _ in score_inputs]
        pop_mean_g, pop_sd_g = compute_expected_score_stats(genotyped_snps)
        z = (log_or_score - pop_mean_g) / pop_sd_g if pop_sd_g > 0 else 0.0
        ptile = score_to_percentile(z)

        ci_str = f"{ci_low:.2f}-{ci_high:.2f}"
        row_data.append(
            (
                panel.short_name,
                str(found),
                allele_str,
                pct,
                f"{combined_or:.3f}",
                ci_str,
                f"{z:+.2f}",
                f"{ptile:.1f}",
            )
        )

    # Column headers
    headers = (
        "Condition",
        "Geno",
        "Alleles",
        "%",
        "Comb.OR",
        "OR 95%CI",
        "Z",
        "%ile",
    )
    # Compute widths
    cols = list(zip(headers, *row_data))
    widths = [max(len(str(v)) for v in col) + 2 for col in cols]

    print()
    total_w = sum(widths) + 4
    print("=" * total_w)
    print("  CROSS-CONDITION SUMMARY")
    print("=" * total_w)
    print()

    def fmt_row(vals: tuple, ws: list[int]) -> str:
        return "  " + "".join(f"{str(v):<{w}}" for v, w in zip(vals, ws))

    print(fmt_row(headers, widths))
    print("  " + "-" * (sum(widths)))
    for rd in row_data:
        print(fmt_row(rd, widths))

    # Shared variants
    shared = find_shared_snps(panels)
    shared_in_data = {r: c for r, c in shared.items() if r in snp_data}
    if shared_in_data:
        print()
        print("  Cross-condition SNPs found in your data:")
        for rsid, conditions in sorted(shared_in_data.items()):
            _, _, genotype = snp_data[rsid]
            gene = "?"
            for p in panels:
                for s in p.snps:
                    if s.rsid == rsid:
                        gene = s.gene
                        break
            print(f"    {rsid} ({gene})")
            print(f"      Genotype: {genotype}  |  {', '.join(conditions)}")

    print()


# ---------------------------------------------------------------------------
# Composite Indices: "Living with Intensity" & "Rainforest Mind"
# ---------------------------------------------------------------------------

# Each composite index is a weighted blend of panel Z-scores.
# Weights reflect how central each panel is to the concept.


def _compute_panel_zscore(
    panel: ConditionPanel,
    snp_data: dict[str, tuple[str, str, str]],
) -> tuple[float, float, float]:
    """Compute population Z-score for a panel. Returns (z, z_lo, z_hi)."""
    score_inputs = []
    for snp in panel.snps:
        if snp.rsid in snp_data:
            _, _, genotype = snp_data[snp.rsid]
            n = count_risk_alleles(genotype, snp.risk_allele)
            if n >= 0:
                score_inputs.append((snp, n))

    if not score_inputs:
        return 0.0, 0.0, 0.0

    _, _, _, log_or_score, comb_se = compute_polygenic_score(score_inputs)
    genotyped = [s for s, _ in score_inputs]
    pop_mean, pop_sd = compute_expected_score_stats(genotyped)

    lor_lo = log_or_score - 1.96 * comb_se
    lor_hi = log_or_score + 1.96 * comb_se

    if pop_sd <= 0:
        return 0.0, 0.0, 0.0

    z = (log_or_score - pop_mean) / pop_sd
    z_lo = (lor_lo - pop_mean) / pop_sd
    z_hi = (lor_hi - pop_mean) / pop_sd

    return z, z_lo, z_hi


def print_composite_indices(snp_data: dict[str, tuple[str, str, str]]) -> None:
    """
    Print composite book-relevance indices using analytical rG-corrected Z-test.

    For each composite, computes a weighted sum of panel Z-scores and
    standardizes by the composite standard deviation, which accounts for
    cross-panel genetic correlations (rG):

        Z_composite = sum(w_i * Z_i) / sqrt(sum_ij w_i * w_j * rG(i,j))

    This is equivalent to a Welch's t-test comparing the user's multi-panel
    profile against the population null, with non-independent panels handled
    via the covariance structure.
    """

    print()
    print("=" * 100)
    print("  BOOK RELEVANCE INDEX")
    print("=" * 100)
    print()
    print("  For each book, we compute a weighted composite Z across its relevant")
    print("  trait panels, corrected for genetic correlations between panels.")
    print("  The percentile tells you: 'You are more likely to resonate with this")
    print("  book than X% of the general population, based on these markers.'")
    print()

    # Collect results for summary table
    book_results = []

    for index_name, config in COMPOSITE_INDICES.items():
        cross_rg = config.get("cross_rg", {})

        # Compute per-panel Z-scores
        panel_z_data = []
        for panel_key, weight, dimension in config["components"]:
            panel = ALL_PANELS.get(panel_key)
            if panel is None:
                continue
            z, _, _ = _compute_panel_zscore(panel, snp_data)
            panel_z_data.append((panel_key, dimension, weight, z))

        if not panel_z_data:
            continue

        # Weighted sum of Z-scores
        weighted_z_sum = sum(w * z for _, _, w, z in panel_z_data)

        # Composite variance: Var = sum_ij w_i * w_j * rG(i,j), rG(i,i)=1
        n = len(panel_z_data)
        composite_var = 0.0
        for i in range(n):
            ki, _, wi, _ = panel_z_data[i]
            for j in range(n):
                kj, _, wj, _ = panel_z_data[j]
                if i == j:
                    rg = 1.0
                else:
                    rg = cross_rg.get((ki, kj), cross_rg.get((kj, ki), 0.0))
                composite_var += wi * wj * rg

        composite_sd = math.sqrt(composite_var) if composite_var > 0 else 1.0
        composite_z = weighted_z_sum / composite_sd
        composite_pct = score_to_percentile(composite_z)

        book_results.append(
            (index_name, config, panel_z_data, composite_z, composite_pct)
        )

    # Sort by percentile descending (most relevant first)
    book_results.sort(key=lambda x: -x[4])

    # Print ranked results
    print(f"  {'#':<4} {'Book':<55} {'Z':<7} {'%ile':<7} {'Relevance':<16}")
    print(f"  {'-' * 89}")
    for rank, (name, config, _, cz, cpct) in enumerate(book_results, 1):
        if cpct >= 80:
            rel = "HIGHLY RELEVANT"
        elif cpct >= 60:
            rel = "RELEVANT"
        elif cpct >= 40:
            rel = "SOMEWHAT"
        else:
            rel = "LESS RELEVANT"
        author = config.get("author", "")
        label = f"{name} ({author})" if author else name
        print(f"  {rank:<4} {label:<60} {cz:<+7.2f} {cpct:<7.1f} {rel:<16}")

    # Detailed breakdown for top books
    print()
    print("  --- DETAILED BREAKDOWNS ---")

    for name, config, panel_z_data, composite_z, composite_pct in book_results:
        author = config.get("author", "")
        print()
        print(f"  {name}")
        if author:
            print(f"  by {author}")
        for line in config["description"].split("\n"):
            print(f"    {line}")
        print()
        print(f"    {'Dimension':<40} {'Panel':<12} {'Wt':<6} {'Z':<7} {'%ile':<6}")
        print(f"    {'-' * 71}")
        for pkey, dimension, weight, z in panel_z_data:
            pct = score_to_percentile(z)
            print(
                f"    {dimension:<40} {pkey.upper():<12} {weight:<6.2f} {z:<+7.2f} {pct:<6.1f}"
            )

        print()
        print(f"    Composite Z:    {composite_z:+.2f}")
        print(f"    Percentile:     {composite_pct:.1f}th")
        print("    -> You are more likely to resonate with this book than")
        print(f"       {composite_pct:.1f}% of the general population.")

        # Visual bar
        bar_len = 50
        pos = max(0, min(bar_len - 1, int(composite_pct / 100 * bar_len)))
        bar = list("." * bar_len)
        bar[pos] = "#"
        print(f"    0%  |{''.join(bar)}|  100%")
        print(f"    {' ' * (pos + 9)}^ {composite_pct:.1f}%")

    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    condition_keys = list(ALL_PANELS.keys())

    parser = argparse.ArgumentParser(
        description="Analyze 23andMe raw data for neurodiversity-associated SNPs.",
        epilog="DISCLAIMER: Educational use only. Not a diagnostic tool.",
    )
    parser.add_argument("datafile", help="Path to your 23andMe raw data text file.")
    parser.add_argument(
        "--condition",
        "-c",
        choices=condition_keys + ["all"],
        default="all",
        help=f"Condition to analyze. Choices: {', '.join(condition_keys)}, all (default: all).",
    )
    args = parser.parse_args()

    print(f"Loading 23andMe data from: {args.datafile}")
    snp_data = parse_23andme_file(args.datafile)
    print(f"Loaded {len(snp_data)} SNPs from file.")

    if args.condition == "all":
        panels = list(ALL_PANELS.values())
    else:
        panels = [ALL_PANELS[args.condition]]

    print()
    print("*" * 100)
    print("  NEURODIVERSITY SNP ANALYSIS REPORT (Quantified Edition)")
    print("*" * 100)
    print()
    print("  DISCLAIMER: This is NOT a diagnostic tool. These conditions are")
    print("  highly polygenic and multifactorial. No SNP panel can predict or")
    print("  diagnose any condition. Consult a genetic counselor for interpretation.")
    print()
    print("  Methodology:")
    print("    * Odds ratios (OR) and 95% CIs from published GWAS/meta-analyses")
    print(
        "    * Log-additive (multiplicative) polygenic model: Combined OR = Product(OR_i ^ n_i)"
    )
    print(
        "    * Population percentiles estimated via Hardy-Weinberg expected distribution"
    )
    print("    * Risk allele frequencies (RAF) from reference populations")
    print()
    print(f"  Panels to analyze: {', '.join(p.short_name for p in panels)}")

    shared = find_shared_snps(panels)
    for panel in panels:
        analyze_panel(panel, snp_data, shared)

    if len(panels) > 1:
        print_cross_condition_summary(panels, snp_data)

    # Composite indices (only when running all panels)
    if args.condition == "all":
        print_composite_indices(snp_data)

    print("=" * 100)
    print("  GENERAL DISCLAIMER")
    print("=" * 100)
    print("""
  1. Each SNP has a VERY SMALL individual effect size (OR ~1.05-1.3).
     Carrying risk alleles does NOT mean you have or will develop any condition.

  2. The combined OR from this panel captures only a TINY fraction of total genetic
     risk. True polygenic risk scores (PRS) use thousands to millions of SNPs and
     require population-matched LD reference panels for proper calibration.

  3. 95% confidence intervals reflect uncertainty in the published effect estimates,
     NOT prediction intervals for individual outcomes.

  4. Population percentiles compare your score to the EXPECTED distribution for these
     specific SNPs only — not a clinically validated PRS percentile.

  5. These conditions involve hundreds to thousands of variants, rare mutations,
     CNVs, and de novo events NOT captured by consumer genotyping chips.

  6. Environmental factors, epigenetics, and gene-environment interactions play
     major roles in all neurodevelopmental conditions.

  7. This analysis should NOT be used for self-diagnosis. Consult a board-certified
     genetic counselor or clinical geneticist for proper interpretation.
""")


if __name__ == "__main__":
    main()
