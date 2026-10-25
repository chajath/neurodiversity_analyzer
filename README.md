# Neurodiversity SNP Analyzer

Analyzes downloaded 23andMe raw data for genetic variants (SNPs) associated with neurodevelopmental traits and conditions, using published GWAS data with quantified odds ratios, confidence intervals, and population percentiles.

## Disclaimer

**This is NOT a diagnostic tool.** These conditions are highly polygenic and multifactorial. No SNP panel can predict or diagnose any of them. Always consult a genetic counselor for interpretation.

## Features

- **10 trait/condition panels**: 117 verified SNPs from published GWAS (up to 2026), plus 166 legacy entries flagged unverified and excluded from scoring by default
- **Quantified scoring**: per-SNP odds ratios, 95% CIs, population Z-scores and percentiles
- **Gene-gene interaction detection**: 10 candidate epistatic pairs (all currently flagged unverified)
- **34-book relevance index**: ranked reading recommendations based on your genetic profile
- **YAML-driven data**: all panels, interactions, and books are in editable YAML files
- Accepts both `.txt` and `.zip` 23andMe raw data files

## Panels

| Panel | Key | Verified SNPs (total) | Main sources |
|-------|-----|------|---------|
| Autism Spectrum Disorder | `asd` | 8 (30) | Grove 2019, Matoba 2020 |
| ADHD | `adhd` | 13 (40) | Demontis 2023 |
| Dyslexia | `dyslexia` | 11 (23) | Doust 2022 (42 loci) |
| OCD | `ocd` | 13 (29) | Strom 2025 (PGC-OCD, 30 loci) |
| Tourette Syndrome | `tourette` | 8 (15) | TSAICG 2026 preprint, Tsetsos 2024, Yu 2019 |
| Bipolar Disorder | `bipolar` | 17 (29) | O'Connell 2025 (PGC BD4, 298 loci) |
| Schizophrenia | `scz` | 14 (29) | Trubetskoy 2022 (PGC3) |
| Cognitive Ability / Giftedness | `gifted` | 11 (30) | van den Berg 2026 (N=455k) |
| Sensory Processing Sensitivity | `sps` | 10 (30) | Schwaba 2026 Big Five (neuroticism, openness) |
| Psychomotor Intensity | `psychomotor` | 12 (28) | Schwaba 2026 (extraversion), Wang 2022 (physical activity), Jones 2019 |

The SPS and Psychomotor panels include **non-pathological trait markers** from personality GWAS (Openness, Extraversion, Empathy, Sensation Seeking, Chronotype) — not just disorder markers.

## Data Verification

In October 2026 every entry was checked against PubMed (does the cited PMID match the study?), Ensembl (is the rsid on the stated chromosome, and is the risk allele a real plus-strand allele?) and the source papers' supplementary tables or the GWAS Catalog. Entries that failed are kept but marked:

```yaml
  unverified: true
  unverified_reason: labelled chr6 but rsid is on chr15; cited PMID is a different paper (...)
```

Unverified SNPs and interactions are **excluded from scoring by default**, and listed under each panel in the report. Pass `--include-unverified` to score them anyway. Most of the classic candidate-gene SNPs (DRD4, DAT1, COMT, SLC1A1, KIAA0319, ...) are flagged: either their citations don't check out, or they never replicated in the large GWAS.

## Requirements

- Python 3.10+
- PyYAML (`pip install pyyaml`)

## Getting Your 23andMe Data

1. Log in at [23andMe](https://you.23andme.com/)
2. Go to **Settings → 23andMe Data → Download Raw Data**
3. Download the `.zip` or `.txt` file

## Usage

```bash
# Analyze all panels
python neurodiversity_analyzer.py genome_data.zip

# Analyze a specific condition
python neurodiversity_analyzer.py genome_data.txt -c asd
python neurodiversity_analyzer.py genome_data.txt -c gifted

# Also score entries flagged as unverified
python neurodiversity_analyzer.py genome_data.zip --include-unverified
```

Available panels: `asd`, `adhd`, `dyslexia`, `ocd`, `tourette`, `bipolar`, `scz`, `gifted`, `sps`, `psychomotor`

## Output

The report includes:

### Per-Panel Analysis
- SNP table with your genotype, published OR, 95% CI, your genotype OR, and risk allele frequency
- Population Z-score with 95% CI and percentile
- Gene-gene interaction effects when multiple risk genes co-occur
- Flagged variants with literature references

### Cross-Condition Summary
- All 10 panels compared side-by-side with Z-scores and percentiles

### Book Relevance Index
- 34 neurodiversity books ranked by how likely you are to resonate with them
- Each book has a composite Z-score computed from weighted panel contributions
- Genetic correlations between panels (rG) are accounted for in the composite variance
- Detailed per-book breakdowns showing which dimensions drive relevance

## Project Structure

```
neurodiversity_analyzer.py   # Analysis engine (loads data from YAML)
data/
  panels.yaml                # 10 panels, 283 SNPs (117 verified) with OR/CI/RAF
  interactions.yaml           # 10 gene-gene interaction effects (unverified)
  books.yaml                  # 34 books with weighted components and cross-rG
```

## Adding Data

All data lives in `data/*.yaml` — no Python changes needed.

**Add a SNP** to a panel: edit `data/panels.yaml`, add an entry under the panel's `snps` list.

**Add a book**: edit `data/books.yaml`, add a new entry with `author`, `description`, `components` (panel/weight/dimension tuples), and `cross_rg` (genetic correlations between component panels).

**Add an interaction**: edit `data/interactions.yaml`, add an entry with two rsids, genes, interaction OR, and description.

## Methodology

- **Odds ratios and 95% CIs** from published GWAS and meta-analyses
- **Quantitative traits** (personality, intelligence, activity) use pseudo-ORs: exp(beta) per SD for personality/activity, exp(beta*pi/sqrt(3)) for intelligence
- **Alleles are plus-strand** (GRCh37/38), matching 23andMe raw data
- **Log-additive polygenic model**: Combined OR = Product(OR_i ^ n_risk_alleles_i)
- **Population Z-scores**: your log-OR score vs. Hardy-Weinberg expected mean and SD
- **Gene interactions**: epistatic effects beyond multiplicative model
- **Book composites**: rG-corrected weighted Z across relevant panels:
  Z_composite = sum(w_i * Z_i) / sqrt(sum_ij(w_i * w_j * rG_ij))

## Limitations

- Individual SNP effect sizes are very small (OR ~1.02–1.1 in modern GWAS)
- The most impactful variants (rare CNVs, de novo mutations) are **not detectable** from consumer chips
- Environmental, epigenetic, and developmental factors are not captured
- SNP coverage depends on 23andMe chip version
- Cross-panel genetic correlations are LDSC estimates (Grotzinger 2026, Doust 2022, Demontis 2023, Schwaba 2026); trait panels use proxy-trait correlations
- This tool captures a tiny fraction of total genetic variance for any condition
