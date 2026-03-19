# Neurodiversity SNP Analyzer

Analyzes downloaded 23andMe raw data for genetic variants (SNPs) associated with neurodevelopmental traits and conditions, using published GWAS data with quantified odds ratios, confidence intervals, and population percentiles.

## Disclaimer

**This is NOT a diagnostic tool.** These conditions are highly polygenic and multifactorial. No SNP panel can predict or diagnose any of them. Always consult a genetic counselor for interpretation.

## Features

- **10 trait/condition panels** with 182 SNPs from published GWAS
- **Quantified scoring**: per-SNP odds ratios, 95% CIs, population Z-scores and percentiles
- **Gene-gene interaction detection**: 10 known epistatic pairs checked
- **34-book relevance index**: ranked reading recommendations based on your genetic profile
- **YAML-driven data**: all panels, interactions, and books are in editable YAML files
- Accepts both `.txt` and `.zip` 23andMe raw data files

## Panels

| Panel | Key | SNPs | Sources |
|-------|-----|------|---------|
| Autism Spectrum Disorder | `asd` | 27 | Grove 2019, Matoba 2020, Wang 2009 |
| ADHD | `adhd` | 28 | Demontis 2023, candidate genes |
| Dyslexia | `dyslexia` | 12 | DCDC2, KIAA0319, DYX1C1 |
| OCD | `ocd` | 17 | Burton 2021, SLC1A1, GRIN2B |
| Tourette Syndrome | `tourette` | 8 | Yu 2019, HDC, TPH2 |
| Bipolar Disorder | `bipolar` | 15 | Mullins 2021 (PGC3-BD) |
| Schizophrenia | `scz` | 18 | Trubetskoy 2022 (PGC3) |
| Cognitive Ability / Giftedness | `gifted` | 19 | Savage 2018, Hill 2019, Lam 2022 |
| Sensory Processing Sensitivity | `sps` | 20 | Aron HSP + Nagel 2018, Warrier 2018, Lo 2017 |
| Psychomotor Intensity | `psychomotor` | 18 | Dopamine/NE + Lo 2017, Klimentidis 2018, Jones 2019 |

The SPS and Psychomotor panels include **non-pathological trait markers** from personality GWAS (Openness, Extraversion, Empathy, Sensation Seeking, Chronotype) — not just disorder markers.

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
  panels.yaml                # 10 panels, 182 SNPs with OR/CI/RAF
  interactions.yaml           # 10 gene-gene interaction effects
  books.yaml                  # 34 books with weighted components and cross-rG
```

## Adding Data

All data lives in `data/*.yaml` — no Python changes needed.

**Add a SNP** to a panel: edit `data/panels.yaml`, add an entry under the panel's `snps` list.

**Add a book**: edit `data/books.yaml`, add a new entry with `author`, `description`, `components` (panel/weight/dimension tuples), and `cross_rg` (genetic correlations between component panels).

**Add an interaction**: edit `data/interactions.yaml`, add an entry with two rsids, genes, interaction OR, and description.

## Methodology

- **Odds ratios and 95% CIs** from published GWAS and meta-analyses
- **Log-additive polygenic model**: Combined OR = Product(OR_i ^ n_risk_alleles_i)
- **Population Z-scores**: your log-OR score vs. Hardy-Weinberg expected mean and SD
- **Gene interactions**: epistatic effects beyond multiplicative model
- **Book composites**: rG-corrected weighted Z across relevant panels:
  Z_composite = sum(w_i * Z_i) / sqrt(sum_ij(w_i * w_j * rG_ij))

## Limitations

- Individual SNP effect sizes are very small (OR ~1.05–1.3)
- The most impactful variants (rare CNVs, de novo mutations) are **not detectable** from consumer chips
- Environmental, epigenetic, and developmental factors are not captured
- SNP coverage depends on 23andMe chip version
- Cross-panel genetic correlations are approximate estimates from published LDSC analyses
- This tool captures a tiny fraction of total genetic variance for any condition
